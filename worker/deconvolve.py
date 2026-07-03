"""This scripts executes lollipop by building the wrangeling the input data
and configurations into foram on the fly.

Its a wrapper around the lollipop deconvolute command.

requires the command line tools:
- lollipop
- xsv
- gawk

"""

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict
import yaml
import subprocess
import json
import pandas as pd
import logging

# Configure logging for the worker
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def devconvolve(
    mutation_counts_df: pd.DataFrame,
    mutation_variant_matrix_df: pd.DataFrame,
    bootstraps: int = 0,
    bandwidth: int = 30,
    regressor: str = "robust",
    regressor_params: dict = {"f_scale": 0.01},
    deconv_params: dict = {"min_tol": 1e-3},
) -> Dict:
    """
    This function runs lollipop on the input data and returns the deconvoluted data.
    """

    with TemporaryDirectory() as tmpdir:
        # make the temporary directory
        tmpdir = Path(tmpdir)
        # make the output directory
        output_dir = tmpdir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        # make the input directory
        input_dir = tmpdir / "input"
        input_dir.mkdir(parents=True, exist_ok=True)

        # name the data to the input directory
        mutation_counts = Path("mutation_counts_coverage.csv")
        mutations_variant_matrix = Path("mutation_variant_matrix.csv")

        # paths to the input files
        mutation_counts_fp = input_dir / mutation_counts.name
        mutations_variant_matrix_fp = input_dir / mutations_variant_matrix.name

        # Reset the index to columns if it's a MultiIndex
        if isinstance(mutation_counts_df.index, pd.MultiIndex):
            logger.debug("Detected MultiIndex in mutation_counts_df, resetting to columns")
            mutation_counts_df = mutation_counts_df.reset_index()


        # Save the dataframes to CSV files in the input directory
        pd.DataFrame.to_csv(
            mutation_counts_df,
            mutation_counts_fp,
            index=False,
            sep=",",
        )
        pd.DataFrame.to_csv(
            mutation_variant_matrix_df,
            mutations_variant_matrix_fp,
            index=False,
            sep=",",
        )

        ###########################
        # make the variants_config.yaml file yaml
        # as
        """ variants_pangolin:
        "KP.3": "KP.3"
        "KP.2": "KP.2"
        "LP.8": "LP.8"

        """
        # get the variants from the mutation_variant_matrix_df, i.e. the columns
        variants_config = input_dir / "variants_config.yaml"
        with open(variants_config, "w") as f:
            # Get column names excluding the first column (which is Mutation)
            variants = list(mutation_variant_matrix_df.columns)[1:]
            variants_dict = {variant: variant for variant in variants}
            yaml.dump({"variants_pangolin": variants_dict}, f)

        ############################
        #  make the deconv_bootstrap_cowwid.yaml config
        """bootstrap: 0

        kernel_params:
        bandwidth: 30

        regressor: robust
        regressor_params:
        f_scale: 0.01

        deconv_params:
        min_tol: 1e-3
        """

        deconv_config_fp = input_dir / "deconv_bootstrap_cowwid.yaml"

        with open(deconv_config_fp, "w") as f:
            deconv_config = {
                "bootstrap": bootstraps,
                "kernel_params": {"bandwidth": bandwidth},
                "regressor": regressor,
                "regressor_params": regressor_params,
                "deconv_params": deconv_params,
                "remove_deletions": True,
            }
            yaml.dump(deconv_config, f)

        ############################
        #  prepare the input data i.e. make the tallymut

        # Generate separate "pos"ition and variant "base" columns from the mutation data
        gawk_command = [
            "gawk",
            r'BEGIN{FS=","};NR==1{print "pos,base," $0};NR>1{match($1,/[ATGC]?([[:digit:]]+)([ATCG\-])/, ary); print ary[1] "," ary[2] "," $0 }',
            str(mutations_variant_matrix_fp),
        ]

        # Create output filename with descriptive suffix
        matrix_pos_base_file = output_dir / (
            mutations_variant_matrix_fp.stem + "_pos_base.csv"
        )

        try:
            with open(matrix_pos_base_file, "w") as f:
                subprocess.run(
                    gawk_command,
                    check=True,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            logger.debug(f"Successfully parsed mutation positions and bases: {matrix_pos_base_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running gawk command: {e}")
            if e.stderr:
                logger.error(e.stderr)
            exit(1)

        # # add the  above-update matrix to the mutation counts+coverage table
        # xsv join --left mutation  /tmp/mutation_counts_coverage.csv  Mutation /tmp/mutation_variant_matrix2.csv |
        #  xsv select samplingDate,count,coverage,frequency,mutation,pos,base,9-  |
        #  xsv fmt --out-delimiter '\t'  | sed '1s/samplingDate/date/;1s/coverage/cov/;1s/frequency/frac/' >  /tmp/tallymut.tsv

        # Create output filename with descriptive suffix
        tallymut_file = output_dir / (mutation_counts.stem + "_tallymut.tsv")

        # do head in the terminal
        head_command = ["head", "-n", "5", str(mutation_counts_fp)]
        try:
            subprocess.run(
                head_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running head command: {e}")
            if e.stderr:
                logger.error(e.stderr)
            exit(1)



        try:
            # First subprocess: xsv join
            join_command = [
                "xsv",
                "join",
                "--left",
                "mutation",
                str(mutation_counts_fp),
                "Mutation",
                str(matrix_pos_base_file),
            ]
            
            # Debug: Print more details about the files
            logger.debug(f"Running command: {' '.join(join_command)}")
            logger.debug(f"Checking if files exist:")
            logger.debug(f"  - mutation_counts_fp exists: {Path(mutation_counts_fp).exists()}")
            logger.debug(f"  - matrix_pos_base_file exists: {Path(matrix_pos_base_file).exists()}")
            
            # Debug: Print file contents
            logger.debug(f"Contents of mutation_counts_fp (first 5 lines):")
            try:
                with open(mutation_counts_fp, 'r') as f:
                    for i, line in enumerate(f):
                        if i < 5:
                            logger.debug(f"  {line.strip()}")
                        else:
                            break
            except Exception as e:
                logger.error(f"Error reading mutation_counts_fp: {e}")
                
            logger.debug(f"Contents of matrix_pos_base_file (first 5 lines):")
            try:
                with open(matrix_pos_base_file, 'r') as f:
                    for i, line in enumerate(f):
                        if i < 5:
                            logger.debug(f"  {line.strip()}")
                        else:
                            break
            except Exception as e:
                logger.error(f"Error reading matrix_pos_base_file: {e}")
            
            # Try basic xsv command first to verify xsv works
            try:
                logger.debug("Testing basic xsv command...")
                test_result = subprocess.run(
                    ["xsv", "count", str(mutation_counts_fp)],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                logger.debug(f"xsv count result: {test_result.stdout.strip()}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Basic xsv command failed: {e}")
                logger.error(f"stderr: {e.stderr}")
            
            # Now try the actual join command
            join_result = subprocess.run(
                join_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            # Second subprocess: xsv select
            select_command = [
                "xsv",
                "select",
                "samplingDate,count,coverage,frequency,mutation,pos,base,9-",
            ]
            select_result = subprocess.run(
                select_command,
                input=join_result.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            # Third subprocess: xsv fmt
            fmt_command = ["xsv", "fmt", "--out-delimiter", "\t"]
            fmt_result = subprocess.run(
                fmt_command,
                input=select_result.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            # Fourth subprocess: sed
            sed_command = [
                "sed",
                "1s/samplingDate/date/;1s/coverage/cov/;1s/frequency/frac/",
            ]
            with open(tallymut_file, "w") as f:
                subprocess.run(
                    sed_command,
                    input=fmt_result.stdout,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
            logger.debug(f"Successfully created tally mutation file: {tallymut_file}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running command: {e}")
            if e.stderr:
                logger.error(f"stderr: {e.stderr}")
            
            # Debug - print file existence and content
            logger.debug(f"\nDebug info:")
            for file_path in [mutation_counts_fp, matrix_pos_base_file]:
                if Path(file_path).exists():
                    logger.debug(f"File exists: {file_path}")
                    try:
                        # Print first few lines of the file
                        with open(file_path, 'r') as f:
                            content = [next(f).strip() for _ in range(3)]
                            logger.debug(f"Content of {file_path}:\n{content}")
                    except Exception as read_err:
                        logger.error(f"Error reading file: {read_err}")
                else:
                    logger.error(f"File does not exist: {file_path}")
                    
            exit(1)

        ##################### run deconvolute - run Lollipop
        # lollipop deconvolute --output /tmp/deconvolved.csv --out-json /tmp/deconvolved.json \
        #   -c /tmp/lolli_config.yaml \
        #   --deconv-config /home/dryak/project/LolliPop/presets/deconv_bootstrap_cowwid.yaml \
        #   --namefield mutation /tmp/tallymut.tsv

        output_json_fp = output_dir / "deconvolved.json"
        output_csv_fp = output_dir / "deconvolved.csv"

        run_lollipop_command = [
            "lollipop",
            "deconvolute",
            "--output",
            str(output_csv_fp),
            "--out-json",
            str(output_json_fp),
            "-c",
            str(variants_config),
            "--deconv-config",
            str(deconv_config_fp),
            "--namefield",
            "mutation",
            str(tallymut_file),
            "--seed",
            str(42),
        ]

        try:
            subprocess.run(
                run_lollipop_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            logger.info(f"Successfully deconvoluted: {output_csv_fp}")

            # Read and format the JSON file properly
            with open(output_json_fp, "r") as f:
                deconvolved_data = json.loads(f.read())

            # Write back with proper indentation
            with open(output_json_fp, "w") as f:
                json.dump(deconvolved_data, f, indent=4)

            logger.debug(f"Formatted JSON output: {output_json_fp}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running lollipop command: {e}")
            if e.stderr:
                logger.error(e.stderr)
            exit(1)

        ###################################
        # read in the deconvoluted.json file
        with open(output_json_fp, "r") as f:
            deconvoluted_data = json.load(f)

    return deconvoluted_data

