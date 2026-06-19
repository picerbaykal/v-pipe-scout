"""This scripts executes deconvolutin with no smoothing and then covvfit infer.

requires the command line tool:
- covvfit
"""

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict
import subprocess
import base64
import pandas as pd
import logging
from deconvolve import devconvolve

# Configure logging for the worker
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_covvfit_inference(
        location_data: Dict,
        matrix_df: pd.DataFrame,
        max_days: int = 180,
        horizon: int = 90,
) -> Dict:
    """
    Runs no-smoothing deconvolution per location, then covvfit inference on the
    pooled result, and returns the generated figures.

    Args:
        location_data: raw mutation counts per city, keyed by location name.
        matrix_df: the mutation-variant matrix, shared across all locations.
        max_days: number of past days to restrict the analysis to.
        horizon: number of future days to forecast.

    Returns:
        dict with:
            "figure_png": base64-encoded PNG of the fit figure (for display).
            "figure_pdf": base64-encoded PDF of the fit figure (for download).
    """
    # no smooth deconvolution per location
    combined_results = {}
    for location, counts_df in location_data.items():
        logger.info(f"Running no-smoothing deconvolution for {location}")
        deconv_result = devconvolve(
            mutation_counts_df=counts_df,
            mutation_variant_matrix_df=matrix_df,
            bandwidth=0.1,
            bootstraps=0
        )
        combined_results[location] = deconv_result["location"]

    # reshape into flat table to prepare for covvfit run
    rows = []
    variant_names = set()
    for location, result_data in combined_results.items():
        variants = result_data.get(location, result_data)

        for variant, data in variants.items():
            for point in data.get("timeseriesSummary", []):
                rows.append({
                    "location": location,
                    "variant": variant,
                    "date": point["date"],
                    "proportion": point["proportion"],
                })
            # "undetermined" stays in the data (covvfit folds it into the
            # "other" baseline) but is never passed as a -v variant, so we
            # don't request a fitness estimate for the leftover bucket.
            if variant != "undetermined":
                variant_names.add(variant)

    if not rows:
        raise ValueError("No data found in deconvolution results.")

    variant_names = sorted(variant_names)

    df = pd.DataFrame(rows)

    logger.info(
        f"Prepared covvfit input with {len(df)} rows across "
        f"{df['location'].nunique()} locations and {len(variant_names)} variants."
    )

    with TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        input_csv = tmpdir / "deconvolved_nosmooth.csv"
        df.to_csv(input_csv, index=False)

        output_dir = tmpdir / "output"

        covvfit_command = [
            "covvfit", "infer",
            "--input", str(input_csv),
            "--output", str(output_dir),
            "--separator", ",",
            "--max-days", str(max_days),
            "--horizon", str(horizon),
        ]

        for variant in variant_names:
            covvfit_command += ["-v", variant]

        logger.info(f"Running command: {' '.join(covvfit_command)}")

        try:
            subprocess.run(
                covvfit_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            logger.info("Successfully ran covvfit inference")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running covvfit command: {e}")
            if e.stderr:
                logger.error(e.stderr)
            raise
        figure_png_path = output_dir / "figure.png"
        figure_pdf_path = output_dir / "figure.pdf"

        if not figure_png_path.exists():
            raise FileNotFoundError(
                f"covvfit did not produce figure.png at {figure_png_path}"
            )

        with open(figure_png_path, "rb") as f:
            figure_png_b64 = base64.b64encode(f.read()).decode("utf-8")

        result = {"figure_png": figure_png_b64}

        if figure_pdf_path.exists():
            with open(figure_pdf_path, "rb") as f:
                result["figure_pdf"] = base64.b64encode(f.read()).decode("utf-8")

        # pairwise fitness table (relative fitness between each pair of variants)
        pairwise_csv_path = output_dir / "pairwise_fitnesses.csv"
        if pairwise_csv_path.exists():
            with open(pairwise_csv_path, "r") as f:
                result["pairwise_fitnesses_csv"] = f.read()

    return result











