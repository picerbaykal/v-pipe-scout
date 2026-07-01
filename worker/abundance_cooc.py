"""
worker/abundance_cooc.py

LAPIS-sourced deconvolution worker for the Abundance & Co-occurrence tab.

Replaces the gawk/xsv tallymut-building chain in deconvolve.py with a
direct LAPIS fetch — our tallymut already has pos/base natively, so no
subprocess preprocessing is needed before calling lollipop deconvolute.

Pipeline:
    WiseLoculusLapis.get_tallymut()
    → build binary variant-membership columns from PangoLoader signatures
    → write tallymut.tsv + variants_config.yaml + deconv_config.yaml
    → subprocess lollipop deconvolute (no --namefield flag needed)
    → parse deconvolved.json
    → return dict matching devconvolve() output shape
"""

import asyncio
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# Add app/ to sys.path so worker container can import from api/ and utils/
sys.path.insert(0, str(Path(__file__).parent.parent / "app"))


def _build_variant_membership_columns(
        df_tally: pd.DataFrame,
        variants: List[str],
        pango_loader,
) -> pd.DataFrame:
    """
    Add binary variant-membership columns to the tallymut DataFrame.

    For each variant, adds a column where:
        1 = this mutation position is in the variant's signature
        0 = it is not

    This must be done manually before writing the tallymut TSV because
    the lollipop CLI reads whatever columns already exist in the file —
    it does NOT derive membership from variants_pangolin at read time.
    (Gotcha #5: first symptom is "Warning, variants_list's {'XFG'} is
    not present in columns".)
    """
    df = df_tally.copy()
    for variant in variants:
        signature = pango_loader.get_signature(variant)
        df[variant] = df["pos"].isin(signature).astype(int)
    return df

def run_deconv_lapis(
        location: str,
        start_date: datetime,
        end_date: datetime,
        variants: List[str],
        bootstraps: int = 100,
        bandwidth: int = 10,
        regressor: str = "robust",
        regressor_params: dict = None,
        deconv_params: dict = None,
) -> Dict:
    """
    Run LolliPop deconvolution sourced live from LAPIS.

    Args:
        location: Location name e.g. "Lugano (TI)"
        start_date: Start of date range
        end_date: End of date range
        variants: Pango lineage names to include in the panel.
            Must be >= 2 (gotcha #1: LolliPop's general_preprocess
            silently empties results with exactly 1 variant).
        bootstraps: Bootstrap iterations.
            0 or 1 = point estimate only, no confidence intervals.
            > 1 = resampling with confidence intervals.
            Default 100 matches UI's Standard preset.
            Note: v0.5.3 behavior — older versions may differ.
        bandwidth: Gaussian kernel bandwidth. Default 10 = Narrow.
        regressor: Regressor type, default "robust".
        regressor_params: Regressor parameters.
        deconv_params: Deconvolution parameters.

    Returns:
        dict: Parsed deconvolved.json — same shape as devconvolve().

    Raises:
        ValueError: if < 2 variants or bootstraps < 1.
        RuntimeError: if LAPIS fetch or lollipop subprocess fails.
    """

    if regressor_params is None:
        regressor_params = {"f_scale": 0.01}
    if deconv_params is None:
        deconv_params = {"min_tol": 1e-3}

    if len(variants) < 2:
        raise ValueError(
            f"LolliPop requires at least 2 variants — got {len(variants)}. "
            "Gotcha #1: general_preprocess silently empties results with 1 variant."
        )

    if bootstraps < 0:
        raise ValueError(
            f"bootstraps must be >= 0, got {bootstraps}. "
            "bootstrap=0 or 1 gives a point estimate (no confidence intervals). "
            "bootstrap > 1 enables resampling and produces confidence intervals."
        )

    # ── 1. Fetch tallymut from LAPIS ─────────────────────────────────────────

    from api.wiseloculus import WiseLoculusLapis
    from utils.config import get_wiseloculus_url
    from api.pango_loader import PangoLoader, get_pango_summary_path

    wise_url = get_wiseloculus_url()
    client = WiseLoculusLapis(wise_url)
    pango_loader = PangoLoader(get_pango_summary_path())

    logger.info(f"Fetching tallymut: {location} {start_date.date()} → {end_date.date()}")

    try:
        df_tally = asyncio.run(
            client.get_tallymut(
                locationName=location,
                date_range=(start_date, end_date),
                variants=variants,
                pango_loader=pango_loader,
            )
        )
    except Exception as e:
        raise RuntimeError(f"LAPIS tallymut fetch failed for {location}: {e}") from e

    if df_tally.empty:
        raise RuntimeError(
            f"No mutation data returned from LAPIS for {location} "
            f"{start_date.date()} → {end_date.date()}. "
            "Check that this location and date range have data."
        )

    logger.info(f"Tallymut: {len(df_tally)} rows, {df_tally['date'].nunique()} dates")

    # ── 2. Build binary variant-membership columns ────────────────────────────
    df_tally = _build_variant_membership_columns(df_tally, variants, pango_loader)

    # ── 3. Write files + run lollipop ────────────────────────────────────────
    with TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        output_dir = tmpdir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write tallymut TSV - already in the shape, no gawk/xsv needed
        tallymut_fp = tmpdir / "tallymut.tsv"
        df_tally.to_csv(tallymut_fp, sep="\t", index=False)

        # Write variants_config.yaml
        # NOTe: start_date goes here, NOT in deconv_config
        variants_config_fp = tmpdir / "variants_config.json"
        with open(variants_config_fp, "w") as f:
            yaml.dump({
                "variants_pangolin": {v: v for v in variants},
                "var_dates": {
                    v: [[start_date.strftime("%Y-%m-%d"), None]]
                    for v in variants
                },
            }, f)

        # Write deconv_config.yaml
        deconv_config_fp = tmpdir / "deconv_config.yaml"
        with open(deconv_config_fp, "w") as f:
            yaml.dump({
                "bootstrap": bootstraps,
                "kernel_params": {"bandwidth": bandwidth},
                "regressor": regressor,
                "regressor_params": regressor_params,
                "deconv_params": deconv_params,
            }, f)

            # Output files
            output_json_fp = output_dir / "deconvolved.json"
            output_csv_fp = output_dir / "deconvolved.csv"

            # Run lollipop deconvolute
            # No --namefield flag — our tallymut already has pos/base natively
            # so lollipop auto-builds the mutations name field (gotcha #4)
            run_command = [
                "lollipop", "deconvolute",
                "--output", str(output_csv_fp),
                "--out-json", str(output_json_fp),
                "-c", str(variants_config_fp),
                "--deconv-config", str(deconv_config_fp),
                "--seed", "42",
                str(tallymut_fp),
            ]

            try:
                result = subprocess.run(
                    run_command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                logger.info(f"lollipop deconvolute completed")
                if result.stderr:
                    logger.debug(f"lollipop stderr: {result.stderr}")
            except subprocess.CalledProcessError as e:
                logger.error(f"lollipop deconvolute failed: {e}")
                if e.stderr:
                    logger.error(f"stderr: {e.stderr}")
                raise RuntimeError(
                    f"lollipop deconvolute failed: {e.stderr}"
                ) from e

            # ── 4. Parse and return results ───────────────────────────────────────
            with open(output_json_fp, "r") as f:
                deconvolved_data = json.load(f)

        return deconvolved_data
