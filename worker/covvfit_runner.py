"""CovvFit for the Abundance & Co-occurrence page (Release 1, step 3b).

Relative growth advantages of the panel variants, from the run's locations,
panel and dates:

1. LolliPop per location WITHOUT smoothing (bandwidth 0.1 day, no bootstraps):
   covvfit needs each sample's own estimate, not the smoothed curve.
2. One table for covvfit: location, date, variant, proportion. LolliPop's
   "undetermined" stays out; covvfit adds its own "other" = 1 - sum of the
   listed variants, which is that same share.
3. `covvfit infer` from the run's start date to its end date, plus `horizon`
   days of prediction, with the page's variant colours.

Named covvfit_runner.py (not covvfit.py, as in PR #227) so it can never shadow
the covvfit package itself.

covvfit is installed from GitHub at a fixed commit (see the worker Dockerfile).
"""
from __future__ import annotations

import base64
import csv
import io
import json
import logging
import subprocess
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, Dict, List, Optional

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# LolliPop without smoothing: a Gaussian kernel this narrow (in days) gives the
# other sampling dates ~0 weight, so each date is estimated on its own.
NO_SMOOTHING_BANDWIDTH = 0.1
# covvfit's own name for 1 - (sum of the listed variants)
OTHER = "other"


def deconv_rows(location: str, deconvolved: Dict) -> List[dict]:
    """Rows {location, date, variant, proportion} from one location's LolliPop
    JSON. LolliPop names the location "location" when the input has no
    location column (ours doesn't), so that key is read, with the real name
    as a fallback. "undetermined" and missing values are left out."""
    per_variant = deconvolved.get("location") or deconvolved.get(location) or {}
    rows = []
    for variant, data in per_variant.items():
        if variant == "undetermined":
            continue
        for point in (data or {}).get("timeseriesSummary", []) or []:
            p = point.get("proportion")
            if p is None or p != p:          # None or NaN
                continue
            rows.append({"location": location, "variant": variant,
                         "date": str(point.get("date", ""))[:10],
                         "proportion": float(p)})
    return rows


def pairwise_rows(csv_text: str) -> List[dict]:
    """pairwise_fitnesses.csv -> [{variant, reference, estimate, lower, upper,
    involves_other}], without the trivial other-vs-other row. Rows with
    "other" are flagged: when the panel explains almost everything, "other"
    is ~0 and its advantage can't be estimated (very wide ranges)."""
    out = []
    for r in csv.DictReader(io.StringIO(csv_text), delimiter="\t"):
        v, ref = r.get("Variant", ""), r.get("Reference_Variant", "")
        if v == ref:
            continue
        try:
            est, lo, hi = (float(r["Estimate"]), float(r["Lower_CI"]), float(r["Upper_CI"]))
        except (KeyError, TypeError, ValueError):
            continue
        out.append({"variant": v, "reference": ref, "estimate": est, "lower": lo,
                    "upper": hi, "involves_other": OTHER in (v, ref)})
    return out


# A location's deconvolution without smoothing is kept this long (seconds) in
# the cache, so re-running CovvFit with other settings or locations doesn't
# fetch it from LAPIS again (2026-10-09: ~40 s per location over 180 days).
CACHE_SECONDS = 3600


def cache_key(location: str, start_date: datetime, end_date: datetime,
              variants: List[str]) -> str:
    return (f"covvfit:deconv:{location}:{start_date.date()}:{end_date.date()}:"
            f"{','.join(sorted(variants))}:{NO_SMOOTHING_BANDWIDTH}")


def run_covvfit_lapis(locations: List[str], start_date: datetime, end_date: datetime,
                      variants: List[str], horizon: int = 60,
                      colors: Optional[Dict[str, str]] = None,
                      progress: Optional[Callable] = None,
                      deconvolve: Optional[Callable] = None,
                      cache=None) -> Dict:
    """CovvFit on the run's locations, panel and dates. Returns
    {figure_png, figure_pdf (base64), pairwise: [...], pairwise_csv,
    predictions_csv, n_rows, locations, variants, date_min, date_max, horizon,
    skipped: {location: reason}}.

    deconvolve(location, start, end, variants) -> LolliPop JSON; defaults to
    abundance_cooc.run_deconv_lapis without smoothing (a parameter so tests
    can pass their own). cache: an object with get(key) and set(key, value,
    ex=seconds), e.g. the worker's redis client, or None."""
    if deconvolve is None:
        from abundance_cooc import run_deconv_lapis

        def deconvolve(loc, d0, d1, vs):
            return run_deconv_lapis(location=loc, start_date=d0, end_date=d1, variants=vs,
                                    bootstraps=0, bandwidth=NO_SMOOTHING_BANDWIDTH)
    n = len(locations)

    def _p(step, msg, frac=None):
        if progress:
            progress(step, msg, frac)

    # 1-2. per-location deconvolution without smoothing -> one table
    rows, skipped, from_cache = [], {}, []
    for i, loc in enumerate(locations):
        _p(1, f"Deconvolution without smoothing: {loc}...", i / max(n, 1))
        key = cache_key(loc, start_date, end_date, variants)
        try:
            raw = None
            if cache is not None:
                try:
                    raw = cache.get(key)
                except Exception:
                    raw = None
            if raw:
                deconvolved = json.loads(raw)
                from_cache.append(loc)
            else:
                deconvolved = deconvolve(loc, start_date, end_date, variants)
                if cache is not None:
                    try:
                        cache.set(key, json.dumps(deconvolved), ex=CACHE_SECONDS)
                    except Exception as e:
                        logger.warning(f"[covvfit] could not cache {loc}: {e}")
            got = deconv_rows(loc, deconvolved)
        except Exception as e:                       # one city must not sink the others
            logger.warning(f"[covvfit] deconvolution failed in {loc}: {e}")
            skipped[loc] = f"deconvolution failed: {e}"[:300]
            continue
        if not got:
            skipped[loc] = "no deconvolved samples"
            continue
        rows += got
    if not rows:
        raise RuntimeError("No deconvolved samples in any location"
                           + (f" ({'; '.join(f'{k}: {v}' for k, v in skipped.items())})"
                              if skipped else ""))
    df = pd.DataFrame(rows)
    used_vars = [v for v in variants if v in set(df["variant"])]
    used_locs = [l for l in locations if l in set(df["location"])]
    logger.info(f"[covvfit] input: {len(df)} rows, {len(used_locs)} location(s), "
                f"{len(used_vars)} variant(s); {len(from_cache)} from the cache")

    # 3. covvfit infer
    _p(2, f"covvfit: fitting {len(used_vars)} variants in {len(used_locs)} location(s)...")
    d_min, d_max = start_date.date().isoformat(), end_date.date().isoformat()
    with TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inp, out = tmp / "deconvolved_nosmooth.tsv", tmp / "output"
        df.to_csv(inp, sep="\t", index=False)
        cmd = ["covvfit", "infer", "--input", str(inp), "--output", str(out),
               "--separator", "\t", "--date-min", d_min, "--date-max", d_max,
               "--horizon", str(int(horizon))]
        if colors:
            cfg = tmp / "covvfit_config.yaml"
            with open(cfg, "w") as f:
                yaml.safe_dump({"plot": {"variant_colors": {v: c for v, c in colors.items()
                                                            if v in used_vars}}}, f)
            cmd += ["--config", str(cfg)]
        for v in used_vars:
            cmd += ["-v", v]
        logger.info(f"[covvfit] {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True)
        except subprocess.CalledProcessError as e:
            tail = (e.stderr or e.stdout or "").strip().splitlines()[-5:]
            logger.error(f"[covvfit] failed: {e}\n" + "\n".join(tail))
            raise RuntimeError("covvfit failed: " + (" | ".join(tail) or str(e))) from e

        _p(3, "Reading covvfit results...")
        png = out / "figure.png"
        if not png.exists():
            raise RuntimeError(f"covvfit did not write figure.png in {out}")
        result = {"figure_png": base64.b64encode(png.read_bytes()).decode("ascii"),
                  "n_rows": int(len(df)), "locations": used_locs, "variants": used_vars,
                  "date_min": d_min, "date_max": d_max, "horizon": int(horizon),
                  "skipped": skipped, "from_cache": from_cache}
        pdf = out / "figure.pdf"
        if pdf.exists():
            result["figure_pdf"] = base64.b64encode(pdf.read_bytes()).decode("ascii")
        pw = out / "pairwise_fitnesses.csv"
        if pw.exists():
            result["pairwise_csv"] = pw.read_text()
            result["pairwise"] = pairwise_rows(result["pairwise_csv"])
        pred = out / "predictions.csv"
        if pred.exists():
            result["predictions_csv"] = pred.read_text()
    return result
