"""Tests for covvfit_runner (2026-10-08): the reshaping of LolliPop's output
for covvfit, the pairwise table, and (when the covvfit command is installed,
as in the worker image) one real `covvfit infer` run on made-up data."""
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# worker/ (this file's parent folder) holds covvfit_runner.py
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from covvfit_runner import deconv_rows, pairwise_rows, run_covvfit_lapis


def _lollipop_json(dates, shares, key="location"):
    """LolliPop --out-json shape: {location: {variant: {timeseriesSummary}}}."""
    out = {}
    for v, ps in shares.items():
        out[v] = {"timeseriesSummary": [{"date": d, "proportion": p,
                                         "proportionLower": None, "proportionUpper": None}
                                        for d, p in zip(dates, ps)]}
    return {key: out}


def test_rows_leave_out_undetermined_and_missing():
    j = _lollipop_json(["2026-01-20", "2026-01-24"],
                       {"KP.2": [0.6, 0.5], "XFG": [0.3, float("nan")],
                        "undetermined": [0.1, 0.2]})
    rows = deconv_rows("Lugano (TI)", j)
    assert {(r["variant"], r["date"]) for r in rows} == {
        ("KP.2", "2026-01-20"), ("KP.2", "2026-01-24"), ("XFG", "2026-01-20")}
    assert all(r["location"] == "Lugano (TI)" for r in rows)


def test_rows_read_the_real_location_name_too():
    j = _lollipop_json(["2026-01-20"], {"KP.2": [0.6]}, key="Lugano (TI)")
    assert len(deconv_rows("Lugano (TI)", j)) == 1


def test_pairwise_flags_other_and_drops_self_rows():
    csv_text = ("Variant\tReference_Variant\tEstimate\tLower_CI\tUpper_CI\n"
                "other\tother\t0.0\t0.0\t0.0\n"
                "XFG\tother\t10.7\t-15.9\t37.4\n"
                "XFG\tKP.2\t0.37\t0.35\t0.40\n")
    rows = pairwise_rows(csv_text)
    assert [(r["variant"], r["reference"], r["involves_other"]) for r in rows] == [
        ("XFG", "other", True), ("XFG", "KP.2", False)]
    assert rows[1]["estimate"] == pytest.approx(0.37)


def _fake_deconvolve(loc, d0, d1, vs):
    if loc == "Broken":
        raise RuntimeError("LAPIS timeout")
    rng = np.random.default_rng(0)
    dates = [d.date().isoformat() for d in pd.date_range(d0, d1, freq="4D")]
    shares = {v: [] for v in vs + ["undetermined"]}
    for i, _ in enumerate(dates):
        t = i / len(dates)
        lg = np.array([2 - 4 * t, 0.0, -2 + 5 * t, -3.0]) + rng.normal(0, .2, 4)
        p = np.exp(lg) / np.exp(lg).sum()
        for v, x in zip(shares, p):
            shares[v].append(float(x))
    return _lollipop_json(dates, shares)


@pytest.mark.skipif(shutil.which("covvfit") is None, reason="covvfit not installed")
def test_covvfit_end_to_end_and_a_failed_location_is_skipped():
    r = run_covvfit_lapis(["Lugano (TI)", "Broken"], datetime(2026, 1, 20),
                          datetime(2026, 4, 5), ["KP.2", "NB.1.8.1", "XFG"], horizon=30,
                          colors={"XFG": "#d62728"}, deconvolve=_fake_deconvolve)
    assert r["locations"] == ["Lugano (TI)"] and "Broken" in r["skipped"]
    assert r["figure_png"] and r["predictions_csv"]
    assert (r["date_min"], r["date_max"]) == ("2026-01-20", "2026-04-05")
    xfg_kp2 = [x for x in r["pairwise"] if (x["variant"], x["reference"]) == ("XFG", "KP.2")]
    assert xfg_kp2 and xfg_kp2[0]["estimate"] > 0      # XFG grows against KP.2


def test_no_samples_anywhere_raises():
    with pytest.raises(RuntimeError, match="No deconvolved samples"):
        run_covvfit_lapis(["Broken"], datetime(2026, 1, 20), datetime(2026, 4, 5),
                          ["KP.2", "XFG"], deconvolve=_fake_deconvolve)


class _DictCache(dict):
    def set(self, k, v, ex=None):
        self[k] = v


def test_second_run_uses_the_cache():
    calls = []

    def counting(loc, d0, d1, vs):
        calls.append(loc)
        return _fake_deconvolve(loc, d0, d1, vs)
    cache = _DictCache()
    args = (["Lugano (TI)"], datetime(2026, 1, 20), datetime(2026, 4, 5), ["KP.2", "XFG"])
    for _ in range(2):
        try:
            run_covvfit_lapis(*args, deconvolve=counting, cache=cache)
        except RuntimeError:
            pass                      # covvfit itself may be missing here; the cache is what's tested
    assert calls == ["Lugano (TI)"] and len(cache) == 1
