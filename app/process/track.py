"""Is a lineage big enough to track? (2026-11-03)

The cooc rule says whether a lineage is there (process.national: confirmed,
+ Add after weeks in a row). A trial deconvolution — the panel plus the found
lineages — says how much of the virus it is. + Add needs both: tracking a
lineage at 0.1 % adds a near-duplicate signature to the deconvolution and
makes the other shares less stable, without telling anything.

Checked against the cooc cells on 18 Sep – 18 Oct 2026 (6 cities, panel + 8
found lineages): Spearman 0.76 in the well-covered cities; lineages with coocs
on 1 sample per city sat at 0.1–0.3 % (SW.2.1, TM.1), the noise floor set by
lineages not there (LP.8, XEC, BA.3.2) was up to ~0.8 %; TA.1, PJ.2.1, ST.1
reached 8–43 %. Thin cities (fewer than check.thin_well_min well-covered
samples) gave shares without coocs (ST.1 Chur 12.7 %) and are left out.

Pure logic: no Streamlit, no IO.
"""
from datetime import date, timedelta
from typing import Dict, Iterable, Optional

TRACK_DEFAULTS = {
    "add_min_share": 0.02,   # + Add: at least this in one well-covered city
    "small_share": 0.01,     # under add_min_share and at least this: "small"
    "last_days": 14,         # the shares are the mean of the last days of the trial
}


def track_cfg(cfg: Optional[dict] = None) -> dict:
    out = dict(TRACK_DEFAULTS)
    try:
        from utils.config import get_cooc_setting
        for k, v in TRACK_DEFAULTS.items():
            out[k] = type(v)(get_cooc_setting(f"track.{k}", v))
    except Exception:
        pass
    out.update(cfg or {})
    return out


def _series(result: dict, city: str, variant: str):
    """[(date, proportion)] of one variant in one city's deconvolution result
    (the worker's {city: {variant: {timeseriesSummary: [...]}}} or unwrapped)."""
    d = result
    if isinstance(d, dict) and city in d and isinstance(d[city], dict):
        d = d[city]
    ts = ((d or {}).get(variant) or {}).get("timeseriesSummary", []) or []
    out = []
    for e in ts:
        try:
            out.append((date.fromisoformat(str(e.get("date"))[:10]), float(e.get("proportion") or 0)))
        except Exception:
            continue
    return sorted(out)


def trial_shares(results: Dict[str, dict], lineages: Iterable[str],
                 well_cities: Iterable[str], cfg: Optional[dict] = None) -> Dict[str, Dict[str, float]]:
    """{lineage: {city: mean share over the last `last_days` of that city's
    trial}} for the well-covered cities only."""
    c = track_cfg(cfg)
    well = set(well_cities)
    out: Dict[str, Dict[str, float]] = {}
    for v in lineages:
        for city, res in (results or {}).items():
            if city not in well:
                continue
            s = _series(res, city, v)
            if not s:
                continue
            last = s[-1][0]
            vals = [p for d, p in s if d > last - timedelta(days=int(c["last_days"]))]
            if vals:
                out.setdefault(v, {})[city] = sum(vals) / len(vals)
    return out


def trial_stats(results: Dict[str, dict], lineages: Iterable[str],
                well_cities: Iterable[str], cfg: Optional[dict] = None) -> Dict[str, Dict[str, dict]]:
    """{lineage: {city: {mean, last, last_date, peak, peak_date}}} over the
    last `last_days` (2026-11-04: the mean alone hid TA.1 in Zürich going
    from 41 % on 3 Oct to 1.3 % on 16 Oct)"""
    c = track_cfg(cfg)
    well = set(well_cities)
    out: Dict[str, Dict[str, dict]] = {}
    for v in lineages:
        for city, res in (results or {}).items():
            if city not in well:
                continue
            s = _series(res, city, v)
            if not s:
                continue
            last = s[-1][0]
            win = [(d, p) for d, p in s if d > last - timedelta(days=int(c["last_days"]))]
            if not win:
                continue
            pk = max(win, key=lambda x: x[1])
            out.setdefault(v, {})[city] = {
                "mean": sum(p for _d, p in win) / len(win), "last": win[-1][1],
                "last_date": win[-1][0].isoformat(), "peak": pk[1], "peak_date": pk[0].isoformat()}
    return out


def size(shares: Optional[Dict[str, float]], cfg: Optional[dict] = None):
    """("track" | "small" | "negligible", best city, best share), or None
    without shares. track: >= add_min_share in one city; small: >= small_share."""
    if not shares:
        return None
    c = track_cfg(cfg)
    city, best = max(shares.items(), key=lambda kv: kv[1])
    word = ("track" if best >= c["add_min_share"] else
            "small" if best >= c["small_share"] else "negligible")
    return word, city, best
