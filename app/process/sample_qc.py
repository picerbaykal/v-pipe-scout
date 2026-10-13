"""Coverage quality of each wastewater sample (2026-10-13).

In low-virus periods (summer 2026) most samples cover only part of the
genome: a median depth near 0, with a few amplicons copied up to hundreds of
thousands of times. Each of those pieces comes from few virus molecules, so
its mutation shares are snapshots of those molecules, not the city's mix,
and they flip between ~0 % and ~100 % from one sample to the next (Lugano
2026-09-04: median depth 2, 724,267 reads at position 8350). The sample is
there; it just says much less than its read count suggests.

Coverage is read from LAPIS at PROBES, one position every 250 bases (about
one per amplicon): /component/nucleotideMutationsOverTime reports the reads
covering a position whatever base is asked for, so each probe asks for a
base other than the reference. Checked against the sample page for Lugano
2026-09-04: median 2 · 30 % >= 10x · 23 % >= 100x · 18 % >= 1000x (page: 2 ·
32 % · 22 % · 17 %).

A sample is "well covered" when its median depth over the probes is >=
quality.min_median_depth and >= quality.min_share_100x of the probes have
>= 100 reads (config)."""
from __future__ import annotations

import asyncio
import logging
import statistics
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# position + a base other than the reference (Wuhan-Hu-1), every 250 bases
PROBES: List[str] = [
    "100G", "350G", "600C", "850G", "1100A", "1350T", "1600A", "1850T", "2100C",
    "2350A", "2600C", "2850A", "3100A", "3350T", "3600T", "3850A", "4100A", "4350T",
    "4600A", "4850T", "5100G", "5350T", "5600A", "5850G", "6100A", "6350T", "6600T",
    "6850C", "7100C", "7350G", "7600G", "7850C", "8100A", "8350A", "8600C", "8850A",
    "9100A", "9350T", "9600A", "9850A", "10100A", "10350G", "10600A", "10850C",
    "11100A", "11350C", "11600T", "11850A", "12100G", "12350T", "12600C", "12850T",
    "13100C", "13350C", "13600A", "13850C", "14100A", "14350C", "14600A", "14850T",
    "15100C", "15350A", "15600G", "15850C", "16100T", "16350T", "16600A", "16850T",
    "17100C", "17350C", "17600A", "17850G", "18100T", "18350T", "18600T", "18850T",
    "19100T", "19350A", "19600T", "19850A", "20100A", "20350G", "20600A", "20850A",
    "21100C", "21350C", "21600C", "21850C", "22100C", "22350G", "22600T", "22850T",
    "23100A", "23350A", "23600A", "23850A", "24100A", "24350C", "24600A", "24850A",
    "25100G", "25350G", "25600G", "25850C", "26100A", "26350C", "26600G", "26850A",
    "27100A", "27350G", "27600G", "27850A", "28100T", "28350C", "28600A", "28850T",
    "29100C", "29350T", "29600T",
]

QUALITY_DEFAULTS = {"min_median_depth": 20, "min_share_100x": 0.4}


def quality_cfg(cfg: Optional[dict] = None) -> dict:
    out = dict(QUALITY_DEFAULTS)
    try:
        from utils.config import get_cooc_setting
        for k, v in QUALITY_DEFAULTS.items():
            out[k] = type(v)(get_cooc_setting(f"quality.{k}", v))
    except Exception:
        pass
    out.update(cfg or {})
    return out


def summarize(cov: Dict[str, int], cfg: Optional[dict] = None) -> dict:
    """One sample: {median, share_10, share_100, share_1000, top5, max,
    well_covered} from {probe: reads covering it} (missing probe = 0)."""
    c = quality_cfg(cfg)
    xs = sorted(int(cov.get(p, 0) or 0) for p in PROBES)
    n, tot = len(xs), sum(xs)
    q = {"median": float(statistics.median(xs)) if xs else 0.0,
         "share_10": sum(x >= 10 for x in xs) / n if n else 0.0,
         "share_100": sum(x >= 100 for x in xs) / n if n else 0.0,
         "share_1000": sum(x >= 1000 for x in xs) / n if n else 0.0,
         "top5": sum(xs[-5:]) / tot if tot else 0.0,
         "max": xs[-1] if xs else 0}
    q["well_covered"] = bool(q["median"] >= c["min_median_depth"]
                             and q["share_100"] >= c["min_share_100x"])
    return q


def well_covered(q: Optional[dict], cfg: Optional[dict] = None) -> Optional[bool]:
    """The thresholds applied now (the UI may retune them without a re-run);
    None when the sample has no quality data."""
    if not q:
        return None
    c = quality_cfg(cfg)
    return bool(q.get("median", 0) >= c["min_median_depth"]
                and q.get("share_100", 0) >= c["min_share_100x"])


def quality_text(q: Optional[dict], cfg: Optional[dict] = None) -> str:
    """'low coverage: median depth 2 · 23 % of the genome ≥ 100 reads · 65 % of
    the reads in 5 places' (or 'well covered: …')."""
    ok = well_covered(q, cfg)
    if ok is None:
        return ""
    head = "well covered" if ok else "low coverage"
    return (f"{head}: median depth {q.get('median', 0):,.0f} · "
            f"{100 * q.get('share_100', 0):.0f} % of the genome ≥ 100 reads · "
            f"{100 * q.get('top5', 0):.0f} % of the reads in 5 places")


async def fetch_sample_quality(client, location: str, date_range, block: int = 4,
                               waits=(0, 10, 30, 60), cfg: Optional[dict] = None
                               ) -> Dict[str, dict]:
    """{date: summarize(...)} for every sampling date of `location` in
    `date_range`. Small requests (`block` dates each) with retries on
    gateway errors: LAPIS answered 504 to 15 dates at once."""
    import aiohttp
    dates = await client._get_sampling_dates(location, date_range)
    cov: Dict[str, Dict[str, int]] = {}
    async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=300)) as s:
        for i in range(0, len(dates), block):
            part = dates[i:i + block]
            res = None
            for w in waits:
                if w:
                    await asyncio.sleep(w)
                try:
                    res = await client._fetch_tallymut_block(s, location, PROBES, part)
                    break
                except Exception as e:
                    if not any(code in str(e) for code in ("502", "503", "504", "imeout")):
                        raise
                    logger.warning(f"[quality][{location}] retry {part[0]}…: {e}")
            if res is None:
                logger.warning(f"[quality][{location}] gave up on {part[0]}…{part[-1]}")
                continue
            content = (res or {}).get("data", {}) or {}
            names = content.get("mutations", []) or []
            ranges = content.get("dateRanges", []) or []
            for i_m, row in enumerate(content.get("data", []) or []):
                name = names[i_m] if i_m < len(names) else ""
                key = name[1:] if name[:1].isalpha() and name[1:2].isdigit() else name
                for j, cell in enumerate(row or []):
                    if j >= len(ranges):
                        break
                    d = str(ranges[j].get("dateFrom", ""))[:10]
                    cov.setdefault(d, {})[key] = int((cell or {}).get("coverage", 0) or 0)
    return {str(d)[:10]: summarize(cov.get(str(d)[:10], {}), cfg) for d in dates}


def sample_quality(client, location: str, date_range, **kw) -> Dict[str, dict]:
    """Blocking wrapper (the worker runs inside Celery, sometimes with a loop)."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(fetch_sample_quality(client, location, date_range, **kw))
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(1) as ex:
        return ex.submit(asyncio.run,
                         fetch_sample_quality(client, location, date_range, **kw)).result()
