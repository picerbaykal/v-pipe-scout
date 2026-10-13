"""Coverage quality per sample (process.sample_qc, 2026-10-13)."""
import os
import sys

sys.path.insert(0, "/app_shared")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from process.sample_qc import PROBES, quality_text, summarize, well_covered  # noqa: E402


def test_probes_every_250_bases_with_a_non_reference_base():
    pos = [int(p[:-1]) for p in PROBES]
    assert len(PROBES) == 119 and pos[0] == 100 and pos[-1] == 29600
    assert all(b - a == 250 for a, b in zip(pos, pos[1:]))


def test_a_spiky_sample_is_low_coverage():
    # Lugano 4 Sep-like: a few huge amplicons, almost nothing elsewhere
    cov = {p: 0 for p in PROBES}
    for p in PROBES[:20]:
        cov[p] = 200_000
    cov.update({p: 3 for p in PROBES[20:60]})
    q = summarize(cov)
    assert q["median"] == 3 and abs(q["share_100"] - 20 / 119) < 1e-9
    assert not q["well_covered"] and well_covered(q) is False
    assert quality_text(q).startswith("low coverage: median depth 3")


def test_an_even_sample_is_well_covered():
    q = summarize({p: 12_000 for p in PROBES})
    assert q["well_covered"] and q["share_100"] == 1.0 and q["top5"] < 0.05


def test_no_quality_data_is_none():
    assert well_covered(None) is None and quality_text({}) == ""
