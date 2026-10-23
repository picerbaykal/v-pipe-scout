"""All cities pooled (process.national, 2026-10-22): the national answer and
one vote per week, and a sample's square counting its reads like the cell
(process.variant_explorer.check_in_data). Pure functions."""
import sys

sys.path.insert(0, "/app_shared")

from process.national import national              # noqa: E402
from process.variant_explorer import check_in_data  # noqa: E402

MK = ["m1", "m2", "m3", "m4"]


def _city(samples, well=()):
    """samples: {date: {marker: [cov, hit, link_n, link_ok]}} (already counted)."""
    q = {d: {"share_100": 0.9 if d in well else 0.1, "median": 100} for d in samples}
    return {"markers": MK, "samples": samples, "quality": q}


def _cells(cov, hit):
    return {m: [cov, hit, hit, hit] for m in MK}


def test_cities_are_pooled_for_the_answer():
    # each city alone has under min_cov (100) reads at every marker; together 120
    nat = national({"GE": _city({"2026-08-03": _cells(60, 30)}),
                    "ZH": _city({"2026-08-04": _cells(60, 30)})})
    assert nat["answer"]["verdict"] == "present" and nat["answer"]["n_present"] == 4
    assert nat["n_cities"] == 2


def test_one_vote_per_week_monday_first():
    nat = national({"GE": _city({"2026-08-03": _cells(500, 0),      # Mon
                                 "2026-08-09": _cells(500, 0),      # Sun, same week
                                 "2026-08-10": _cells(500, 250)},   # next Mon
                                well={"2026-08-10"})})
    w = nat["weeks"]
    assert [x["start"] for x in w] == ["2026-08-03", "2026-08-10"]
    assert w[0]["answer"]["verdict"] == "absent" and w[0]["n_samples"] == 2
    assert w[1]["answer"]["verdict"] == "present" and w[1]["n_well"] == 1
    assert w[1]["cities"] == {"GE": "present"}


def test_no_samples_no_result():
    assert national({"GE": {"markers": MK, "samples": {}}}) is None


def test_a_thin_sample_square_counts_like_the_cell():
    # median depth 2 → at most 20 reads per marker: a jackpot of 50,000 reads
    # at 30 % no longer makes a firm square (until 2026-10-22 it said present)
    pd = {"2026-08-03": {m: [50_000, 15_000, 15_000, 15_000] for m in MK}}
    q = {"2026-08-03": {"median": 2, "share_100": 0.05}}
    r = check_in_data(MK, pd, ["2026-08-03"], quality=q)
    assert r["timeline"] == [["2026-08-03", "uncovered"]]
    assert all(v[0] == 20 for v in r["samples"]["2026-08-03"].values())


def test_a_well_covered_sample_square_still_answers():
    pd = {"2026-08-03": {m: [50_000, 15_000, 15_000, 15_000] for m in MK}}
    q = {"2026-08-03": {"median": 5_000, "share_100": 0.9}}
    r = check_in_data(MK, pd, ["2026-08-03"], quality=q)
    assert r["timeline"] == [["2026-08-03", "present"]]


def test_a_not_seen_square_is_pale_not_amber():
    # amber means only "in panel, not seen in any city" (XFG Basel 27 Aug was amber)
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "components" / "variants_table.py"
    spec = importlib.util.spec_from_file_location("vt_pale", path)
    vt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vt)
    pd = {"2026-08-27": {m: [5_000, 0, 0, 0] for m in MK}}
    q = {"2026-08-27": {"median": 9_000, "share_100": 0.84}}
    d = check_in_data(MK, pd, ["2026-08-27"], quality=q)
    html = vt._detail_row("g", "XFG", {"Basel (BS)": d}, ["Basel (BS)"], "Basel (BS)", 9)
    assert vt.AMBER not in html and vt.D_PALE in html


def _week(start, verdict, n_well, cities):
    return {"start": start, "answer": {"verdict": verdict}, "n_well": n_well,
            "n_samples": 5, "cities": cities}


def test_sudden_appearance_after_two_trustworthy_quiet_weeks():
    # PJ.2.1, 5 Sep – 5 Oct: not seen on 4 and 5 well-covered samples, then
    # seen in Basel and Lugano in the same week
    from process.national import sudden
    nat = {"weeks": [_week("2026-09-14", "absent", 4, {}), _week("2026-09-21", "absent", 5, {}),
                     _week("2026-09-28", "present", 5, {"Basel (BS)": "present",
                                                        "Lugano (TI)": "present",
                                                        "Zürich (ZH)": "absent"})]}
    s = sudden(nat)
    assert s["cities"] == ["Basel (BS)", "Lugano (TI)"] and s["week"] == "2026-09-28"
    # not when a quiet week rests on thin samples, or only one city sees it,
    # or it was already seen before
    thin = {"weeks": [dict(nat["weeks"][0], n_well=1)] + nat["weeks"][1:]}
    assert sudden(thin) is None
    one = {"weeks": nat["weeks"][:2] + [_week("2026-09-28", "present", 5, {"Basel (BS)": "present"})]}
    assert sudden(one) is None
    seen = {"weeks": [_week("2026-09-14", "present", 4, {})] + nat["weeks"][1:]}
    assert sudden(seen) is None
