"""Tests for the process.cooc module."""

import pandas as pd
import pytest

from process.cooc import (
    row_to_confirmed_sets,
    classify_pattern,
    annotate_cooc_dataframe,
    panel_completeness_by_date,
)


class TestRowToConfirmedSets:
    """Test the row_to_confirmed_sets function."""

    def test_all_reference(self):
        """A row with all reference bases → all tracked alts in absent."""
        row = {"[241]": "C", "[297]": "T"}
        amp_dict = {241: ["T"], 297: ["G"]}
        cp, ca = row_to_confirmed_sets(row, [241, 297], amp_dict)
        assert cp == set()
        assert ca == {"241T", "297G"}

    def test_mutation_at_all_positions(self):
        """A row with tracked mutations at all positions → all in present."""
        row = {"[241]": "T", "[297]": "G"}
        amp_dict = {241: ["T"], 297: ["G"]}
        cp, ca = row_to_confirmed_sets(row, [241, 297], amp_dict)
        assert cp == {"241T", "297G"}
        assert ca == set()

    def test_uncovered_position_excluded(self):
        """N at a position → contributes to neither set."""
        row = {"[241]": "T", "[297]": "N"}
        amp_dict = {241: ["T"], 297: ["G"]}
        cp, ca = row_to_confirmed_sets(row, [241, 297], amp_dict)
        assert cp == {"241T"}
        assert ca == set()

    def test_deletion_excluded_by_default(self):
        """'-' treated as uncovered when include_deletion_states=False."""
        row = {"[241]": "T", "[297]": "-"}
        amp_dict = {241: ["T"], 297: ["G"]}
        cp, ca = row_to_confirmed_sets(row, [241, 297], amp_dict)
        assert cp == {"241T"}
        assert ca == set()

    def test_deletion_included_when_configured(self):
        """When include_deletion_states=True, '-' is a real observation."""
        row = {"[241]": "T", "[297]": "-"}
        amp_dict = {241: ["T"], 297: ["G"]}
        cp, ca = row_to_confirmed_sets(
            row, [241, 297], amp_dict, include_deletion_states=True
        )
        # '-' is not in tracked_alts for 297, so it's treated as reference → absent
        assert cp == {"241T"}
        assert ca == {"297G"}

    def test_allele_conflation(self):
        """Multiple tracked alts at same position → all added when any observed."""
        row = {"[241]": "A"}  # A is one of the tracked alts
        amp_dict = {241: ["A", "T"]}
        cp, ca = row_to_confirmed_sets(row, [241], amp_dict)
        # Since A is tracked, both 241A and 241T go to confirmed_present
        assert cp == {"241A", "241T"}
        assert ca == set()

    def test_position_not_in_amp_dict(self):
        """Positions with no tracked mutations contribute nothing."""
        row = {"[500]": "C"}
        amp_dict = {}  # 500 has no tracked mutations
        cp, ca = row_to_confirmed_sets(row, [500], amp_dict)
        assert cp == set()
        assert ca == set()


class TestClassifyPattern:
    """Test the classify_pattern function."""

    def test_uninformative_empty(self):
        """Empty confirmed_present → uninformative."""
        assert classify_pattern(set(), {"V": {"241T"}}) == "uninformative"

    def test_uninformative_singleton(self):
        """Single mutation → uninformative (filtered as noise)."""
        assert classify_pattern({"241T"}, {"V": {"241T", "297G"}}) == "uninformative"

    def test_matched_by_one_variant(self):
        """Pattern's mutations all in one variant's signature → matched."""
        sigs = {"KP.2": {"241T", "297G", "500C"}}
        assert classify_pattern({"241T", "297G"}, sigs) == "matched"

    def test_matched_by_multiple_variants(self):
        """Pattern explained by more than one variant → still matched."""
        sigs = {
            "KP.2": {"241T", "297G"},
            "NB.1.8.1": {"241T", "297G", "500C"},
        }
        assert classify_pattern({"241T", "297G"}, sigs) == "matched"

    def test_unexplained(self):
        """No variant covers all mutations → unexplained."""
        sigs = {
            "KP.2": {"241T", "297G"},
            "NB.1.8.1": {"500C", "700A"},
        }
        # 241T is in KP.2, 500C is in NB.1.8.1, but no single variant has both
        assert classify_pattern({"241T", "500C"}, sigs) == "unexplained"


class TestAnnotateCoocDataframe:
    """Test the annotate_cooc_dataframe function."""

    def test_empty_dataframe(self):
        """Empty input → empty output."""
        result = annotate_cooc_dataframe(
            pd.DataFrame(), [241], {241: ["T"]}, {"V": {"241T"}}
        )
        assert result.empty

    def test_basic_annotation(self):
        """Rows get classified correctly end-to-end."""
        df = pd.DataFrame([
            {"date": "2025-11-09", "count": 100, "[241]": "T", "[297]": "G"},  # matched
            {"date": "2025-11-09", "count":  10, "[241]": "T", "[297]": "T"},  # singleton
            {"date": "2025-11-09", "count":  50, "[241]": "C", "[297]": "T"},  # reference
        ])
        amp_dict = {241: ["T"], 297: ["G"]}
        sigs = {"KP.2": {"241T", "297G"}}
        result = annotate_cooc_dataframe(df, [241, 297], amp_dict, sigs)
        assert len(result) == 3
        classifications = result["classification"].tolist()
        assert "matched" in classifications
        assert "uninformative" in classifications


class TestPanelCompletenessByDate:
    """Test the panel_completeness_by_date function."""

    def test_empty(self):
        result = panel_completeness_by_date(pd.DataFrame())
        assert result.empty

    def test_full_completeness(self):
        """All patterns matched → completeness = 1.0."""
        annotated = pd.DataFrame([
            {"date": "2025-11-09", "count": 100, "classification": "matched",
             "confirmed_present": [], "confirmed_absent": []},
            {"date": "2025-11-09", "count":  50, "classification": "matched",
             "confirmed_present": [], "confirmed_absent": []},
        ])
        result = panel_completeness_by_date(annotated)
        assert result.iloc[0]["completeness"] == 1.0

    def test_partial_completeness(self):
        """Mixed matched and unexplained → fractional completeness."""
        annotated = pd.DataFrame([
            {"date": "2025-11-09", "count": 75, "classification": "matched",
             "confirmed_present": [], "confirmed_absent": []},
            {"date": "2025-11-09", "count": 25, "classification": "unexplained",
             "confirmed_present": [], "confirmed_absent": []},
        ])
        result = panel_completeness_by_date(annotated)
        assert result.iloc[0]["completeness"] == 0.75

    def test_uninformative_ignored(self):
        """Uninformative patterns don't affect completeness."""
        annotated = pd.DataFrame([
            {"date": "2025-11-09", "count": 100, "classification": "matched",
             "confirmed_present": [], "confirmed_absent": []},
            {"date": "2025-11-09", "count": 9999, "classification": "uninformative",
             "confirmed_present": [], "confirmed_absent": []},
        ])
        result = panel_completeness_by_date(annotated)
        assert result.iloc[0]["completeness"] == 1.0