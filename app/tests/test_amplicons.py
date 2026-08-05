"""Tests for the process.amplicons module."""

import pytest
from pathlib import Path

from process.amplicons import (
    load_amplicons,
    group_positions_by_amplicon,
    split_positions_by_distance,
    build_amp_dict_from_variants,
)


class TestLoadAmplicons:
    """Test the load_amplicons function."""

    def test_loads_real_bed_file(self):
        """The bundled SARS-CoV-2 BED file loads and contains ~99 amplicons."""
        bed_path = Path(__file__).parent.parent / "data" / "SARS-CoV-2.insert.bed"
        amplicons = load_amplicons(bed_path)
        assert len(amplicons) > 50
        # Each entry is (name, start, end) with valid ranges
        for name, start, end in amplicons:
            assert isinstance(name, str) and name
            assert isinstance(start, int) and start >= 1
            assert isinstance(end, int) and end > start


class TestGroupPositionsByAmplicon:
    """Test the group_positions_by_amplicon function."""

    def test_simple_grouping(self):
        """Positions get assigned to their containing amplicon."""
        amplicons = [
            ("amp1", 100, 400),
            ("amp2", 300, 600),
        ]
        positions = {150, 350, 500}
        grouped = group_positions_by_amplicon(positions, amplicons)
        # 150 → amp1 only; 350 → both (overlap); 500 → amp2 only
        assert 150 in grouped["amp1"]
        assert 350 in grouped["amp1"]
        assert 350 in grouped["amp2"]
        assert 500 in grouped["amp2"]

    def test_uncovered_positions_excluded(self):
        """Positions outside all amplicons don't appear in the result."""
        amplicons = [("amp1", 100, 200)]
        positions = {150, 500}  # 500 is outside
        grouped = group_positions_by_amplicon(positions, amplicons)
        assert 150 in grouped["amp1"]
        assert 500 not in grouped.get("amp1", [])

    def test_empty_amplicons_dropped(self):
        """Amplicons with no matching positions are not in the output."""
        amplicons = [
            ("amp1", 100, 200),
            ("amp2", 300, 400),
        ]
        positions = {150}  # only in amp1
        grouped = group_positions_by_amplicon(positions, amplicons)
        assert "amp1" in grouped
        assert "amp2" not in grouped


class TestSplitPositionsByDistance:
    """Test the split_positions_by_distance function."""

    def test_within_distance_stays_together(self):
        """Positions within max_distance form one batch."""
        amp_groups = {"amp1": [60, 100, 150, 200]}
        batches = split_positions_by_distance(amp_groups, max_distance=200)
        assert len(batches) == 1
        assert batches[0] == ("amp1", [60, 100, 150, 200])

    def test_span_exceeded_splits(self):
        """Positions beyond max_distance from group start create a new batch."""
        amp_groups = {"amp1": [60, 100, 260, 300]}
        batches = split_positions_by_distance(amp_groups, max_distance=200)
        assert len(batches) == 2
        assert batches[0] == ("amp1", [60, 100])
        assert batches[1] == ("amp1", [260, 300])

    def test_boundary_case(self):
        """Span exactly equal to max_distance stays in one batch."""
        amp_groups = {"amp1": [60, 259]}  # span = 259-60+1 = 200
        batches = split_positions_by_distance(amp_groups, max_distance=200)
        assert len(batches) == 1

    def test_single_position(self):
        """A single position forms one batch."""
        amp_groups = {"amp1": [100]}
        batches = split_positions_by_distance(amp_groups, max_distance=200)
        assert len(batches) == 1
        assert batches[0] == ("amp1", [100])

    def test_empty_amplicon_ignored(self):
        """Amplicons with empty position lists produce no batches."""
        amp_groups = {"amp1": [], "amp2": [100, 150]}
        batches = split_positions_by_distance(amp_groups, max_distance=200)
        assert len(batches) == 1
        assert batches[0] == ("amp2", [100, 150])


class TestBuildAmpDictFromVariants:
    """Test the build_amp_dict_from_variants function."""

    def test_merges_alleles_across_variants(self):
        """When two variants have different alts at the same position, both are recorded."""
        # Minimal mock pango_loader
        class MockLoader:
            _reconstructed_signatures = set()
            _sigs = {
                "VarA": {"241T", "297G"},
                "VarB": {"241A", "500C"},
            }
            def get_signature(self, name):
                return self._sigs.get(name, set())

        amp_dict = build_amp_dict_from_variants(["VarA", "VarB"], MockLoader())
        assert set(amp_dict[241]) == {"A", "T"}
        assert amp_dict[297] == ["G"]
        assert amp_dict[500] == ["C"]

    def test_deletions_excluded(self):
        """Deletion entries (ending with '-') are skipped."""
        class MockLoader:
            _reconstructed_signatures = set()
            def get_signature(self, name):
                return {"241T", "27469-"}

        amp_dict = build_amp_dict_from_variants(["V"], MockLoader())
        assert 241 in amp_dict
        assert 27469 not in amp_dict