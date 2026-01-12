import pytest
import tempfile
import shutil
import numpy as np

from unittest import mock
from pathlib import Path

from alphafold.analysis.attention_pipeline import run_pipeline
from alphafold.analysis import analyze_residue, plot_difference, process_attention


@pytest.fixture
def test_output_dir():
    """Create a temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def mock_base_args(test_output_dir):
    """Create base mock arguments for pipeline tests."""
    return {
        "query_seq_path": "tests/test_data/test.fasta",
        "query_attn_dir": "tests/test_data/query_attn",
        "query_name": "test_protein",
        "target_seq_path": None,
        "target_attn_dir": None,
        "target_name": None,
        "alignment_path": None,
        "save_path": test_output_dir,
        "query_highlight_indices": None,
        "target_highlight_indices": None,
        "query_highlight_color": "#AE0639",
        "target_highlight_color": "#1f77b4",
    }


class TestPipelineBasic:
    """Tests for basic pipeline functionality."""

    @mock.patch("alphafold.analysis.process_attention.read_sequence_file")
    @mock.patch("alphafold.analysis.process_attention.get_n")
    @mock.patch("alphafold.analysis.process_attention.get_attention")
    @mock.patch("alphafold.analysis.process_attention.average")
    @mock.patch("alphafold.analysis.process_attention.min_max")
    @mock.patch("alphafold.analysis.analyze_residue.find_important")
    @mock.patch("alphafold.analysis.analyze_residue.find_highest_attention")
    @mock.patch("alphafold.analysis.plot_difference.plot_attention")
    def test_run_pipeline_single_protein(
        self,
        mock_plot_attn,
        mock_find_highest,
        mock_find_important,
        mock_min_max,
        mock_average,
        mock_get_attn,
        mock_get_n,
        mock_read_seq,
        mock_base_args,
    ):
        """Test basic pipeline execution without target."""
        mock_read_seq.return_value = "ACDEFGHIKLMNPQRSTVWY"
        mock_get_n.return_value = 1
        mock_get_attn.return_value = np.random.rand(20, 20)
        mock_average.return_value = np.random.rand(20)
        mock_min_max.return_value = np.random.rand(20)
        mock_find_important.return_value = np.array([5, 10, 15])

        run_pipeline(**mock_base_args)

        assert Path(mock_base_args["save_path"]).exists()

        mock_read_seq.assert_called_once()
        mock_get_attn.assert_called_once()
        mock_plot_attn.assert_called_once()

    @mock.patch("alphafold.analysis.process_attention.read_sequence_file")
    @mock.patch("alphafold.analysis.process_attention.get_n")
    @mock.patch("alphafold.analysis.process_attention.get_attention")
    @mock.patch("alphafold.analysis.process_attention.average")
    @mock.patch("alphafold.analysis.process_attention.min_max")
    @mock.patch("alphafold.analysis.analyze_residue.find_important")
    @mock.patch("alphafold.analysis.analyze_residue.find_highest_attention")
    @mock.patch("alphafold.analysis.analyze_residue.blosum_scores")
    @mock.patch("alphafold.analysis.analyze_residue.calculate_differences")
    @mock.patch("alphafold.analysis.plot_difference.plot_attention")
    @mock.patch("alphafold.analysis.plot_difference.plot_difference")
    def test_run_pipeline_with_target_equal_length(
        self,
        mock_plot_diff,
        mock_plot_attn,
        mock_calc_diff,
        mock_blosum,
        mock_find_highest,
        mock_find_important,
        mock_min_max,
        mock_average,
        mock_get_attn,
        mock_get_n,
        mock_read_seq,
        mock_base_args,
    ):
        """Test pipeline with equal-length target sequence."""
        seq = "ACDEFGHIKLMNPQRSTVWY"
        mock_read_seq.return_value = seq
        mock_get_n.return_value = 1
        mock_get_attn.return_value = np.random.rand(20, 20)
        mock_average.return_value = np.random.rand(20)
        mock_min_max.return_value = np.random.rand(20)
        mock_find_important.return_value = np.array([5, 10])
        mock_blosum.return_value = (np.random.rand(20), np.random.rand(20))
        mock_calc_diff.return_value = (np.random.rand(20), np.random.rand(20))

        mock_base_args["target_name"] = "target_protein"
        mock_base_args["target_seq_path"] = "tests/test_data/target.fasta"
        mock_base_args["target_attn_dir"] = "tests/test_data/target_attn"

        run_pipeline(**mock_base_args)

        # Check pipeline is ran twice (for query and target)
        assert mock_read_seq.call_count >= 2
        mock_blosum.assert_called_once()
        mock_calc_diff.assert_called_once()

    @mock.patch("alphafold.analysis.process_attention.read_sequence_file")
    @mock.patch("alphafold.analysis.process_attention.get_n")
    @mock.patch("alphafold.analysis.process_attention.get_attention")
    @mock.patch("alphafold.analysis.process_attention.average")
    @mock.patch("alphafold.analysis.process_attention.min_max")
    @mock.patch("alphafold.analysis.analyze_residue.find_important")
    def test_run_pipeline_length_mismatch_no_alignment(
        self,
        mock_find_important,
        mock_min_max,
        mock_average,
        mock_get_attn,
        mock_get_n,
        mock_read_seq,
        mock_base_args,
    ):
        """Test pipeline fails with length mismatch and no alignment."""
        mock_read_seq.side_effect = ["ACDEF", "ACDEFGHIKLMNPQRSTVWY"]
        mock_get_n.return_value = 1
        mock_get_attn.return_value = np.random.rand(20, 20)
        mock_average.return_value = np.random.rand(20)
        mock_min_max.return_value = np.random.rand(20)

        mock_base_args["target_name"] = "target_protein"
        mock_base_args["target_seq_path"] = "tests/test_data/target.fasta"
        mock_base_args["target_attn_dir"] = "tests/test_data/target_attn"

        # Should exit with error
        with pytest.raises(SystemExit):
            run_pipeline(**mock_base_args)


class TestFindImportant:
    """Tests for residue importance analysis."""

    def test_find_important_basic(self):
        """Test finding important residues with basic input."""
        attention = np.array([0.1, 0.5, 0.2, 0.9, 0.4])
        zscores = np.array([1.0, 2.0, 1.5, 3.0, 2.5])

        important = analyze_residue.find_important(attention, zscores)

        assert isinstance(important, np.ndarray), "Should return numpy array"
        assert important.shape == attention.shape, "Output shape should match input"

        non_zero_count = np.count_nonzero(important)
        assert non_zero_count <= 5, "Should have at most 5 non-zero values"
        assert non_zero_count == min(5, len(attention)), "Should identify top indices"

        assert important[3] > 0, "Highest attention index should be non-zero"

    def test_find_important_empty(self):
        """Test finding important residues with empty arrays."""
        attention = np.array([])
        zscores = np.array([])

        important = analyze_residue.find_important(attention, zscores)

        assert isinstance(important, np.ndarray), "Should return numpy array"
        assert len(important) == 0, "Should return empty array for empty input"

    def test_find_important_single_residue(self):
        """Test finding important residues with single residue."""
        attention = np.array([0.5])
        zscores = np.array([2.0])

        important = analyze_residue.find_important(attention, zscores)

        assert isinstance(important, np.ndarray), "Should return numpy array"
        assert important.shape == (1,), "Should preserve shape"
        assert important[0] == 0.5, "Single residue should be marked as important"

    def test_find_important_all_below_threshold(self):
        """Test with low attention values."""
        attention = np.array([0.01, 0.01, 0.01, 0.01])
        zscores = np.array([0.1, 0.1, 0.1, 0.1])

        important = analyze_residue.find_important(attention, zscores)

        non_zero_count = np.count_nonzero(important)
        assert non_zero_count == 4, "Should select all 4 residues (less than top-5)"
        np.testing.assert_array_almost_equal(important, attention)

    def test_find_important_large_array(self):
        """Test with array larger than top-5 threshold."""
        attention = np.array([0.1, 0.5, 0.2, 0.9, 0.4, 0.3, 0.7, 0.6, 0.8, 0.1])

        important = analyze_residue.find_important(attention)

        non_zero_count = np.count_nonzero(important)
        assert non_zero_count == 5, "Should select exactly 5 residues"

        expected_top_indices = {1, 3, 6, 7, 8}
        actual_top_indices = set(np.where(important > 0)[0])
        assert (
            actual_top_indices == expected_top_indices
        ), "Should identify correct top-5 indices"

    def test_find_important_values_preserved(self):
        """Test that original attention values are preserved at selected indices."""
        attention = np.array([0.1, 0.9, 0.2, 0.8, 0.3])

        important = analyze_residue.find_important(attention)

        # Non-zero entries should have original attention values
        for idx in np.where(important > 0)[0]:
            assert (
                important[idx] == attention[idx]
            ), f"Value at index {idx} should be preserved"


class TestBlosumScores:
    """Tests for BLOSUM score calculation."""

    def test_blosum_scores_identical(self):
        """Test BLOSUM score for identical sequences."""
        seq = "ACDEFGHIK"

        scores = analyze_residue.blosum_scores(seq, seq)

        assert len(scores) == 2, "Should return tuple of two score arrays"
        assert isinstance(scores[0], np.ndarray), "Scores should be numpy arrays"

    def test_blosum_scores_different(self):
        """Test BLOSUM score for different sequences."""
        seq1 = "ACDEFGHIK"
        seq2 = "ACDFEGHIK"

        scores = analyze_residue.blosum_scores(seq1, seq2)

        assert len(scores) == 2, "Should return tuple of (query_scores, target_scores)"

    def test_blosum_scores_empty(self):
        """Test BLOSUM score with empty sequences."""
        scores = analyze_residue.blosum_scores("", "")

        assert len(scores) == 2, "Should handle empty sequences"

    def test_blosum_scores_different_lengths(self):
        """Test BLOSUM score with different length sequences."""
        seq1 = "ACDEF"
        seq2 = "ACDEFGHIK"

        with pytest.raises((ValueError, ValueError)):
            analyze_residue.blosum_scores(seq1, seq2)


class TestCalculateDifference:
    """Tests for attention difference calculation."""

    def test_calculate_difference_basic(self):
        """Test calculating attention difference with basic input."""
        attention1 = np.array([0.1, 0.4, 0.3, 0.8])
        attention2 = np.array([0.2, 0.3, 0.5, 0.7])

        diff1, diff2 = analyze_residue.calculate_differences(
            scores1=np.array([1, 2, 3, 4]),
            scores2=np.array([1, 2, 3, 4]),
            attention1=attention1,
            attention2=attention2,
            gaps1=[],
            gaps2=[],
            bool_alignment=False,
        )

        assert isinstance(diff1, list), "Should return list"
        assert isinstance(diff2, list), "Should return list"
        assert len(diff1) == len(attention1), "Output length should match input"
        assert len(diff2) == len(attention2), "Output length should match input"

        # When attention arrays differ, differences should be non-zero
        assert not np.allclose(
            diff1, 0, atol=1e-6
        ), "Should have positive differences where attention1 > attention2"
        assert not np.allclose(
            diff2, 0, atol=1e-6
        ), "Should have positive differences where attention2 > attention1"

    def test_calculate_difference_zero(self):
        """Test calculating difference when arrays are identical."""
        attention = np.array([0.1, 0.4, 0.3, 0.8])
        scores = np.array([1, 2, 3, 4])

        diff1, diff2 = analyze_residue.calculate_differences(
            scores1=scores,
            scores2=scores,
            attention1=attention,
            attention2=attention,
            gaps1=[],
            gaps2=[],
            bool_alignment=False,
        )

        assert isinstance(diff1, list), "Should return list"
        assert isinstance(diff2, list), "Should return list"
        assert np.allclose(
            diff1, 0, atol=1e-6
        ), "Identical arrays should give zero difference"
        assert np.allclose(
            diff2, 0, atol=1e-6
        ), "Identical arrays should give zero difference"

    def test_calculate_difference_shape_mismatch(self):
        """Test calculating difference with mismatched array shapes."""
        attention1 = np.array([0.1, 0.4, 0.3])
        attention2 = np.array([0.2, 0.3])
        scores1 = np.array([1, 2, 3])
        scores2 = np.array([1, 2])

        with pytest.raises((ValueError, IndexError)):
            analyze_residue.calculate_differences(
                scores1=scores1,
                scores2=scores2,
                attention1=attention1,
                attention2=attention2,
                gaps1=[],
                gaps2=[],
                bool_alignment=False,
            )


class TestPlotAttention:
    """Tests for plotting attention scores."""

    def test_plot_attention_basic(self, test_output_dir):
        """Test plotting attention with basic input."""
        attention_scores = np.array([0.1, 0.5, 0.2, 0.9, 0.4])
        highlighted_scores = np.array([0.0, 0.5, 0.0, 0.9, 0.0])
        protein_name = "test_protein"
        sequence = "ACDEF"

        plot_difference.plot_attention(
            attention_scores=attention_scores,
            highlighted_scores=highlighted_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence=sequence,
        )

        output_file = Path(test_output_dir) / f"{protein_name}_average_attention.png"
        assert output_file.exists(), "Output PNG should be created"

    def test_plot_attention_without_sequence(self, test_output_dir):
        """Test plotting attention without sequence labels."""
        attention_scores = np.array([0.1, 0.5, 0.2, 0.9, 0.4])
        highlighted_scores = None
        protein_name = "test_protein_no_seq"

        plot_difference.plot_attention(
            attention_scores=attention_scores,
            highlighted_scores=highlighted_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence=None,
        )

        output_file = Path(test_output_dir) / f"{protein_name}_average_attention.png"
        assert output_file.exists(), "Output PNG should be created without sequence"

    def test_plot_attention_empty_array(self, test_output_dir):
        """Test plotting attention with empty array."""
        attention_scores = np.array([])
        highlighted_scores = np.array([])
        protein_name = "test_empty"

        plot_difference.plot_attention(
            attention_scores=attention_scores,
            highlighted_scores=highlighted_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence="",
        )

        # Should still create output without error
        output_file = Path(test_output_dir) / f"{protein_name}_average_attention.png"
        assert output_file.exists(), "Output should be created even for empty arrays"

    def test_plot_attention_creates_directory(self):
        """Test that plot_attention creates output directory if it doesn't exist."""
        attention_scores = np.array([0.1, 0.5, 0.2])
        protein_name = "test_create_dir"
        output_dir = Path(tempfile.mkdtemp()) / "nested" / "dir"

        plot_difference.plot_attention(
            attention_scores=attention_scores,
            highlighted_scores=None,
            protein_name=protein_name,
            output_dir=output_dir,
            sequence=None,
        )

        assert output_dir.exists(), "Output directory should be created"
        output_file = output_dir / f"{protein_name}_average_attention.png"
        assert output_file.exists(), "PNG should be saved in created directory"

        # Cleanup
        shutil.rmtree(output_dir.parent, ignore_errors=True)


class TestPlotDifference:
    """Tests for plotting attention differences."""

    def test_plot_difference_basic(self, test_output_dir):
        """Test plotting attention difference with basic input."""
        attn_diff_scores = np.array([0.1, -0.4, 0.3, -0.8, 0.2])
        protein_name = "test_diff"
        sequence = "ACDEF"

        plot_difference.plot_difference(
            attn_diff_scores=attn_diff_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence=sequence,
        )

        output_file = Path(test_output_dir) / f"{protein_name}_attention_difference.png"
        assert output_file.exists(), "Output PNG should be created"

    def test_plot_difference_with_highlights(self, test_output_dir):
        """Test plotting with highlight positions."""
        attn_diff_scores = np.array([-0.1, -0.5, -0.2, -0.9, -0.4])
        protein_name = "test_diff_highlight"
        sequence = "ACDEF"
        query_highlights = [1, 3]
        target_highlights = [2, 4]

        plot_difference.plot_difference(
            attn_diff_scores=attn_diff_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence=sequence,
            query_highlight_positions=query_highlights,
            target_highlight_positions=target_highlights,
            query_highlight_color="#AE0639",
            target_highlight_color="#1f77b4",
        )

        output_file = Path(test_output_dir) / f"{protein_name}_attention_difference.png"
        assert output_file.exists(), "Output PNG with highlights should be created"

    def test_plot_difference_no_sequence(self, test_output_dir):
        """Test plotting difference without sequence labels."""
        attn_diff_scores = np.array([-0.1, -0.5, -0.2, -0.9])
        protein_name = "test_diff_no_seq"

        plot_difference.plot_difference(
            attn_diff_scores=attn_diff_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence=None,
        )

        output_file = Path(test_output_dir) / f"{protein_name}_attention_difference.png"
        assert output_file.exists(), "Output should be created without sequence"

    def test_plot_difference_empty_array(self, test_output_dir):
        """Test plotting difference with empty array."""
        attn_diff_scores = np.array([])
        protein_name = "test_diff_empty"

        plot_difference.plot_difference(
            attn_diff_scores=attn_diff_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence="",
        )

        output_file = Path(test_output_dir) / f"{protein_name}_attention_difference.png"
        assert output_file.exists(), "Output should be created for empty arrays"

    def test_plot_difference_positive_values_zeroed(self, test_output_dir):
        """Test that positive values are converted to zero in plot."""
        attn_diff_scores = np.array([0.5, -0.5, 0.3, -0.3])
        protein_name = "test_positive_zero"

        plot_difference.plot_difference(
            attn_diff_scores=attn_diff_scores,
            protein_name=protein_name,
            output_dir=test_output_dir,
            sequence="ACDE",
        )

        output_file = Path(test_output_dir) / f"{protein_name}_attention_difference.png"
        assert (
            output_file.exists()
        ), "Output should be created with positive values zeroed"
