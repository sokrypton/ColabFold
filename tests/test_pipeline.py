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


class TestGetN:
    """Tests for sequence length inference from attention files."""

    def test_get_n_single_shape(self, test_output_dir):
        """Test get_n with uniform array shapes."""
        # Create test .npy files with shape (20, 4, 20, 20)
        for i in range(3):
            arr = np.random.rand(20, 4, 20, 20)
            np.save(Path(test_output_dir) / f"attention_{i}.npy", arr)

        n = process_attention.get_n(test_output_dir)

        assert n == 20, "Should return first dimension (20)"

    def test_get_n_multiple_shapes_returns_most_frequent(self, test_output_dir):
        """Test get_n returns most frequent shape when multiple exist."""
        # Create 3 files with shape (20, 4, 20, 20)
        for i in range(3):
            arr = np.random.rand(20, 4, 20, 20)
            np.save(Path(test_output_dir) / f"attention_{i}.npy", arr)

        # Create 1 file with shape (30, 4, 30, 30)
        arr = np.random.rand(30, 4, 30, 30)
        np.save(Path(test_output_dir) / f"attention_outlier.npy", arr)

        n = process_attention.get_n(test_output_dir)

        assert n == 20, "Should return most frequent shape's first dimension"

    def test_get_n_no_npy_files_exits(self, test_output_dir):
        """Test get_n exits when no .npy files found."""
        # Create empty directory
        assert not any(Path(test_output_dir).glob("*.npy"))

        with pytest.raises(SystemExit):
            process_attention.get_n(test_output_dir)

    def test_get_n_corrupted_file_skipped(self, test_output_dir):
        """Test get_n skips corrupted files and uses valid ones."""
        # Create valid file
        arr = np.random.rand(20, 4, 20, 20)
        np.save(Path(test_output_dir) / "attention_valid.npy", arr)

        # Create corrupted file (just write invalid bytes)
        with open(Path(test_output_dir) / "attention_corrupt.npy", "wb") as f:
            f.write(b"invalid npy data")

        n = process_attention.get_n(test_output_dir)

        assert n == 20, "Should skip corrupted file and use valid one"


class TestGetAttention:
    """Tests for attention spectrum loading and processing."""

    def test_get_attention_basic(self, test_output_dir):
        """Test loading attention from .npy files."""
        n = 20
        # Create attention files with shape (n, 4, n, n)
        for i in range(2):
            arr = np.random.rand(n, 4, n, n).astype(np.float16)
            np.save(Path(test_output_dir) / f"attention_{i}.npy", arr)

        spectrum = process_attention.get_attention(test_output_dir, n)

        assert spectrum.shape == (2, n), "Should return (num_files, n)"
        assert spectrum.dtype == np.float16, "Should preserve data type"

    def test_get_attention_filters_wrong_shape(self, test_output_dir):
        """Test that arrays with wrong shape are filtered out."""
        n = 20
        # Create valid file
        arr_valid = np.random.rand(n, 4, n, n).astype(np.float16)
        np.save(Path(test_output_dir) / f"attention_0.npy", arr_valid)

        # Create invalid shape file
        arr_invalid = np.random.rand(n, 2, n, n).astype(np.float16)
        np.save(Path(test_output_dir) / f"attention_1.npy", arr_invalid)

        spectrum = process_attention.get_attention(test_output_dir, n)

        assert spectrum.shape[0] == 1, "Should only include correct shape files"

    def test_get_attention_files_sorted_numerically(self, test_output_dir):
        """Test that files are processed in numerical order."""
        n = 5
        # Create files in non-sequential order
        for i in [2, 0, 1]:
            arr = np.full((n, 4, n, n), i, dtype=np.float16)
            np.save(Path(test_output_dir) / f"attention_{i}.npy", arr)

        spectrum = process_attention.get_attention(test_output_dir, n)

        # Files should be ordered 0, 1, 2
        assert spectrum.shape[0] == 3, "Should process all 3 files"


class TestAverage:
    """Tests for attention spectrum averaging."""

    def test_average_basic(self):
        """Test averaging attention spectrum."""
        spectrum = np.array([[0.1, 0.2, 0.3], [0.3, 0.2, 0.1], [0.2, 0.2, 0.2]])

        avg = process_attention.average(spectrum)

        expected = np.array([0.2, 0.2, 0.2])
        np.testing.assert_array_almost_equal(avg, expected)

    def test_average_single_file(self):
        """Test averaging with single file."""
        spectrum = np.array([[0.5, 0.7, 0.3]])

        avg = process_attention.average(spectrum)

        np.testing.assert_array_almost_equal(avg, spectrum[0])

    def test_average_preserves_shape(self):
        """Test that averaging preserves residue dimension."""
        spectrum = np.random.rand(10, 20)

        avg = process_attention.average(spectrum)

        assert avg.shape == (20,), "Should return 1D array of residues"


class TestMinMax:
    """Tests for min-max normalization."""

    def test_min_max_basic(self):
        """Test min-max normalization."""
        data = np.array([0.2, 0.4, 0.6, 0.8])

        normalized = process_attention.min_max(data)

        expected = np.array([0, 1/3, 2/3, 1.0])
        np.testing.assert_array_almost_equal(normalized, expected)

    def test_min_max_with_zeros(self):
        """Test that zero values stay zero."""
        data = np.array([0, 0.5, 0, 1.0])

        normalized = process_attention.min_max(data)

        assert normalized[0] == 0, "Zero should remain zero"
        assert normalized[2] == 0, "Zero should remain zero"
        assert normalized[3] == 1, "Max should be 1"

    def test_min_max_constant_array(self):
        """Test min-max with constant non-zero array."""
        data = np.array([0.5, 0.5, 0.5, 0.5])

        normalized = process_attention.min_max(data)

        # When all non-zero values are equal, division by zero is avoided
        # Function should handle gracefully
        assert len(normalized) == 4, "Should return same length"

    def test_min_max_single_nonzero(self):
        """Test with single non-zero value."""
        data = np.array([0, 0.5, 0, 0])

        normalized = process_attention.min_max(data)
        print(f"Normalized single non-zero: {normalized}")

        assert normalized[1] == 1, "Single non-zero should normalize to 1"
        assert np.sum(normalized) == 1, "Only that value should be non-zero"


class TestReadSequenceFile:
    """Tests for FASTA sequence reading."""

    def test_read_sequence_single_record(self, test_output_dir):
        """Test reading single FASTA record."""
        fasta_path = Path(test_output_dir) / "test.fasta"
        fasta_content = """>seq1\nACDEFGHIKLMNPQRSTVWY"""
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        seq = process_attention.read_sequence_file(str(fasta_path))

        assert seq == "ACDEFGHIKLMNPQRSTVWY"

    def test_read_sequence_multiline(self, test_output_dir):
        """Test reading sequence split across multiple lines."""
        fasta_path = Path(test_output_dir) / "test_multiline.fasta"
        fasta_content = """>seq1\nACDEF\nGHIKLM\nNPQRS\nTVW\nY"""
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        seq = process_attention.read_sequence_file(str(fasta_path))

        assert seq == "ACDEFGHIKLMNPQRSTV" + "WY"

    def test_read_sequence_multiple_records(self, test_output_dir):
        """Test reading first sequence from multi-record FASTA."""
        fasta_path = Path(test_output_dir) / "test_multi.fasta"
        fasta_content = """>seq1\nACDEFGHIK\n>seq2\nLMNPQRSTVWY"""
        with open(fasta_path, "w") as f:
            f.write(fasta_content)

        seq = process_attention.read_sequence_file(str(fasta_path))

        assert seq == "ACDEFGHIK", "Should read first sequence only"


class TestReadAlignment:
    """Tests for alignment file reading."""

    def test_read_alignment_basic(self, test_output_dir):
        """Test reading aligned sequences from file."""
        align_path = Path(test_output_dir) / "alignment.txt"
        align_content = """protein1\nACDEF-GHIK\nprotein2\nACDE--FHIK"""
        with open(align_path, "w") as f:
            f.write(align_content)

        seq1, seq2 = process_attention.read_alignment(
            "protein1", "protein2", str(align_path)
        )

        assert seq1 == "ACDEF-GHIK"
        assert seq2 == "ACDE--FHIK"

    def test_read_alignment_case_insensitive(self, test_output_dir):
        """Test that alignment matching is case-insensitive."""
        align_path = Path(test_output_dir) / "alignment_case.txt"
        align_content = """PROTEIN1\nACDEFGHIK\nProtein2\nLMNPQRSTVWY"""
        with open(align_path, "w") as f:
            f.write(align_content)

        seq1, seq2 = process_attention.read_alignment(
            "protein1", "protein2", str(align_path)
        )

        assert seq1 == "ACDEFGHIK"
        assert seq2 == "LMNPQRSTVWY"

    def test_read_alignment_not_found(self, test_output_dir):
        """Test when protein not found in alignment."""
        align_path = Path(test_output_dir) / "alignment_partial.txt"
        align_content = """protein1\nACDEFGHIK"""
        with open(align_path, "w") as f:
            f.write(align_content)

        seq1, seq2 = process_attention.read_alignment(
            "protein1", "missing", str(align_path)
        )

        assert seq1 == "ACDEFGHIK"
        assert seq2 == "", "Missing protein should return empty string"


class TestAlignAttention:
    """Tests for attention alignment mapping."""

    def test_align_attention_no_gaps(self):
        """Test attention alignment with no gaps."""
        attention1 = np.array([0.1, 0.2, 0.3, 0.4])
        attention2 = np.array([0.2, 0.3, 0.4, 0.5])
        seq1 = "ACDE"
        seq2 = "ACDE"

        aligned1, aligned2, gaps1, gaps2 = process_attention.align_attention(
            attention1, attention2, seq1, seq2
        )

        np.testing.assert_array_almost_equal(aligned1, attention1)
        np.testing.assert_array_almost_equal(aligned2, attention2)
        assert gaps1 == []
        assert gaps2 == []

    def test_align_attention_with_gaps_seq1(self):
        """Test alignment when seq1 has gaps (aligned sequences same length)."""
        attention1 = np.array([0.1, 0.2, 0.3, 0.4])
        attention2 = np.array([0.2, 0.3, 0.4, 0.5, 0.6])
        seq1 = "AC-DE"
        seq2 = "ACDEF"

        aligned1, aligned2, gaps1, gaps2 = process_attention.align_attention(
            attention1, attention2, seq1, seq2
        )

        assert aligned1[2] == 0, "Gap position should be zero"
        assert 2 in gaps1, "Gap index should be recorded"
        assert len(aligned1) == len(seq1), "aligned1 length should match seq1"
        assert len(aligned2) == len(seq2), "aligned2 length should match seq2"

    def test_align_attention_with_gaps_seq2(self):
        """Test alignment when seq2 has gaps (aligned sequences same length)."""
        attention1 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        attention2 = np.array([0.2, 0.3, 0.4, 0.5])
        seq1 = "ACDEF"
        seq2 = "AC-DE"

        aligned1, aligned2, gaps1, gaps2 = process_attention.align_attention(
            attention1, attention2, seq1, seq2
        )

        assert aligned2[2] == 0, "Gap position in seq2 should be zero"
        assert 2 in gaps2, "Gap index should be recorded in gaps2"
        assert len(aligned1) == len(seq1), "aligned1 length should match seq1"
        assert len(aligned2) == len(seq2), "aligned2 length should match seq2"

    def test_align_attention_length_scaling(self):
        """Test that attention is scaled by length ratio."""
        attention1 = np.array([0.5, 0.5])
        attention2 = np.array([0.5, 0.5, 0.5, 0.5])
        seq1 = "AC"
        seq2 = "ACDE"

        aligned1, aligned2, _, _ = process_attention.align_attention(
            attention1, attention2, seq1, seq2
        )

        # seq1 is shorter, so ratio = 2/4 = 0.5
        # aligned1 should be scaled down
        expected1 = np.array([0.25, 0.25])
        np.testing.assert_array_almost_equal(aligned1, expected1)

    def test_align_attention_output_shapes(self):
        """Test that output shapes match input sequences."""
        attention1 = np.array([0.1, 0.2, 0.3, 0.5])
        attention2 = np.array([0.4, 0.5, 0.6, 0.7, 0.8])
        seq1 = "A-CDE"
        seq2 = "AC-DEF"

        aligned1, aligned2, gaps1, gaps2 = process_attention.align_attention(
            attention1, attention2, seq1, seq2
        )

        assert len(aligned1) == len(seq1), "aligned1 length should match seq1"
        assert len(aligned2) == len(seq2), "aligned2 length should match seq2"
        assert isinstance(gaps1, list), "gaps1 should be a list"
        assert isinstance(gaps2, list), "gaps2 should be a list"
