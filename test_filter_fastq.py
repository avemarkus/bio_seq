import os
import pytest
import numpy as np
from filter_fastq import filter_fastq

class TestFilterFastq:
    """
    Test suite for the filter_fastq function in filter_fastq module.
    Tests cover functionality, error handling, logging, and edge cases.
    """
    
    @pytest.fixture
    def input_file_path(self, tmp_path):
        """
        Fixture to create a temporary FastQ file for testing.
        """
        content = (
            "@seq1\nGACCTTTCCGCAAGCTGTCGC\n+\nIIIIIIIIIIIIIIIIIIIII\n"
            "@seq2\nCATGGTGGCG\n+\nIIIIIIIIII\n"
            "@seq3\nC\n+\nI\n"
            "@seq4\nATCG\n+\n!!!!\n"
        )
        file_path = tmp_path / "example_fastq.fastq"
        file_path.write_text(content)
        return str(file_path)

    @pytest.fixture
    def log_file_path(self, tmp_path):
        """
        Fixture to clean up the log file after the test.
        """
        yield str(tmp_path / "filter.log")
        #if os.path.exists("filter.log"):
            #os.remove("filter.log")

    def test_output_exists(self, input_file_path):
        """
        Test that filter_fastq returns a non-empty dictionary with default parameters.
        """
        result = filter_fastq(input_file_path)
        assert isinstance(result, dict), "Result should be a dictionary"
        assert len(result) > 0, "Result dictionary should not be empty"

    def test_quality_threshold_output(self, input_file_path):
        """
        Test the correctness of output sequences with quality_threshold=38.
        """
        target_seqs = ['GACCTTTCCGCAAGCTGTCGC', 'CATGGTGGCG', 'C']
        result = filter_fastq(input_file_path, quality_threshold=38)
        assert isinstance(result, dict), "Result should be a dictionary"
        result_seqs = [seq for seq, _ in result.values()]
        assert result_seqs == target_seqs, f"Expected sequences {target_seqs}, but got {result_seqs}"

    def test_gc_bounds_filtering(self, input_file_path):
        """
        Test filtering by GC content (gc_bounds=(50.1, 100)).
        """
        target_seqs = ['GACCTTTCCGCAAGCTGTCGC', 'CATGGTGGCG', 'C']
        result = filter_fastq(input_file_path, gc_bounds=(50.1, 100))
        assert isinstance(result, dict), "Result should be a dictionary"
        result_seqs = [seq for seq, _ in result.values()]
        assert result_seqs == target_seqs, f"Expected sequences {target_seqs}, but got {result_seqs}"

    def test_length_bounds_filtering(self, input_file_path):
        """
        Test filtering by sequence length (length_bounds=(5, 15)).
        """
        target_seqs = ['CATGGTGGCG']
        result = filter_fastq(input_file_path, length_bounds=(5, 15))
        assert isinstance(result, dict), "Result should be a dictionary"
        result_seqs = [seq for seq, _ in result.values()]
        assert result_seqs == target_seqs, f"Expected sequences {target_seqs}, but got {result_seqs}"

    def test_empty_result(self, input_file_path):
        """
        Test that filter_fastq returns a string when no sequences pass the filter.
        """
        result = filter_fastq(input_file_path, quality_threshold=100)
        assert isinstance(result, str), "Result should be a string"
        assert "No sequences was found" in result, "Result should indicate no sequences found"

    def test_file_not_found_error(self):
        """
        Test that filter_fastq raises FileNotFoundError for a non-existent file.
        """
        with pytest.raises(FileNotFoundError, match="Fastq file 'non_existent.fastq' was not found"):
            filter_fastq("non_existent.fastq")

    def test_logging_output(self, input_file_path, log_file_path):
        """
        Test that filter_fastq creates a log file and writes expected messages.
        """
        result = filter_fastq(input_file_path)
        log_file = "filter.log"
        assert os.path.exists(log_file), "Log file 'filter.log' should be created"
        with open(log_file, 'r') as f:
            log_content = f.read()
            assert "Starting filtering of fastq file" in log_content, "Log should contain start message"
            assert f"Filtered {len(result)} sequences" in log_content, "Log should contain result message"

    def test_invalid_fastq_format(self, tmp_path):
        """
        Test that filter_fastq handles invalid FastQ format gracefully.
        """
        invalid_content = (
            "@seq1\nATCG\n+\n!!\n"
        )
        file_path = tmp_path / "invalid_fastq.fastq"
        file_path.write_text(invalid_content)
        result = filter_fastq(str(file_path))
        assert isinstance(result, str), "Result should be a string"
        assert "Invalid format" in result, "Result should indicate invalid format"