import os
import tempfile
import pytest
import pandas as pd
from finaletoolkit.utils._sort_bed_file import sort_bed_file

class TestSortBedFile:
    
    @pytest.fixture
    def temp_bed_file(self):
        """Create a temporary BED file for testing."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as f:
            # Write test data with mixed chromosome order
            f.write("chr2\t150\t200\tfeature1\t1\t+\n")
            f.write("chr1\t100\t200\tfeature2\t2\t-\n")
            f.write("chr10\t50\t100\tfeature3\t3\t+\n")
            f.write("chrX\t10\t20\tfeature4\t4\t-\n")
            f.write("chr1\t100\t150\tfeature5\t5\t+\n")  # Same start as feature2, earlier end
            f.write("chr1\t100\t300\tfeature6\t6\t-\n")  # Same start as feature2, later end
            f.write("# This is a comment\n")              # Comment line
            f.write("track name=test\n")                  # Track line
            filename = f.name
        
        yield filename
        # Clean up the temp file after the test
        os.unlink(filename)
    
    @pytest.fixture
    def temp_chrom_sizes_file(self):
        """Create a temporary chrom.sizes file for testing."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.sizes') as f:
            # Custom chromosome order: chr10, chr2, chr1, chrX
            f.write("chr10\t135534747\n")
            f.write("chr2\t243199373\n")
            f.write("chr1\t248956422\n")
            f.write("chrX\t156040895\n")
            filename = f.name
        
        yield filename
        # Clean up the temp file after the test
        os.unlink(filename)
    
    def test_default_sorting(self, temp_bed_file):
        """Test default sorting by chromosome, start, and end."""
        sorted_entries = list(sort_bed_file(temp_bed_file))
        
        # Expected order: chr1 (sorted by start, then end), chr2, chr10, chrX
        assert len(sorted_entries) == 6  # Excluding comment and track lines
        
        # Check chromosome order: chr1, chr2, chr10, chrX
        assert sorted_entries[0][0] == "chr1"
        assert sorted_entries[1][0] == "chr1"
        assert sorted_entries[2][0] == "chr1"
        assert sorted_entries[3][0] == "chr2"
        assert sorted_entries[4][0] == "chr10"
        assert sorted_entries[5][0] == "chrX"
        
        # Check tie-breaking by end coordinate for entries with same chrom and start
        # The three chr1 entries all have start=100, should be ordered by end
        assert sorted_entries[0][0] == "chr1" and sorted_entries[0][1] == "100" and sorted_entries[0][2] == "150"
        assert sorted_entries[1][0] == "chr1" and sorted_entries[1][1] == "100" and sorted_entries[1][2] == "200"
        assert sorted_entries[2][0] == "chr1" and sorted_entries[2][1] == "100" and sorted_entries[2][2] == "300"
    
    def test_chrom_sizes_sorting(self, temp_bed_file, temp_chrom_sizes_file):
        """Test sorting using a custom chrom.sizes file."""
        sorted_entries = list(sort_bed_file(temp_bed_file, temp_chrom_sizes_file))
        
        # Expected order based on chrom.sizes: chr10, chr2, chr1, chrX
        assert len(sorted_entries) == 6  # Excluding comment and track lines
        
        # Check chromosome order according to the chrom.sizes file
        assert sorted_entries[0][0] == "chr10"
        assert sorted_entries[1][0] == "chr2"
        assert sorted_entries[2][0] == "chr1"
        assert sorted_entries[3][0] == "chr1"
        assert sorted_entries[4][0] == "chr1"
        assert sorted_entries[5][0] == "chrX"
        
        # Check tie-breaking by end coordinate for entries with same chrom and start
        # The three chr1 entries all have start=100, should still be ordered by end
        assert sorted_entries[2][0] == "chr1" and sorted_entries[2][1] == "100" and sorted_entries[2][2] == "150"
        assert sorted_entries[3][0] == "chr1" and sorted_entries[3][1] == "100" and sorted_entries[3][2] == "200"
        assert sorted_entries[4][0] == "chr1" and sorted_entries[4][1] == "100" and sorted_entries[4][2] == "300"
    
    def test_complete_entry_preservation(self, temp_bed_file):
        """Test that all fields from the original BED file are preserved."""
        sorted_entries = list(sort_bed_file(temp_bed_file))
        
        # Check that a complete entry with all fields is preserved
        sample_entry = next(entry for entry in sorted_entries if entry[3] == "feature1")
        assert len(sample_entry) == 6
        assert sample_entry == ("chr2", "150", "200", "feature1", "1", "+")
    
    def test_empty_file(self):
        """Test behavior with an empty file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as f:
            filename = f.name
        
        # Should not raise an exception and yield no entries
        sorted_entries = list(sort_bed_file(filename))
        assert len(sorted_entries) == 0
        
        os.unlink(filename)
    
    def test_malformed_entries(self):
        """Test handling of malformed BED entries."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as f:
            f.write("chr1\t100\t200\tgood_entry\t1\t+\n")
            f.write("chr2\tnot_a_number\t200\tbad_entry\t2\t-\n")  # Bad start coordinate
            f.write("chr3\t300\tnot_a_number\tbad_entry\t3\t+\n")  # Bad end coordinate
            f.write("chr4\t400\t500\tgood_entry2\t4\t-\n")
            filename = f.name
        
        # Should skip the malformed entries
        sorted_entries = list(sort_bed_file(filename))
        assert len(sorted_entries) == 2
        assert sorted_entries[0][3] == "good_entry"
        assert sorted_entries[1][3] == "good_entry2"
        
        os.unlink(filename)
    
    def test_integration_with_pandas(self, temp_bed_file):
        """Test integration with pandas for further data manipulation."""
        sorted_entries = list(sort_bed_file(temp_bed_file))
        
        # Convert to pandas DataFrame for further analysis
        columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        df = pd.DataFrame(sorted_entries, columns=columns)
        
        # Verify DataFrame operations work correctly
        assert len(df) == 6
        assert df['chrom'].iloc[0] == "chr1"
        assert int(df['start'].iloc[0]) == 100
        assert int(df['end'].iloc[0]) == 150