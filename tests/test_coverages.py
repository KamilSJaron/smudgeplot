# tests/test_coverages.py
import pytest
import pandas as pd
import numpy as np
from smudgeplot.smudgeplot import Coverages
from smudgeplot.exceptions import InvalidCoverageDataError

class TestCoveragesInit:
    """Test Coverages initialization."""
    
    def test_valid_initialization(self):
        """Test with valid data."""
        df = pd.DataFrame({
            'covA': [10, 20, 30],
            'covB': [5, 10, 15],
            'freq': [100, 200, 300]
        })
        cov = Coverages(df)
        assert len(cov.cov_tab) == 3
        assert cov.total_kmers is None
    
    def test_missing_columns(self):
        """Test rejection of invalid columns."""
        df = pd.DataFrame({'covA': [10], 'freq': [100]})
        with pytest.raises(InvalidCoverageDataError, match="must contain columns"):
            Coverages(df)
    
    def test_empty_dataframe(self):
        """Test rejection of empty data."""
        df = pd.DataFrame(columns=['covA', 'covB', 'freq'])
        with pytest.raises(InvalidCoverageDataError, match="cannot be empty"):
            Coverages(df)
    
    def test_negative_values(self):
        """Test rejection of negative coverages."""
        df = pd.DataFrame({
            'covA': [10, -5],
            'covB': [5, 10],
            'freq': [100, 200]
        })
        with pytest.raises(InvalidCoverageDataError, match="non-negative"):
            Coverages(df)

class TestCoveragesCountKmers:
    """Test k-mer counting functionality."""
    
    def test_count_with_errors(self):
        """Test counting with error k-mers."""
        df = pd.DataFrame({
            'covA': [10, 20, 30],
            'covB': [5, 10, 15],
            'freq': [100, 200, 300],
            'smudge': [1, 1, -1]  # Last one is error
        })
        cov = Coverages(df)
        cov.count_kmers()
        
        assert cov.total_kmers == 600
        assert cov.total_error_kmers == 300
        assert cov.error_fraction == pytest.approx(0.5)

# Run tests
if __name__ == '__main__':
    pytest.main([__file__, '-v'])