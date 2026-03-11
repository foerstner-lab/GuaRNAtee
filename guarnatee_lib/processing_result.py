"""
Data classes for processing results.
Encapsulates results from different processing stages.
"""
from dataclasses import dataclass
from typing import Optional
import pandas as pd


@dataclass
class ProcessingResult:
    """
    Results from processing a library group.

    Attributes:
        candidates_df: DataFrame with RNA candidates
        stats_df: DataFrame with processing statistics
        library_type: Type of library processed
        sample_count: Number of samples processed
    """
    candidates_df: pd.DataFrame
    stats_df: pd.DataFrame
    library_type: str
    sample_count: int

    def is_empty(self) -> bool:
        """Check if result contains any candidates"""
        return self.candidates_df.empty

    def get_candidate_count(self) -> int:
        """Get total number of candidates"""
        return len(self.candidates_df)

    def get_cluster_count(self) -> int:
        """Get number of unique clusters if cluster_id exists"""
        if "cluster_id" in self.candidates_df.columns:
            return self.candidates_df["cluster_id"].nunique()
        return 0
