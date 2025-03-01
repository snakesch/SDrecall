import pandas as pd
import logging
import numpy as np

from src.log import logger

class Pair:
    """Represents a pair of genomic intervals."""

    def __init__(self, row, chrA, startA, endA, strandA, chrB, startB, endB, strandB, overlap_len_col=None):
        """
        Initializes a Pair object from a pandas Series (row).

        Args:
            row (pandas.Series): A row from the DataFrame representing the pair.
            chrA, startA, endA, strandA, chrB, startB, endB, strandB (str): Column names.
            overlap_len_col (str, optional): Column name for the overlap length.  Defaults to None.
        """
        self.chrA = row[chrA]
        self.startA = row[startA]
        self.endA = row[endA]
        self.strandA = row[strandA]
        self.chrB = row[chrB]
        self.startB = row[startB]
        self.endB = row[endB]
        self.strandB = row[strandB]
        self.overlap_len = row[overlap_len_col] if overlap_len_col else None
        self.row_index = row.name  # Store the original row index

    def calculate_interval_overlaps(self, other, fraction_select=None):
        """
        Calculates the overlap between this pair and another Pair.

        Args:
            other (Pair): The other Pair object to calculate overlap with.
            fraction_select (int, optional): If 0, returns overlap fraction relative to this Pair.
                                            If 1, returns overlap fraction relative to the other Pair.
                                            If None (default), returns the maximum overlap fraction.

        Returns:
            float: The overlap fraction (0 if no overlap).
        """
        overlap_span_A = max(0, min(self.endA, other.endA) - max(self.startA, other.startA))
        overlap_span_B = max(0, min(self.endB, other.endB) - max(self.startB, other.startB))

        if overlap_span_A == 0 or overlap_span_B == 0:
            return 0

        size1_A = self.endA - self.startA
        size2_A = other.endA - other.startA
        size1_B = self.endB - self.startB
        size2_B = other.endB - other.startB

        if fraction_select is not None:
            if fraction_select == 0:
                return max(overlap_span_A / size1_A, overlap_span_B / size1_B)
            elif fraction_select == 1:
                return max(overlap_span_A / size2_A, overlap_span_B / size2_B)
            else:
                raise ValueError("fraction_select must be 0, 1, or None.")
        else:  # Return the *maximum* overlap fraction
            return max(overlap_span_A / size1_A, overlap_span_B / size1_B,
                       overlap_span_A / size2_A, overlap_span_B / size2_B)

    def is_same_pair(self, other):
        """
        Checks if this pair is the same as another pair (or its swapped version).

        Args:
            other (Pair): The other Pair object to compare with.

        Returns:
            bool: True if the pairs are the same, False otherwise.
        """
        # Check for both original and swapped order
        return (
            (self.chrA == other.chrA and self.startA == other.startA and
             self.endA == other.endA and self.strandA == other.strandA and
             self.chrB == other.chrB and self.startB == other.startB and
             self.endB == other.endB and self.strandB == other.strandB)
            or
            (self.chrA == other.chrB and self.startA == other.startB and
             self.endA == other.endB and self.strandA == other.strandB and
             self.chrB == other.chrA and self.startB == other.startA and
             self.endB == other.endA and self.strandB == other.strandA)
        )

    def is_umbrella_pair(self, other, coverage_threshold):
        """
        Checks if this pair is an umbrella pair enclosing another pair.

        Args:
            other (Pair): The other Pair object to compare with.
            coverage_threshold (float): The minimum overlap fraction required.

        Returns:
            bool: True if this pair is an umbrella pair, False otherwise.
        """
        if self.chrA != other.chrA or self.chrB != other.chrB:
            return False

        # Check for strand consistency
        if (self.strandA == other.strandA and self.strandB != other.strandB) or \
           (self.strandA != other.strandA and self.strandB == other.strandB):
            return False

        # Calculate overlap fractions
        coverage_A = self.calculate_interval_overlaps(other)
        if coverage_A < coverage_threshold:
            return False

        coverage_B = other.calculate_interval_overlaps(self) #Correct implementation
        if coverage_B < coverage_threshold:
            return False

        # Check if this pair's overlap length is greater than or equal to the other's
        return self.overlap_len >= other.overlap_len

    def __repr__(self):
        """Returns a string representation of the Pair."""
        return (f"Pair({self.chrA}:{self.startA}-{self.endA}:{self.strandA}, "
                f"{self.chrB}:{self.startB}-{self.endB}:{self.strandB})")

def _find_umbrella_pairs(groupdf, coverage_threshold, overlap_len_col):
    """
    Internal function to find umbrella pairs within a grouped DataFrame.

    Args:
        groupdf (pd.DataFrame): Grouped DataFrame.
        coverage_threshold (float): Coverage threshold.
        overlap_len_col (str): Column name for overlap length.

    Returns:
        list: List of indices of umbrella pairs.
    """
    umbrella_indices = []
    chrA, startA, endA, strandA, chrB, startB, endB, strandB = groupdf.columns[:8]
    # Create Pair objects for each row in the group
    pairs = [Pair(row, chrA, startA, endA, strandA, chrB, startB, endB, strandB, overlap_len_col)
             for _, row in groupdf.iterrows()]

    # Compare each pair to every other pair
    for i in range(len(pairs) - 1):
        for j in range(i + 1, len(pairs)):
            if pairs[i].is_umbrella_pair(pairs[j], coverage_threshold):
                umbrella_indices.append(pairs[i].row_index)  # Store the *original index*

    if umbrella_indices:
        logger.info(f"Found {len(umbrella_indices)} umbrella pairs.")
    return umbrella_indices


def filter_umbrella_pairs(groupdf, coverage_threshold=0.9, overlap_len_col="overlap_len"):
    """
    Filters out umbrella pairs from a grouped DataFrame.

    Args:
        groupdf (pd.DataFrame): Grouped DataFrame.
        coverage_threshold (float): Coverage threshold.  Defaults to 0.9.
        overlap_len_col (str): Column name for overlap length. Defaults to "overlap_len".

    Returns:
        pd.DataFrame: Filtered DataFrame with umbrella pairs removed.
    """
    if overlap_len_col not in groupdf.columns:
        raise ValueError(f"Overlap length column '{overlap_len_col}' not found in DataFrame.")

    # 1. Preprocessing
    chrA, startA, endA, strandA, chrB, startB, endB, strandB = groupdf.columns[:8]
    groupdf.drop_duplicates(inplace=True)
    groupdf.reset_index(drop=True, inplace=True)  # Reset index after dropping duplicates
    groupdf["size"] = (groupdf[endA] - groupdf[startA]) + (groupdf[endB] - groupdf[startB])
    groupdf.sort_values(by="size", ascending=False, inplace=True) # Sort by the size of the pairs, so that the larger pairs are more likely to be kept.
    groupdf.drop(columns=["size"], inplace=True) # Remove the temp "size" column
    groupdf.reset_index(inplace=True) # Keep track of original index before dropping rows.

    # 2. Find umbrella pairs (using the internal function)
    umbrella_indices = _find_umbrella_pairs(groupdf, coverage_threshold, overlap_len_col)

    # 3. Filter the DataFrame
    filtered_df = groupdf.drop(index=umbrella_indices) # Drop based on original indices.
    filtered_df.reset_index(drop=True, inplace=True) # Reset index after drop

    if len(filtered_df) < len(groupdf):
        logger.info(f"Removed {len(groupdf) - len(filtered_df)} umbrella pairs.")
    return filtered_df
