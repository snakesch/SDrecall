import sys
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
            fraction_select (int, optional): If self, returns overlap fraction relative to this Pair.
                                             If other, returns overlap fraction relative to the other Pair.
                                             If None (default), returns the maximum overlap fraction.

        Returns:
            float: The overlap fraction (0 if no overlap).
        """
        overlap_span_A = max(0, min(self.endA, other.endA) - max(self.startA, other.startA))
        overlap_span_B = max(0, min(self.endB, other.endB) - max(self.startB, other.startB))

        size1_A = self.endA - self.startA
        size2_A = other.endA - other.startA
        size1_B = self.endB - self.startB
        size2_B = other.endB - other.startB

        fraction_1 = (overlap_span_A / size1_A, overlap_span_B / size1_B)
        fraction_2 = (overlap_span_A / size2_A, overlap_span_B / size2_B)

        if fraction_select is not None:
            if fraction_select == "self":
                return fraction_1
            elif fraction_select == "other":
                return fraction_2
            else:
                raise ValueError("fraction_select must be 'self' or 'other' or None.")
        else:  # Return the *maximum* overlap fraction
            return max(fraction_1, fraction_2, key=lambda x: sum(x))
        

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
            (self.chrA == other.chrA and \
             self.startA == other.startA and \
             self.endA == other.endA and \
             self.chrB == other.chrB and \
             self.startB == other.startB and \
             self.endB == other.endB and \
             (self.strandB == self.strandA) == (other.strandB == other.strandA))
            or
            (self.chrA == other.chrB and \
             self.startA == other.startB and \
             self.endA == other.endB and \
             self.chrB == other.chrA and \
             self.startB == other.startA and \
             self.endB == other.endA and \
             (self.strandB == self.strandA) == (other.strandB == other.strandA))
        )

    def is_umbrella_pair(self, other, coverage_threshold=0.95):
        """
        Checks if this pair is an umbrella pair enclosing another pair.
        Input order matters for the same pair of SDs.

        Args:
            other (Pair): The other Pair object to compare with.
            coverage_threshold (float): The minimum overlap fraction required.

        Returns:
            bool: True if this pair is an umbrella pair, False otherwise.
        """
        if self.chrA != other.chrA or self.chrB != other.chrB:
            return False

        # Check for strand consistency
        if (self.strandA == self.strandB) != (other.strandA == other.strandB):
            return False

        # Calculate overlap fractions
        coverage_fraction_A, coverage_fraction_B = self.calculate_interval_overlaps(other, fraction_select="other")
        # print(f"Coverage fractions: {coverage_fraction_A}, {coverage_fraction_B} between {self} and {other}", file=sys.stderr)
        if coverage_fraction_A < coverage_threshold or coverage_fraction_B < coverage_threshold:
            return False

        # Check if this pair's overlap length is smaller than or equal to the other's
        return self.overlap_len <= other.overlap_len


    def __repr__(self):
        """Returns a string representation of the Pair."""
        return (f"Pair({self.chrA}:{self.startA}-{self.endA}:{self.strandA}, {self.chrB}:{self.startB}-{self.endB}:{self.strandB})")

def _find_umbrella_pairs(groupdf, overlap_len_col, coverage_threshold=0.95):
    """
    Internal function to find umbrella pairs within a grouped DataFrame.
    Umbrella pairs (more general pairs that encompass other pairs) are identified for removal.

    Args:
        groupdf (pd.DataFrame): Grouped DataFrame.
        overlap_len_col (str): Column name for overlap length.
        coverage_threshold (float): Coverage threshold.

    Returns:
        list: List of indices of umbrella pairs to be removed.
    """
    to_remove = set()  # Pairs marked for removal
    chrA, startA, endA, strandA, chrB, startB, endB, strandB = groupdf.columns[:8]
    
    # Create Pair objects for each row in the group
    pairs = [Pair(row, chrA, startA, endA, strandA, chrB, startB, endB, strandB, overlap_len_col)
             for _, row in groupdf.iterrows()]

    # print(f"Pairs: {pairs}\n", file=sys.stderr)
    
    # Compare each pair to every other pair
    for i in range(len(pairs)):
        if i in to_remove:
            continue  # Skip pairs already marked for removal
        
        for j in range(i + 1, len(pairs)):
            if j in to_remove:
                continue  # Skip pairs already marked for removal
            
            # Check umbrella relationship (if one pair encompasses the other)
            if pairs[i].is_umbrella_pair(pairs[j], coverage_threshold):
                # print(f"Comparing {pairs[i]} to {pairs[j]}", file=sys.stderr)
                # If i is an umbrella pair for j, remove i and stop comparing i
                to_remove.add(i)
                break  # No need to compare i with other pairs
            elif pairs[j].is_umbrella_pair(pairs[i], coverage_threshold):
                # print(f"Comparing {pairs[j]} to {pairs[i]}", file=sys.stderr)
                # If j is an umbrella pair for i, remove j and continue with i
                to_remove.add(j)
                # Continue comparing i with other pairs
    
    # Convert to list of original row indices
    return [pairs[i].row_index for i in to_remove]


def filter_umbrella_pairs(groupdf, coverage_threshold=0.95, overlap_len_col="overlap_len"):
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

    # 2. Find umbrella pairs (using the internal function)
    # print(f"Grouped DataFrame: \n{groupdf.head(10).to_string(index=False)}", file=sys.stderr)
    umbrella_indices = _find_umbrella_pairs(groupdf, overlap_len_col, coverage_threshold)

    # 3. Filter the DataFrame
    filtered_df = groupdf.drop(index=umbrella_indices) # Drop based on original indices.
    filtered_df.reset_index(drop=True, inplace=True) # Reset index after drop

    if len(filtered_df) < len(groupdf):
        logger.debug(f"Removed {len(groupdf) - len(filtered_df)} umbrella pairs.")
    return filtered_df
