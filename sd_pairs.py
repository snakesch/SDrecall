import pandas
import logging
logger = logging.getLogger("SDrecall")

'''
A pair of SD is considered an umbrella pair if it completely encloses another pair of SDs.
'''

def calculate_interval_overlaps(interval1, interval2, fraction_select = None):
    overlap_span = min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])
    if overlap_span < 0:
        return 0
    else:
        interval1_size = interval1[1] - interval1[0]
        interval2_size = interval2[1] - interval2[0]
        if fraction_select:
            return [overlap_span/interval1_size, overlap_span/interval2_size][fraction_select]
        else:
            return max(overlap_span/interval1_size, overlap_span/interval2_size)

def is_same_pair(row, pair,
                 chrA, startA, endA, strandA, chrB, startB, endB, strandB):
    """
    Check if a row represents the same pair as the given pair.

    Args:
        row (pandas.Series): Row from the dataframe.
        pair (tuple): Pair represented as (chrA, startA, endA, strandA, chrB, startB, endB, strandB).
        chr*, start*, end*, strand*: Column names.

    Returns:
        bool: True if the row represents the same pair, False otherwise.
    """
    return (
        row[chrA] == pair[0] and
        row[startA] == pair[1] and
        row[endA] == pair[2] and
        row[strandA] == pair[3] and
        row[chrB] == pair[4] and
        row[startB] == pair[5] and
        row[endB] == pair[6] and
        row[strandB] == pair[7]
    )

def find_umbrella_pairs(groupdf, coverage_threshold):
    """
    Find the umbrella pairs within a grouped dataframe.

    Args:
        groupdf (pandas.DataFrame): Grouped dataframe containing SD pairs.
        coverage_threshold (float): Minimum coverage threshold for umbrella pairs.

    Returns:
        list: List of umbrella pairs represented as tuples of (chrA, startA, endA, strandA, chrB, startB, endB, strandB).
    """
    umbrella_pairs = []
    chrA, startA, endA, strandA, chrB, startB, endB, strandB = groupdf.columns[:8]

    for i, row1 in groupdf.iloc[:-1].iterrows():
        for j, row2 in groupdf.iloc[i+1:].iterrows():
            if is_umbrella_pair(row1, row2, coverage_threshold, row1["overlap_len"], row2["overlap_len"],
                                chrA, startA, endA, strandA, chrB, startB, endB, strandB):
                # If the condition returns True, then we consider pair1 (row1) is an umbrella pair enclosing pair2 (row2), so we can rule out pair1 record from the group df.
                umbrella_pairs.append((
                    row1[chrA], row1[startA], row1[endA], row1[strandA],
                    row1[chrB], row1[startB], row1[endB], row1[strandB]
                ))

    if len(umbrella_pairs) > 0:
        logger.info(f"Found {len(umbrella_pairs)} umbrella pairs: {umbrella_pairs}")

    return umbrella_pairs


def filter_umbrella_pairs(groupdf, coverage_threshold=0.9):
    """
    Filter out rows containing umbrella pairs from the grouped dataframe.
    If a pair of SDs is an umbrella to other pairs of SDs, then the umbrella pair might be redundant to be considered for further analysis.

    Remember that the groupdf contains a list of SD pairs where the first interval in each pair is overlapping with the same target region interval.

    Args:
        groupdf (pandas.DataFrame): Grouped dataframe containing SD pairs.
        coverage_threshold (float, optional): Minimum coverage threshold for umbrella pairs (default: 0.9).

    Returns:
        pandas.DataFrame: Filtered dataframe with rows containing umbrella pairs removed.
    """
    # Need to sort the groupdf by interval pair total sizes.
    chrA, startA, endA, strandA, chrB, startB, endB, strandB = groupdf.columns[:8]
    groupdf.drop_duplicates(inplace=True)
    groupdf["size"] = groupdf.apply(lambda row: (row[endA] - row[startA]) + (row[endB] - row[startB]), axis=1)
    groupdf = groupdf.sort_values(by="size", ascending=False).drop(columns=["size"])

    umbrella_pairs = find_umbrella_pairs(groupdf, coverage_threshold)

    filtered_df = groupdf[
        ~groupdf.apply(lambda row: any(
            is_same_pair(row, pair, chrA, startA, endA, strandA, chrB, startB, endB, strandB) or is_same_pair(row, (pair[4:] + pair[:4]), chrA, startA, endA, strandA, chrB, startB, endB, strandB)
            for pair in umbrella_pairs
        ), axis=1)
    ]

    if len(filtered_df) < len(groupdf):
        logger.info(f"Removed {len(groupdf) - len(filtered_df)} umbrella pairs from the grouped dataframe. Before filtering, the table looks like this: \n{groupdf.to_string(index=False)}\n Now the table looks like this:\n{filtered_df.to_string(index=False)}\n")

    # logger.info(f"Removed {len(groupdf) - len(filtered_df)} umbrella pairs from the grouped dataframe. Before filtering, the table looks like this: \n{groupdf.to_string(index=False)}\n Now the table looks like this:\n{filtered_df.to_string(index=False)}\n")
    return filtered_df



def is_umbrella_pair(pair1, pair2, coverage_threshold, overlap1, overlap2,
                     chrA, startA, endA, strandA, chrB, startB, endB, strandB):
    """
    Basically the input pairs follow the rules below:
    the first interval in both pairs are overlapping with the same target region interval.
    We want to identify whether the pair of SDs in pair1 is almost enclosing the pair of SDs in pair2.
    If pair1 is an umbrella pair of pair2, then pair1 might be redundant to be considered for further analysis.

    Args:
        pair1 (pandas.Series): First pair.
        pair2 (pandas.Series): Second pair.
        coverage_threshold (float): Minimum coverage threshold for umbrella pairs.

    Returns:
        bool: True if the pairs form an umbrella pair, False otherwise.
    """
    if pair1[chrA] != pair2[chrA] or pair1[chrB] != pair2[chrB]:
        return False

    if pair1[strandA] == pair2[strandA] and pair1[strandB] != pair2[strandB]:
        return False

    if pair1[strandA] != pair2[strandA] and pair1[strandB] == pair2[strandB]:
        return False

    # The returned coverageA is actually the overlapping coefficient between the first intervals in pair1 and pair2
    coverage_A = calculate_interval_overlaps(
        (pair1[startA], pair1[endA]),
        (pair2[startA], pair2[endA])
    )

    # The returned coverageB is actually the overlapping coefficient between the second intervals in pair1 and pair2
    coverage_B = calculate_interval_overlaps(
        (pair1[startB], pair1[endB]),
        (pair2[startB], pair2[endB])
    )

    # If the first intervals between pair1 and pair2 are largely overlapping and the second intervals are also largely overlapping,
    # And the two intervals in pair1 are larger than the pair of intervals in pair2,
    # Then we consider pair1 is an "umbrella pair" enclosing pair2.
    result = coverage_A >= coverage_threshold and coverage_B >= coverage_threshold and overlap1 >= overlap2

    if result:
        logger.info(f"Found an umbrella pair: {pair1[chrA]}:{pair1[startA]}-{pair1[endA]}-{pair1[strandA]} and {pair1[chrB]}:{pair1[startB]}-{pair1[endB]}-{pair1[strandB]}. \n This is by comparison with {pair2[chrA]}:{pair2[startA]}-{pair2[endA]}-{pair2[strandA]} and {pair2[chrB]}:{pair2[startB]}-{pair2[endB]}-{pair2[strandB]}")

    return result
