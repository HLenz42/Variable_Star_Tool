# Uses a star as a comparison star only if it is tracked through 100% of the images and has the flag is_variable='no'
# If the target star's magnitude is NaN in an image, the differential magnitude will be output as NaN

import pandas as pd
import numpy as np
from scipy.signal import savgol_filter  # Used for potential smoothing (not used in this version)
import sys
import os
import heapq  # For efficient nearest-neighbor selection in KDTree
import json
from tqdm import tqdm  # For progress bar during processing

# === KDTree CLASS ===
# KDTree is used to efficiently find the nearest comparison stars based on sky coordinates (RA, DEC)
class KDTree(object):
    def __init__(self, points, dim, dist_sq_func=None):
        # Initialize KDTree with a list of points (RA, DEC, star_index), dimension (2 for RA/DEC), and optional distance function
        if dist_sq_func is None:
            # Default distance function: Euclidean distance squared in RA/DEC space
            dist_sq_func = lambda a, b: sum((x - b[i]) ** 2 for i, x in enumerate(a))

        def make(points, i=0):
            # Recursive function to build the KDTree
            if len(points) > 1:
                # Sort points by the current dimension (alternates between RA and DEC)
                points.sort(key=lambda x: x[i])
                i = (i + 1) % dim  # Switch to next dimension for child nodes
                m = len(points) >> 1  # Find median index
                # Create node with left subtree, right subtree, and median point
                return [make(points[:m], i), make(points[m + 1:], i), points[m]]
            if len(points) == 1:
                # Leaf node with single point
                return [None, None, points[0]]
            return None  # Empty node

        def get_knn(node, point, k, return_dist_sq, heap, i=0, tiebreaker=1):
            # Recursive function to find k nearest neighbors
            if node is not None:
                # Calculate distance squared from query point to current node's point
                dist_sq = dist_sq_func(point, node[2])
                dx = node[2][i] - point[i]  # Distance along current dimension
                # Add to heap if fewer than k points or if closer than farthest in heap
                if len(heap) < k:
                    heapq.heappush(heap, (-dist_sq, tiebreaker, node[2]))
                elif dist_sq < -heap[0][0]:
                    heapq.heappushpop(heap, (-dist_sq, tiebreaker, node[2]))
                i = (i + 1) % dim  # Switch to next dimension
                # Explore subtrees based on distance to splitting plane
                for b in (dx < 0, dx >= 0)[:1 + (dx * dx < -heap[0][0])]:
                    get_knn(node[int(b)], point, k, return_dist_sq, heap, i, (tiebreaker << 1) | b)
            if tiebreaker == 1:
                # Return sorted list of k nearest points (with distances if requested)
                return [(-h[0], h[2]) if return_dist_sq else h[2] for h in sorted(heap)][::-1]

        self._get_knn = get_knn
        self._root = make(points)  # Build the tree from input points

    def get_knn(self, point, k, return_dist_sq=True):
        # Public method to get k nearest neighbors for a given point
        return self._get_knn(self._root, point, k, return_dist_sq, [])

# === User Input for Magnitude Type ===
# Prompt user to choose between MAG_AUTO (automatic aperture) or MAG_APER (fixed aperture)
while True:
    print(">>> Use MAG_AUTO or MAG_APER?")
    print("(1) MAG_AUTO")
    print("(2) MAG_APER\n")

    which_mag = int(input("Answer: "))
    print()

    if which_mag == 1:
        MAG_COLUMN = 'MAG_AUTO'
        MAG_ERROR_COLUMN = 'MAGERR_AUTO'
        break
    if which_mag == 2:
        MAG_COLUMN = 'MAG_APER'
        MAG_ERROR_COLUMN = 'MAGERR_APER'
        break
    else:
        print("Error: Please select either 1, 2" + "\n")

# Display selected magnitude columns
print(">>> Magnitude option selection")
print(f"    MAG_COLUMN: {MAG_COLUMN}")
print(f"    MAG_ERROR_COLUMN: {MAG_ERROR_COLUMN}\n")

# === Configuration ===
# Define column names and parameters used throughout the script
TIME_COLUMN = 'julian_date'  # Column for observation timestamps
RA_COLUMN = 'ALPHA_J2000'   # Right Ascension (degrees)
DEC_COLUMN = 'DELTA_J2000'  # Declination (degrees)
STAR_NAME_COLUMN = 'Simbad_Name'  # Star identifier from Simbad
COMPARISON_THRESHOLD_FACTOR = 1.5  # Not used in this version, but reserved for brightness filtering
SG_WINDOW_LENGTH = 20  # Savitzky-Golay filter window (not used)
SG_POLYNOMIAL_ORDER = 2  # Savitzky-Golay polynomial order (not used)

# === Main Script ===
if __name__ == "__main__":
    # Parse command-line arguments
    target = sys.argv[3]  # Target star identifier
    date = sys.argv[4]    # Observation date
    csv_name = sys.argv[5]  # Input CSV filename
    input_path = os.path.join(sys.argv[1], csv_name)  # Full path to input CSV
    output_dir = sys.argv[2]  # Directory for output files
    output_csv = f"tracked_stars_with_differential_mags_{target}_{date}.csv"  # Output CSV filename
    output_json = f"comparison_star_map_{target}_{date}.json"  # Output JSON filename for comparison star map

    # Load input data
    print("Loading data...")
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        sys.exit(1)

    # Ensure star_index is integer and sort by time and star index
    df['star_index'] = df['star_index'].astype(int)
    df.sort_values(by=[TIME_COLUMN, 'star_index'], inplace=True)

    # Calculate tracking percentage for each star
    print("Calculating tracking percentage...")
    total_images = df[TIME_COLUMN].nunique()  # Number of unique observation times
    # Count non-NaN measurements for each star (must have valid mag, RA, DEC)
    non_nan_counts = df.groupby('star_index').apply(
        lambda x: x[[MAG_COLUMN, RA_COLUMN, DEC_COLUMN]].notna().all(axis=1).sum(),
        include_groups=False
    )
    tracking_percentage = (non_nan_counts / total_images) * 100  # Percentage of images where star is tracked

    # Input number of comparison stars
    while True:
        try:
            num_comps = int(input("Enter number of comparison stars: "))
            if num_comps > 0:
                break
        except ValueError:
            print("Please enter a positive integer.")
            continue

    # Pivot data for efficient access: rows are times, columns are stars
    print("Pivoting data...")
    pivot_mags = df.pivot_table(index=TIME_COLUMN, columns='star_index', values=MAG_COLUMN)
    pivot_errs = df.pivot_table(index=TIME_COLUMN, columns='star_index', values=MAG_ERROR_COLUMN)

    # Extract unique metadata for each star (e.g., RA, DEC, is_variable)
    unique_meta = df.drop_duplicates(subset=['star_index']).set_index('star_index')
    unique_meta['tracking_percentage'] = tracking_percentage

    # Select comparison stars: must have is_variable='no' and 100% tracking
    print("Selecting comparison stars...")
    if 'is_variable' not in unique_meta.columns:
        print("Error: 'is_variable' column not found in data.")
        sys.exit(1)
    pot_comps = unique_meta[
        (unique_meta['is_variable'].str.lower() == 'no') &  # Non-variable stars
        (unique_meta['tracking_percentage'] == 100)        # Tracked in all images
    ].copy()

    print(f"Found {len(pot_comps)} potential comparison stars (is_variable='no', 100% tracking).")
    if len(pot_comps) < num_comps:
        print(f"Warning: Fewer comparison stars ({len(pot_comps)}) than requested ({num_comps}). Using available stars.")

    # Build KDTree for comparison stars using RA, DEC coordinates
    tree_coords = list({(unique_meta.loc[i, RA_COLUMN], unique_meta.loc[i, DEC_COLUMN], i) for i in pot_comps.index})
    if not tree_coords:
        print("Error: No valid comparison stars found.")
        sys.exit(1)
    comp_tree = KDTree(tree_coords, 2)

    # Initialize dictionaries to store results for each target star
    selected_comps_map = {}  # Maps target star index to list of comparison star indices
    unw_diff_mags = {}      # Unweighted differential magnitudes
    unw_diff_errs = {}      # Unweighted differential magnitude errors
    w_diff_mags = {}        # Weighted differential magnitudes
    w_diff_errs = {}        # Weighted differential magnitude errors

    print("\nStarting differential photometry loop...")
    for target_idx in tqdm(unique_meta.index, desc="Processing targets"):
        # Get target star's coordinates
        target_coords = (unique_meta.loc[target_idx, RA_COLUMN], unique_meta.loc[target_idx, DEC_COLUMN])
        # Find nearest comparison stars (request extra to account for self-exclusion)
        neighbors = comp_tree.get_knn(target_coords, num_comps + 10, False)

        # Select up to num_comps unique comparison stars, excluding the target itself
        seen = set()
        comps = []
        for ra, dec, idx in neighbors:
            if idx != target_idx and idx not in seen:
                comps.append(idx)
                seen.add(idx)
            if len(comps) >= num_comps:
                break

        selected_comps_map[target_idx] = comps

        # If no comparison stars are found, output NaN for all results
        if not comps:
            unw_diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
            unw_diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
            w_diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
            w_diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
            continue

        # Get target magnitudes and errors
        t_mag = pivot_mags[target_idx]
        t_err = pivot_errs[target_idx]
        diffs = []    # List of differential magnitudes (m_T - m_n) for each comparison star
        vars_list = []  # List of variances (σ_T^2 + σ_n^2) for each comparison star

        # Calculate differential magnitudes and variances for each comparison star
        for c_idx in comps:
            c_mag = pivot_mags[c_idx]
            c_err = pivot_errs[c_idx]
            # Align target and comparison data to ensure matching timestamps
            t_aligned, c_aligned = t_mag.align(c_mag, join='inner')
            te_aligned, ce_aligned = t_err.align(c_err, join='inner')
            if not t_aligned.empty:
                # Compute differences and variances only where both mags are non-NaN
                valid_mask = t_aligned.notna() & c_aligned.notna()
                diff = pd.Series(np.nan, index=t_aligned.index)
                diff[valid_mask] = t_aligned[valid_mask] - c_aligned[valid_mask]  # m_T - m_n
                var = pd.Series(np.nan, index=te_aligned.index)
                var[valid_mask] = np.square(te_aligned[valid_mask]) + np.square(ce_aligned[valid_mask])  # σ_T^2 + σ_n^2
                diffs.append(diff)
                vars_list.append(var)

        # If no valid differences, output NaN for all results
        if not diffs:
            unw_diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
            unw_diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
            w_diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
            w_diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
            continue

        # Convert lists to DataFrames: rows are times, columns are comparison stars
        df_diffs = pd.DataFrame(diffs).T
        df_vars = pd.DataFrame(vars_list).T

        # Unweighted differential photometry
        # m_diff = (1/N) * Σ (m_T - m_n)
        unw_diff = df_diffs.mean(axis=1, skipna=True)  # Mean of differences across comparison stars
        num_valid = df_diffs.notna().sum(axis=1)  # Number of valid comps per time
        # σ_diff = sqrt(Σ (σ_T^2 + σ_n^2)) / N
        unw_err = np.sqrt(df_vars.sum(axis=1, skipna=True)) / num_valid

        # Weighted differential photometry
        # w_n = 1 / (σ_T^2 + σ_n^2)
        df_weights = 1 / df_vars  # Weights are inverse variances (NaN where var is NaN or 0)
        # m_diff = Σ (w_n * (m_T - m_n)) / Σ w_n
        weighted_sum = (df_diffs * df_weights).sum(axis=1, skipna=True)
        total_weight = df_weights.sum(axis=1, skipna=True)
        w_diff = weighted_sum / total_weight  # Weighted mean of differences
        # σ_diff = 1 / sqrt(Σ w_n)
        w_err = np.sqrt(1 / total_weight)

        # Reindex to full time series to ensure all timestamps are included
        unw_diff_mags[target_idx] = unw_diff.reindex(pivot_mags.index)
        unw_diff_errs[target_idx] = unw_err.reindex(pivot_errs.index)
        w_diff_mags[target_idx] = w_diff.reindex(pivot_mags.index)
        w_diff_errs[target_idx] = w_err.reindex(pivot_errs.index)

    print("\nDifferential photometry complete. Mapping and saving results...")

    # Convert result dictionaries to DataFrames
    df_unw_mags = pd.DataFrame(unw_diff_mags)
    df_unw_errs = pd.DataFrame(unw_diff_errs)
    df_w_mags = pd.DataFrame(w_diff_mags)
    df_w_errs = pd.DataFrame(w_diff_errs)

    # Reshape to long format for merging with original DataFrame
    long_unw_mags = df_unw_mags.stack().rename('differential_mag').reset_index()
    long_unw_mags.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)
    long_unw_errs = df_unw_errs.stack().rename('differential_mag_err').reset_index()
    long_unw_errs.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)
    long_w_mags = df_w_mags.stack().rename('weighted_differential_mag').reset_index()
    long_w_mags.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)
    long_w_errs = df_w_errs.stack().rename('weighted_differential_mag_err').reset_index()
    long_w_errs.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)

    # Merge results into original DataFrame
    df = df.merge(long_unw_mags, on=['star_index', TIME_COLUMN], how='left')
    df = df.merge(long_unw_errs, on=['star_index', TIME_COLUMN], how='left')
    df = df.merge(long_w_mags, on=['star_index', TIME_COLUMN], how='left')
    df = df.merge(long_w_errs, on=['star_index', TIME_COLUMN], how='left')

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, output_csv), index=False, mode = 'w')

    # Save comparison star map as JSON
    with open(os.path.join(output_dir, output_json), 'w') as f:
        json.dump(selected_comps_map, f, indent=4)

    print(f"\nSaved results to {output_csv} and {output_json}")

# # Uses a star as a comparison star only if it is tracked though 100% of the images, and has the flag is_variable=no
# # If the target star itself is a NaN then the diff mag will be outputted as a NaN

# import pandas as pd
# import numpy as np
# from scipy.signal import savgol_filter
# import sys
# import os
# import heapq
# import json
# from tqdm import tqdm  

# # === KDTree CLASS ===
# class KDTree(object):
#     def __init__(self, points, dim, dist_sq_func=None):
#         if dist_sq_func is None:
#             dist_sq_func = lambda a, b: sum((x - b[i]) ** 2 for i, x in enumerate(a))

#         def make(points, i=0):
#             if len(points) > 1:
#                 points.sort(key=lambda x: x[i])
#                 i = (i + 1) % dim
#                 m = len(points) >> 1
#                 return [make(points[:m], i), make(points[m + 1:], i), points[m]]
#             if len(points) == 1:
#                 return [None, None, points[0]]
#             return None

#         def get_knn(node, point, k, return_dist_sq, heap, i=0, tiebreaker=1):
#             if node is not None:
#                 dist_sq = dist_sq_func(point, node[2])
#                 dx = node[2][i] - point[i]
#                 if len(heap) < k:
#                     heapq.heappush(heap, (-dist_sq, tiebreaker, node[2]))
#                 elif dist_sq < -heap[0][0]:
#                     heapq.heappushpop(heap, (-dist_sq, tiebreaker, node[2]))
#                 i = (i + 1) % dim
#                 for b in (dx < 0, dx >= 0)[:1 + (dx * dx < -heap[0][0])]:
#                     get_knn(node[int(b)], point, k, return_dist_sq, heap, i, (tiebreaker << 1) | b)
#             if tiebreaker == 1:
#                 return [(-h[0], h[2]) if return_dist_sq else h[2] for h in sorted(heap)][::-1]

#         self._get_knn = get_knn
#         self._root = make(points)

#     def get_knn(self, point, k, return_dist_sq=True):
#         return self._get_knn(self._root, point, k, return_dist_sq, [])

# while True:
#     print(">>> Use MAG_AUTO or MAG_APER?")
#     print("(1) MAG_AUTO")
#     print("(2) MAG_APER\n")

#     which_mag = int(input("Answer: "))
#     print()

#     if which_mag == 1:
#         MAG_COLUMN = 'MAG_AUTO'
#         MAG_ERROR_COLUMN = 'MAGERR_AUTO'
#         break

#     if which_mag == 2:
#         MAG_COLUMN = 'MAG_APER'
#         MAG_ERROR_COLUMN = 'MAGERR_APER'
#         break

#     else:
#         print("Error: Please select either 1, 2" + "\n")

# print(">>> Magnitude option selection")
# print(f"    MAG_COLUMN: {MAG_COLUMN}")
# print(f"    MAG_ERROR_COLUMN: {MAG_ERROR_COLUMN}\n")


# # === Configuration ===

# TIME_COLUMN = 'julian_date'
# RA_COLUMN = 'ALPHA_J2000'
# DEC_COLUMN = 'DELTA_J2000'
# STAR_NAME_COLUMN = 'Simbad_Name'
# COMPARISON_THRESHOLD_FACTOR = 1.5  # For potential brightness filter
# SG_WINDOW_LENGTH = 20
# SG_POLYNOMIAL_ORDER = 2

# # === Main Script ===
# if __name__ == "__main__":
#     target = sys.argv[3] 
#     date = sys.argv[4]
#     csv_name = sys.argv[5]
#     input_path = os.path.join(sys.argv[1], csv_name)
#     output_dir = sys.argv[2]
#     output_csv = f"tracked_stars_with_differential_mags_{target}_{date}.csv"
#     output_json = f"comparison_star_map_{target}_{date}.json"

#     print("Loading data...")
#     try:
#         df = pd.read_csv(input_path)
#     except FileNotFoundError:
#         print(f"Error: Input file not found at {input_path}")
#         sys.exit(1)

#     df['star_index'] = df['star_index'].astype(int)
#     df.sort_values(by=[TIME_COLUMN, 'star_index'], inplace=True)

#     # Calculate tracking percentage
#     print("Calculating tracking percentage...")
#     total_images = df[TIME_COLUMN].nunique()
#     non_nan_counts = df.groupby('star_index').apply(
#         lambda x: x[[MAG_COLUMN, RA_COLUMN, DEC_COLUMN]].notna().all(axis=1).sum(),
#         include_groups=False
#     )
#     tracking_percentage = (non_nan_counts / total_images) * 100

#     # Input number of comparison stars
#     while True:
#         try:
#             num_comps = int(input("Enter number of comparison stars: "))
#             if num_comps > 0:
#                 break
#         except ValueError:
#             print("Please enter a positive integer.")
#             continue

#     print("Pivoting data...")
#     pivot_mags = df.pivot_table(index=TIME_COLUMN, columns='star_index', values=MAG_COLUMN)
#     pivot_errs = df.pivot_table(index=TIME_COLUMN, columns='star_index', values=MAG_ERROR_COLUMN)

#     unique_meta = df.drop_duplicates(subset=['star_index']).set_index('star_index')
#     unique_meta['tracking_percentage'] = tracking_percentage

#     # Select comparison stars: no NaNs (100% tracking) and is_variable == 'no'
#     print("Selecting comparison stars...")
#     if 'is_variable' not in unique_meta.columns:
#         print("Error: 'is_variable' column not found in data.")
#         sys.exit(1)
#     pot_comps = unique_meta[
#         (unique_meta['is_variable'].str.lower() == 'no') &     #CONDITIONS FOR A COMPARISON STAR
#         (unique_meta['tracking_percentage'] == 100)
#     ].copy()

#     print(f"Found {len(pot_comps)} potential comparison stars (is_variable='no', 100% tracking).")
#     if len(pot_comps) < num_comps:
#         print(f"Warning: Fewer comparison stars ({len(pot_comps)}) than requested ({num_comps}). Using available stars.")

#     # Build KD-tree for comparison stars
#     tree_coords = list({(unique_meta.loc[i, RA_COLUMN], unique_meta.loc[i, DEC_COLUMN], i) for i in pot_comps.index})
#     if not tree_coords:
#         print("Error: No valid comparison stars found.")
#         sys.exit(1)
#     comp_tree = KDTree(tree_coords, 2)

#     selected_comps_map = {}
#     diff_mags = {}
#     diff_errs = {}

#     print("\nStarting differential photometry loop...")
#     for target_idx in tqdm(unique_meta.index, desc="Processing targets"):
#         target_coords = (unique_meta.loc[target_idx, RA_COLUMN], unique_meta.loc[target_idx, DEC_COLUMN])
#         neighbors = comp_tree.get_knn(target_coords, num_comps + 10, False)

#         seen = set()
#         comps = []
#         for ra, dec, idx in neighbors:
#             if idx != target_idx and idx not in seen:
#                 comps.append(idx)
#                 seen.add(idx)
#             if len(comps) >= num_comps:
#                 break

#         selected_comps_map[target_idx] = comps

#         if not comps:
#             diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
#             diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
#             continue

#         t_mag = pivot_mags[target_idx]
#         t_err = pivot_errs[target_idx]
#         diffs = []
#         errs = []

#         for c_idx in comps:
#             c_mag = pivot_mags[c_idx]
#             c_err = pivot_errs[c_idx]
#             t_aligned, c_aligned = t_mag.align(c_mag, join='inner')
#             te_aligned, ce_aligned = t_err.align(c_err, join='inner')
#             if not t_aligned.empty:
#                 # Only compute differences where target mag is non-NaN
#                 valid_mask = t_aligned.notna()
#                 diff = pd.Series(np.nan, index=t_aligned.index)
#                 diff[valid_mask] = t_aligned[valid_mask] - c_aligned[valid_mask]
#                 err = pd.Series(np.nan, index=te_aligned.index)
#                 err[valid_mask] = np.sqrt(np.square(te_aligned[valid_mask]) + np.square(ce_aligned[valid_mask]))
#                 diffs.append(diff)
#                 errs.append(err)

#         if not diffs:
#             diff_mags[target_idx] = pd.Series(np.nan, index=pivot_mags.index)
#             diff_errs[target_idx] = pd.Series(np.nan, index=pivot_errs.index)
#             continue

#         df_diffs = pd.DataFrame(diffs).T
#         df_errs = pd.DataFrame(errs).T

#         med_diff = df_diffs.median(axis=1)
#         total_err = np.sqrt(df_errs.sum(axis=1)) / len(comps)

#         diff_mags[target_idx] = med_diff.reindex(pivot_mags.index)
#         diff_errs[target_idx] = total_err.reindex(pivot_errs.index)

#     print("\nDifferential photometry complete. Mapping and saving results...")

#     df_mags = pd.DataFrame(diff_mags)
#     df_errs = pd.DataFrame(diff_errs)

#     long_mags = df_mags.stack().rename('differential_magnitude').reset_index()
#     long_mags.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)
#     long_errs = df_errs.stack().rename('differential_magnitude_error').reset_index()
#     long_errs.rename(columns={'level_1': 'star_index', 'level_0': TIME_COLUMN}, inplace=True)

#     df = df.merge(long_mags, on=['star_index', TIME_COLUMN], how='left')
#     df = df.merge(long_errs, on=['star_index', TIME_COLUMN], how='left')

#     os.makedirs(output_dir, exist_ok=True)
#     df.to_csv(os.path.join(output_dir, output_csv), index=False)

#     with open(os.path.join(output_dir, output_json), 'w') as f:
#         json.dump(selected_comps_map, f, indent=4)

#     print(f"\nSaved results to {output_csv} and {output_json}")
