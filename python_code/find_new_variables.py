import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import os
import sys

# === Configuration ===
input_path = sys.argv[1] + f"/tracked_stars_with_differential_mags_{sys.argv[3]}_{sys.argv[4]}.csv"
output_dir = sys.argv[2]
output_csv = f"tracked_stars_with_updated_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
output_log = f"variability_flag_update_log_{sys.argv[3]}_{sys.argv[4]}.txt"

# Column names (matching diff_phot.py and plot_both_lightcurves.py)
# DIFF_MAG_COLUMN = 'differential_magnitude'
# DIFF_MAG_ERROR_COLUMN = 'differential_magnitude_error'
TIME_COLUMN = 'julian_date'
STAR_INDEX_COLUMN = 'star_index'

while True:
    print(">>> Use weighted or unweighted differential magnitudes to find new variables?")
    print("(1) unweighted")
    print("(2) weighted\n")

    which_mag = int(input("Answer: "))
    print()

    if which_mag == 1:
        DIFF_MAG_COLUMN = 'differential_mag'
        DIFF_MAG_ERROR_COLUMN = 'differential_mag_err'
        break
    if which_mag == 2:
        DIFF_MAG_COLUMN = 'weighted_differential_mag'
        DIFF_MAG_ERROR_COLUMN = 'weighted_differential_mag_err'
        break
    else:
        print("Error: Please select either 1, 2" + "\n")

# Display selected magnitude columns
print(">>> Magnitude option selection")
print(f"    DIFF_MAG_COLUMN: {DIFF_MAG_COLUMN}")
print(f"    DIFF_MAG_ERROR_COLUMN: {DIFF_MAG_ERROR_COLUMN}\n")


# Variability detection parameters (matching plot_both_lightcurves.py)
VARIABILITY_THRESHOLD_FACTOR = 1.5
SG_WINDOW_LENGTH_PLOT = 15
SG_POLYNOMIAL_ORDER_PLOT = 2

# === Main Script ===
if __name__ == "__main__":
    print("Loading data...")
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        sys.exit(1)

    df[STAR_INDEX_COLUMN] = df[STAR_INDEX_COLUMN].astype(int)
    df = df.sort_values(by=[TIME_COLUMN, STAR_INDEX_COLUMN]).reset_index(drop=True)

    # Verify required columns
    if 'is_variable' not in df.columns:
        print("Error: 'is_variable' column not found in data.")
        sys.exit(1)

    # Get unique star indices
    unique_star_indices = df[STAR_INDEX_COLUMN].unique()
    num_unique_stars = len(unique_star_indices)

    # Pivot differential magnitudes
    print("Pivoting differential magnitudes...")
    star_magnitudes_pivot = df.pivot_table(
        index=TIME_COLUMN, columns=STAR_INDEX_COLUMN, values=DIFF_MAG_COLUMN
    )

    # Calculate variability (max - min of smoothed differential magnitudes)
    print("\nCalculating variability for each star...")
    variation_dict = {}
    for i, star_idx in enumerate(unique_star_indices):
        print(f'PROGRESS: {((i+1)/num_unique_stars)*100:.1f}% ({i+1}/{num_unique_stars} stars processed)')
        diff_mags = star_magnitudes_pivot[star_idx].dropna()
        if len(diff_mags) >= SG_WINDOW_LENGTH_PLOT:
            smoothed_mags = savgol_filter(diff_mags.to_numpy(), SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
            variation_dict[star_idx] = smoothed_mags.max() - smoothed_mags.min()
        else:
            variation_dict[star_idx] = np.nan

    # Filter out NaN variations
    variation_dict = {k: v for k, v in variation_dict.items() if not np.isnan(v)}
    if not variation_dict:
        print("ERROR: No valid variability measurements. Cannot proceed.")
        sys.exit(1)

    # Classify new possible variables
    print("\nClassifying new possible variables...")
    variability_std = np.std(list(variation_dict.values()))
    if variability_std == 0:
        print("Warning: All stars have zero variability. No new variables flagged.")
        new_variable_stars = []
    else:
        threshold_value = VARIABILITY_THRESHOLD_FACTOR * variability_std
        print(f"Variability threshold: {threshold_value:.4f} ({VARIABILITY_THRESHOLD_FACTOR} * std dev)")
        new_variable_stars = []
        for star_idx, variation in variation_dict.items():
            is_variable = df[df[STAR_INDEX_COLUMN] == star_idx]['is_variable'].iloc[0].lower()
            if variation > threshold_value and is_variable != 'yes':
                new_variable_stars.append(star_idx)

    print(f"Found {len(new_variable_stars)} new possible variable stars.")

    # Update is_variable column directly
    print("\nUpdating variability flags...")
    df.loc[df[STAR_INDEX_COLUMN].isin(new_variable_stars), 'is_variable'] = 'yes'

    # Log file content
    log_lines = []
    log_lines.append("===== Variability Flag Update Summary =====")
    log_lines.append(f"Input CSV: {input_path}")
    log_lines.append(f"Total unique stars: {num_unique_stars}")
    log_lines.append(f"Total images (julian_date): {df[TIME_COLUMN].nunique()}")
    log_lines.append(f"Variability threshold: {threshold_value:.4f} ({VARIABILITY_THRESHOLD_FACTOR} * std dev)")
    log_lines.append(f"New possible variable stars: {len(new_variable_stars)}")
    log_lines.append(f"New variable star indices: {', '.join(map(str, sorted(new_variable_stars)))}")
    log_lines.append("\nNew variable stars by Simbad_Type:")
    new_variable_types = df[df[STAR_INDEX_COLUMN].isin(new_variable_stars)]['Simbad_Type'].value_counts()
    for v_type, count in new_variable_types.items():
        log_lines.append(f"  {v_type}: {count} stars")
    log_lines.append("\nOriginal variable stars (Simbad, before update):")
    try:
        original_df = pd.read_csv(
            sys.argv[1] + f"/tracked_stars_with_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
        )
        original_variable_stars = original_df[original_df['is_variable'] == 'yes'][STAR_INDEX_COLUMN].unique()
        log_lines.append(f"Total original variable stars: {len(original_variable_stars)}")
        log_lines.append(f"Original variable star indices: {', '.join(map(str, sorted(original_variable_stars)))}")
        original_variable_types = original_df[original_df['is_variable'] == 'yes']['Simbad_Type'].value_counts()
        for v_type, count in original_variable_types.items():
            log_lines.append(f"  {v_type}: {count} stars")
    except FileNotFoundError:
        log_lines.append("Error: Could not load original variability flag CSV for comparison.")
        original_variable_stars = []
    log_lines.append("\nUpdated variable stars (Simbad + photometry):")
    updated_variable_stars = df[df['is_variable'] == 'yes'][STAR_INDEX_COLUMN].unique()
    log_lines.append(f"Total updated variable stars: {len(updated_variable_stars)}")
    log_lines.append(f"Updated variable star indices: {', '.join(map(str, sorted(updated_variable_stars)))}")
    updated_variable_types = df[df['is_variable'] == 'yes']['Simbad_Type'].value_counts()
    for v_type, count in updated_variable_types.items():
        log_lines.append(f"  {v_type}: {count} stars")
    log_lines.append(f"\nOutput CSV: {os.path.join(output_dir, output_csv)}")
    log_lines.append("=========================================\n")

    # Remove differential_magnitude and differential_magnitude_error columns
    df = df.drop(columns=['differential_mag_err','differential_mag','weighted_differential_mag_err','weighted_differential_mag'], errors='ignore')

    # Save updated CSV
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, output_csv), index=False)
    print(f"Updated CSV saved to: {os.path.join(output_dir, output_csv)}")

    # Save log file
    with open(os.path.join(output_dir, output_log), 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Log saved to: {os.path.join(output_dir, output_log)}")





