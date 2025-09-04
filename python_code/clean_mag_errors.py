import pandas as pd
import numpy as np
import os
import sys

# === Configuration ===
target = sys.argv[2]  # Target star identifier
date = sys.argv[3]    # Observation date
input_path = sys.argv[1] + f"/tracked_stars_with_variability_flag_{target}_{date}.csv"
output_dir = sys.argv[1]
output_csv = f"tracked_stars_with_variability_flag_{target}_{date}.csv"  # Overwrite input CSV
output_log = f"mag_error_cleaning_log_{target}_{date}.txt"

# Column names (matching diff_phot_overwrite_csv.py input)
TIME_COLUMN = 'julian_date'
STAR_INDEX_COLUMN = 'star_index'
AUTO_MAG_COLUMN = 'MAG_AUTO'
AUTO_ERROR_COLUMN = 'MAGERR_AUTO'
APER_MAG_COLUMN = 'MAG_APER'
APER_ERROR_COLUMN = 'MAGERR_APER'

# === Main Script ===
if __name__ == "__main__":
    print("Loading data...")
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        sys.exit(1)

    # Ensure star_index is integer and sort by time and star index
    df[STAR_INDEX_COLUMN] = df[STAR_INDEX_COLUMN].astype(int)
    df = df.sort_values(by=[TIME_COLUMN, STAR_INDEX_COLUMN]).reset_index(drop=True)

    # Verify required columns
    required_columns = [STAR_INDEX_COLUMN, TIME_COLUMN, AUTO_MAG_COLUMN, AUTO_ERROR_COLUMN, 
                        APER_MAG_COLUMN, APER_ERROR_COLUMN]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing required columns: {', '.join(missing_columns)}")
        sys.exit(1)

    # Calculate and print statistics for error columns
    print("\nError statistics:")
    auto_err = df[AUTO_ERROR_COLUMN].dropna()
    aper_err = df[APER_ERROR_COLUMN].dropna()
    
    print(f"\nAutomatic aperture magnitude error ({AUTO_ERROR_COLUMN}):")
    print(f"  Maximum: {auto_err.max():.4f}")
    print(f"  Minimum: {auto_err.min():.4f}")
    print(f"  Average: {auto_err.mean():.4f}")
    
    print(f"\nFixed aperture magnitude error ({APER_ERROR_COLUMN}):")
    print(f"  Maximum: {aper_err.max():.4f}")
    print(f"  Minimum: {aper_err.min():.4f}")
    print(f"  Average: {aper_err.mean():.4f}")

    # Prompt user to choose error column for threshold
    while True:
        print("\n>>> Define error limit based on which column?")
        print("(1) Automatic aperture magnitude error (MAGERR_AUTO)")
        print("(2) Fixed aperture magnitude error (MAGERR_APER)\n")
        which_err = int(input("Answer: "))
        print()

        if which_err == 1:
            ERROR_COLUMN = AUTO_ERROR_COLUMN
            break
        if which_err == 2:
            ERROR_COLUMN = APER_ERROR_COLUMN
            break
        else:
            print("Error: Please select either 1 or 2\n")

    # Prompt user for error threshold
    while True:
        try:
            error_threshold = float(input(f"Enter error threshold for {ERROR_COLUMN}: "))
            if error_threshold >= 0:
                break
            else:
                print("Please enter a non-negative number.")
        except ValueError:
            print("Please enter a valid number.")

    print(f"\nSelected error column: {ERROR_COLUMN}")
    print(f"Error threshold: {error_threshold:.4f}")

    # Initialize log
    log_lines = []
    log_lines.append("===== Magnitude Error Cleaning Summary =====")
    log_lines.append(f"Input CSV: {input_path}")
    log_lines.append(f"Total unique stars: {df[STAR_INDEX_COLUMN].nunique()}")
    log_lines.append(f"Total images (julian_date): {df[TIME_COLUMN].nunique()}")
    log_lines.append(f"Selected error column: {ERROR_COLUMN}")
    log_lines.append(f"Error threshold: {error_threshold:.4f}")
    log_lines.append("\nPoints with errors >= threshold set to NaN:")

    # Identify and clean points with errors >= threshold
    print("\nCleaning high-error points...")
    total_cleaned = 0
    for star_idx in df[STAR_INDEX_COLUMN].unique():
        star_data = df[df[STAR_INDEX_COLUMN] == star_idx][[TIME_COLUMN, ERROR_COLUMN]].dropna(subset=[ERROR_COLUMN])
        if len(star_data) == 0:
            log_lines.append(f"Star {star_idx}: Skipped (no valid error data)")
            continue

        # Find points where error >= threshold
        high_error_mask = star_data[ERROR_COLUMN] >= error_threshold
        num_high_errors = high_error_mask.sum()

        if num_high_errors > 0:
            # Get indices in original DataFrame for these points
            high_error_indices = star_data[high_error_mask].index
            # Set all four magnitude columns to NaN
            df.loc[high_error_indices, [AUTO_MAG_COLUMN, AUTO_ERROR_COLUMN, 
                                        APER_MAG_COLUMN, APER_ERROR_COLUMN]] = np.nan
            total_cleaned += num_high_errors
            log_lines.append(f"Star {star_idx}: {num_high_errors} points set to NaN at JD: {', '.join(map(str, star_data[high_error_mask][TIME_COLUMN].values))}")
        else:
            log_lines.append(f"Star {star_idx}: No high-error points detected")

    log_lines.append(f"\nTotal points cleaned: {total_cleaned}")
    log_lines.append(f"Output CSV: {os.path.join(output_dir, output_csv)}")
    log_lines.append("============================================\n")

    # Save updated CSV, overwriting the input file
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_csv)
    df.to_csv(output_path, index=False, mode='w')
    print(f"Updated CSV saved to: {output_path}")

    # Save log file
    with open(os.path.join(output_dir, output_log), 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Log saved to: {os.path.join(output_dir, output_log)}")