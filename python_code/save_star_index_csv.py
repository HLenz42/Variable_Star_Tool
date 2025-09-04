import pandas as pd
import os
import sys

# === Configuration ===
input_path = sys.argv[1] + f"/tracked_stars_with_differential_mags_{sys.argv[3]}_{sys.argv[4]}.csv"
output_dir = sys.argv[2]
star_index = int(sys.argv[5])
output_csv = f"star_data_{star_index}_{sys.argv[3]}_{sys.argv[4]}.csv"
output_log = f"star_data_extraction_log_{star_index}_{sys.argv[3]}_{sys.argv[4]}.txt"

# Column names (matching your pipeline)
STAR_INDEX_COLUMN = 'star_index'
TIME_COLUMN = 'julian_date'

# === Main Script ===
if __name__ == "__main__":
    # Load data
    print("Loading data...")
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        sys.exit(1)

    df[STAR_INDEX_COLUMN] = df[STAR_INDEX_COLUMN].astype(int)

    # Verify star_index exists
    unique_star_indices = df[STAR_INDEX_COLUMN].unique()
    if star_index not in unique_star_indices:
        print(f"Error: Star index {star_index} not found in the dataset.")
        sys.exit(1)

    # Filter data for the specified star_index and sort by julian_date
    print(f"Extracting data for star {star_index}...")
    star_data = df[df[STAR_INDEX_COLUMN] == star_index].copy()
    star_data = star_data.sort_values(by=TIME_COLUMN).reset_index(drop=True)

    # Add image_index column (1, 2, 3, ...)
    star_data['image_index'] = range(1, len(star_data) + 1)

    # Log file content
    log_lines = []
    log_lines.append("===== Star Data Extraction Summary =====")
    log_lines.append(f"Input CSV: {input_path}")
    log_lines.append(f"Requested star index: {star_index}")
    log_lines.append(f"Number of observations (images) for star {star_index}: {len(star_data)}")
    log_lines.append(f"Columns in output: {', '.join(star_data.columns)}")
    if 'Simbad_Type' in star_data.columns and 'Simbad_Name' in star_data.columns:
        simbad_type = star_data['Simbad_Type'].iloc[0]
        simbad_name = star_data['Simbad_Name'].iloc[0]
        is_variable = star_data['is_variable'].iloc[0] if 'is_variable' in star_data.columns else 'N/A'
        log_lines.append(f"Simbad Type: {simbad_type}")
        log_lines.append(f"Simbad Name: {simbad_name}")
        log_lines.append(f"Is Variable: {is_variable}")
    log_lines.append(f"Output CSV: {os.path.join(output_dir, output_csv)}")
    log_lines.append("=======================================")

    # Save filtered CSV
    os.makedirs(output_dir, exist_ok=True)
    star_data.to_csv(os.path.join(output_dir, output_csv), index=False)
    print(f"Star data saved to: {os.path.join(output_dir, output_csv)}")

    # Save log file
    with open(os.path.join(output_dir, output_log), 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Log saved to: {os.path.join(output_dir, output_log)}")