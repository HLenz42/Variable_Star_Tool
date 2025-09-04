import pandas as pd
import os
import sys

def extract_columns_and_rows(input_csv, output_folder, col1, col2, col3, start_row, end_row):
    """
    Extract specified columns and row range from a CSV and save to a new CSV.
    
    Parameters:
    - input_csv: Path to input CSV file (e.g., 'input_folder/tracked_stars.csv')
    - output_folder: Path to output folder
    - col1, col2, col3: Names of the three columns to extract
    - start_row, end_row: Row indices (0-based) to extract (inclusive start, exclusive end)
    """
    # Ensure output directory exists
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, "extracted_star_subsets_4simbad.csv")

    # Check if output file already exists
    if os.path.isfile(output_path):
        print(f"[!] Output file {output_path} already exists. Please remove or specify a different output folder.")
        sys.exit(1)

    # Load the CSV
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"[!] Input CSV {input_csv} not found.")
        sys.exit(1)

    # Validate columns
    required_columns = [col1, col2, col3]
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        print(f"[!] Missing columns: {missing_cols}")
        sys.exit(1)

    # Validate row range
    num_rows = len(df)
    if start_row < 0 or end_row > num_rows or start_row >= end_row:
        print(f"[!] Invalid row range: start_row={start_row}, end_row={end_row}, total_rows={num_rows}")
        sys.exit(1)

    # Extract specified columns and rows
    extracted_df = df.loc[start_row:end_row-1, required_columns]

    # Save to output CSV
    extracted_df.to_csv(output_path, index=False)
    print(f"[✓] Extracted {len(extracted_df)} rows with columns {required_columns} to {output_path}")

if __name__ == "__main__":
    # if len(sys.argv) != 7:
    #     print("Usage: python extract_columns.py <input_folder> <output_folder> <col1> <col2> <col3> <start_row> <end_row>")
    #     sys.exit(1)

    input_folder = sys.argv[1]
    print(f"input_folder: {input_folder}")
    output_folder = sys.argv[2]
    print(f"output_folder: {output_folder}")
    col1 = sys.argv[3]
    print(f"col1: {col1}")
    col2 = sys.argv[4]
    print(f"col2: {col2}")
    col3 = sys.argv[5]
    print(f"col3: {col3}")
    try:
        start_row = float(sys.argv[6])
        print(f"start_row: {start_row}")
        end_row = float(sys.argv[7])
        print(f"end_row: {end_row}")
    except ValueError:
        print("[!] start_row and end_row must be integers")
        sys.exit(1)

    input_csv = os.path.join(input_folder, "tracked_stars.csv")
    extract_columns_and_rows(input_csv, output_folder, col1, col2, col3, start_row, end_row)