import pandas as pd
import os
import sys

# === Configuration ===
input_path = sys.argv[1] + f"/tracked_stars_with_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
output_dir = sys.argv[2]
output_csv = f"tracked_stars_with_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
# output_csv = f"tracked_stars_with_manual_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
output_log = f"manual_variability_update_log_{sys.argv[3]}_{sys.argv[4]}.txt"

STAR_INDEX_COLUMN = 'star_index'

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

    # Verify required columns
    if 'is_variable' not in df.columns:
        print("Error: 'is_variable' column not found in data.")
        sys.exit(1)

    # Get unique star indices
    unique_star_indices = df[STAR_INDEX_COLUMN].unique()
    print(f"Found {len(unique_star_indices)} unique stars.")

    # Initialize lists to track changes
    variable_stars_updated = []
    non_variable_stars_updated = []
    log_lines = []
    log_lines.append("===== Manual Variability Flag Update Summary =====")
    log_lines.append(f"Input CSV: {input_path}")
    log_lines.append(f"Total unique stars: {len(unique_star_indices)}")

    # Function to validate user input for star index
    def validate_star_index(user_input, valid_indices):
        try:
            star_idx = int(user_input)
            if star_idx in valid_indices:
                return star_idx
            else:
                print(f"Error: Star index {star_idx} not found in the dataset.")
                return None
        except ValueError:
            print("Error: Please enter a valid integer for star index.")
            return None

    # Prompt for variable stars
    print("\n=== Assign Variable Stars ===")
    while True:
        response = input("Would you like to assign a star index as variable? (yes/no): ").strip().lower()
        if response != 'yes':
            break
        star_input = input("Enter star index: ").strip()
        star_idx = validate_star_index(star_input, unique_star_indices)
        if star_idx is not None:
            # Check current is_variable value
            current_value = df[df[STAR_INDEX_COLUMN] == star_idx]['is_variable'].iloc[0].lower()
            if current_value != 'yes':
                df.loc[df[STAR_INDEX_COLUMN] == star_idx, 'is_variable'] = 'yes'
                variable_stars_updated.append(star_idx)
                print(f"Star {star_idx} set to variable (is_variable = 'yes').")
            else:
                print(f"Star {star_idx} is already marked as variable.")

    # Prompt for non-variable stars
    print("\n=== Assign Non-Variable Stars ===")
    while True:
        response = input("Would you like to assign a star index as non-variable? (yes/no): ").strip().lower()
        if response != 'yes':
            break
        star_input = input("Enter star index: ").strip()
        star_idx = validate_star_index(star_input, unique_star_indices)
        if star_idx is not None:
            # Check current is_variable value
            current_value = df[df[STAR_INDEX_COLUMN] == star_idx]['is_variable'].iloc[0].lower()
            if current_value != 'no':
                df.loc[df[STAR_INDEX_COLUMN] == star_idx, 'is_variable'] = 'no'
                non_variable_stars_updated.append(star_idx)
                print(f"Star {star_idx} set to non-variable (is_variable = 'no').")
            else:
                print(f"Star {star_idx} is already marked as non-variable.")

    # Log changes
    log_lines.append("\nVariable stars updated (set to 'yes'):")
    log_lines.append(f"Total variable stars updated: {len(variable_stars_updated)}")
    if variable_stars_updated:
        log_lines.append(f"Variable star indices: {', '.join(map(str, sorted(variable_stars_updated)))}")
        variable_types = df[df[STAR_INDEX_COLUMN].isin(variable_stars_updated)]['Simbad_Type'].value_counts()
        log_lines.append("Variable stars by Simbad_Type:")
        for v_type, count in variable_types.items():
            log_lines.append(f"  {v_type}: {count} stars")
    else:
        log_lines.append("No stars set to variable.")

    log_lines.append("\nNon-variable stars updated (set to 'no'):")
    log_lines.append(f"Total non-variable stars updated: {len(non_variable_stars_updated)}")
    if non_variable_stars_updated:
        log_lines.append(f"Non-variable star indices: {', '.join(map(str, sorted(non_variable_stars_updated)))}")
        non_variable_types = df[df[STAR_INDEX_COLUMN].isin(non_variable_stars_updated)]['Simbad_Type'].value_counts()
        log_lines.append("Non-variable stars by Simbad_Type:")
        for v_type, count in non_variable_types.items():
            log_lines.append(f"  {v_type}: {count} stars")
    else:
        log_lines.append("No stars set to non-variable.")

    # Summary of final variable stars
    final_variable_stars = df[df['is_variable'] == 'yes'][STAR_INDEX_COLUMN].unique()
    log_lines.append("\nFinal variable stars (after updates):")
    log_lines.append(f"Total variable stars: {len(final_variable_stars)}")
    log_lines.append(f"Final variable star indices: {', '.join(map(str, sorted(final_variable_stars)))}")
    final_variable_types = df[df['is_variable'] == 'yes']['Simbad_Type'].value_counts()
    log_lines.append("Final variable stars by Simbad_Type:")
    for v_type, count in final_variable_types.items():
        log_lines.append(f"  {v_type}: {count} stars")

    log_lines.append(f"\nOutput CSV: {os.path.join(output_dir, output_csv)}")
    log_lines.append("=========================================")

    # Remove differential_magnitude and differential_magnitude_error columns
    df = df.drop(columns=['differential_magnitude', 'differential_magnitude_error'], errors='ignore')

    # Save updated CSV
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, output_csv), index=False)
    print(f"Updated CSV saved to: {os.path.join(output_dir, output_csv)}")

    # Save log file
    with open(os.path.join(output_dir, output_log), 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Log saved to: {os.path.join(output_dir, output_log)}")