import pandas as pd
import sys
import os

# Command-line arguments
input_dir = sys.argv[1]
output_dir = sys.argv[1]
target = sys.argv[2]
date = sys.argv[3]
star_index = sys.argv[4]
input_csv = os.path.join(input_dir, f"star_data_{star_index}_{target}_{date}.csv")
output_csv = os.path.join(output_dir, f"star_data_filtered_{star_index}_{target}_{date}.csv")
output_log = os.path.join(output_dir, f"remove_points_log_{star_index}_{target}_{date}.txt")

# Column names (matching extract_star_data.py)
IMAGE_INDEX_COLUMN = 'image_index'

# Load data
try:
    data = pd.read_csv(input_csv, sep=',', engine='python')
except FileNotFoundError:
    print(f"Error: Input CSV not found at {input_csv}")
    sys.exit(1)

# Ensure image_index is integer
data[IMAGE_INDEX_COLUMN] = data[IMAGE_INDEX_COLUMN].astype(int)

# Prompt user
while True:
    print(">>> Do you want to remove any data points? (yes/no)")
    response = input("Answer: ").strip().lower()
    if response in ['yes', 'no']:
        break
    print("Error: Please enter 'yes' or 'no'")

# Initialize log
log_lines = []
log_lines.append("===== Remove Points Summary =====")
log_lines.append(f"Input CSV: {input_csv}")
log_lines.append(f"Star Index: {star_index}")
log_lines.append(f"Total rows before: {len(data)}")

if response == 'no':
    # No points to remove, save original data as filtered
    print("No points selected for removal. Saving original data as filtered CSV.")
    log_lines.append("Action: No points removed (user selected 'no')")
    os.makedirs(output_dir, exist_ok=True)
    data.to_csv(output_csv, index=False)
    log_lines.append(f"Output CSV: {output_csv}")
    log_lines.append(f"Total rows after: {len(data)}")
    log_lines.append("================================")
    with open(output_log, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Saved filtered CSV to {output_csv}")
    print(f"Saved log to {output_log}")
    sys.exit(0)

# Prompt for image_index values
while True:
    print(">>> Enter image_index values to remove (comma-separated, e.g., 1,2,3):")
    indices_input = input("Answer: ").strip()
    try:
        indices_to_remove = [int(x.strip()) for x in indices_input.split(',') if x.strip()]
        break
    except ValueError:
        print("Error: Please enter valid integers separated by commas")

# Validate image_index values
valid_indices = data[IMAGE_INDEX_COLUMN].unique()
invalid_indices = [idx for idx in indices_to_remove if idx not in valid_indices]
if invalid_indices:
    print(f"Error: Invalid image_index values: {invalid_indices}")
    log_lines.append(f"Error: Invalid image_index values provided: {invalid_indices}")
    log_lines.append("No changes made to the CSV")
    log_lines.append("================================")
    with open(output_log, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Saved log to {output_log}")
    sys.exit(1)

# Filter data
filtered_data = data[~data[IMAGE_INDEX_COLUMN].isin(indices_to_remove)].copy()

# Log details
log_lines.append(f"Image indices requested for removal: {indices_to_remove}")
log_lines.append(f"Valid indices removed: {[idx for idx in indices_to_remove if idx in valid_indices]}")
log_lines.append(f"Total rows after: {len(filtered_data)}")
if 'Simbad_Name' in data.columns and 'Simbad_Type' in data.columns:
    log_lines.append(f"Simbad Name: {data['Simbad_Name'].iloc[0]}")
    log_lines.append(f"Simbad Type: {data['Simbad_Type'].iloc[0]}")
if 'is_variable' in data.columns:
    log_lines.append(f"Is Variable: {data['is_variable'].iloc[0]}")
log_lines.append(f"Output CSV: {output_csv}")
log_lines.append("================================")

# Save filtered CSV
os.makedirs(output_dir, exist_ok=True)
filtered_data.to_csv(output_csv, index=False)
print(f"Saved filtered CSV to {output_csv}")

# Save log
with open(output_log, 'w') as f:
    f.write('\n'.join(log_lines))
print(f"Saved log to {output_log}")