import pandas as pd
import sys
import os

# --- Configuration ---
input_full_data_path = sys.argv[1] + "/tracked_stars.csv"
output_classified_full_data_path = sys.argv[2] + "/tracked_stars_with_simbad.csv"

# --- Main Script ---
print("Loading full data...")
try:
    full_star_data = pd.read_csv(input_full_data_path)
except FileNotFoundError:
    print(f"Error: Input file not found at {input_full_data_path}")
    sys.exit(1)

# Identify unique stars
unique_star_indices = full_star_data['star_index'].unique()
print(f"Found {len(unique_star_indices)} unique stars.")

# Assign default values for Simbad_Type and Simbad_Name
print("Assigning default classifications (UNCLASSIFIED, No Name)...")
unique_star_types = {idx: "UNCLASSIFIED" for idx in unique_star_indices}
unique_simbad_names = {idx: "No Name" for idx in unique_star_indices}

# Map classifications to full DataFrame
print("Mapping classifications to all observations...")
full_star_data['Simbad_Type'] = full_star_data['star_index'].map(unique_star_types)
full_star_data['Simbad_Name'] = full_star_data['star_index'].map(unique_simbad_names)

# Save final output
os.makedirs(os.path.dirname(output_classified_full_data_path), exist_ok=True)
full_star_data.to_csv(output_classified_full_data_path, index=False)
print(f"Output saved to: {output_classified_full_data_path}")