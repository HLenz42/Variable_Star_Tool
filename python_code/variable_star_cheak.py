import pandas as pd
import sys 
import os

target = sys.argv[3]
date = sys.argv[4]
input_classified_data_path = sys.argv[1] + "/tracked_stars_with_simbad.csv"
output_dir = sys.argv[2]
output_data_with_variability_flag_path = output_dir + f"/tracked_stars_with_variability_flag_{target}_{date}.csv"


# --- Define the variable star types that are considered "confirmed" ---
# This should align with the 'Simbad_Type' values found in your data.
CONFIRMED_VARIABLE_SIMBAD_TYPES = [
    'EB*',       # Eclipsing Binary
    'RR*',       # RR Lyrae Variable
    'BY*',       # BY Draconis Variable
    'V*',        # General Variable Star
    'LP?',       # Low-Luminosity Pulsar Candidate (if you count candidates as 'Yes')
    'LP*',       # Low-Luminosity Pulsar
    'RS*',       # Radio Star (some can be variable due to flares/activity)
    'EXOPLANET', # Exoplanets
    'Ro*',       # Rotating variable star
    'Er*',       # Eruptive Variable
    'cC*',       # Classical Cepheid 
    'Cepheid',   # general cepheid variables stars
    'cCeph',     # classical cepheids
    'dCep',      # anomalous cepheids
    'dSct',      # delta cepheids
    'Mira',      # mira type long period variables
    'S*',        # semiregular pulsating stars
    'Sr*',       # semiregular varible stars
    'RV*',       # RV Tauri variables
    'bCep',      # beta cephei variables
    'gDor',      # gamma doradus variables
    'CV*',       # cataclysmic variables
    'No*',       # novae
    'AM*',       # AM herculies stars
    'FU*',       # FU orionis stars
    'TT*',       # T Tauri stars
    'Or*',       # Orion variables
    'UV*',       # UV ceti stars
    'a2CVn',     # alpha^2 canum venaticorum variables
    'SX*',       # SX arietis variables
    'PulsV*',    # pulsating variable stars
    'Ir*',       # Irregular variables
    'Be*',       # Be stars
]

# ADD MORE TO THE LIST ABOVE WHEN FOUND
print("\n>>> Current list of star types to be defined as variable")
print(CONFIRMED_VARIABLE_SIMBAD_TYPES)
print()

print("--- Loading data ---")
try:
    full_star_data = pd.read_csv(input_classified_data_path)
except FileNotFoundError:
    print(f"Error: Input file not found at {input_classified_data_path}")
    sys.exit(1)

# Ensure 'star_index' is treated as a numerical type if it's not already
full_star_data['star_index'] = full_star_data['star_index'].astype(int)

print("--- Categorizing Stars by Variability ---")

# --- 1. Create the 'is_variable' column ---
# Initialize the new column with 'no' for all stars.
# This operation automatically applies to all rows in the DataFrame.
full_star_data['is_variable'] = 'no'
print(">>> Initialize script by setting all stars as non-variable")

# Set 'yes' for stars whose 'Simbad_Type' is in our confirmed variable list.
# The .loc[] operation with the boolean mask will set 'Yes' for every row
# where the condition (Simbad_Type is in CONFIRMED_VARIABLE_SIMBAD_TYPES) is True.
full_star_data.loc[full_star_data['Simbad_Type'].isin(CONFIRMED_VARIABLE_SIMBAD_TYPES), 'is_variable'] = 'yes'
print(">>> Detect if Simbad_Type is in the list of variable stars")
print(">>> Set those types to variable")

# --- 2. Identify specific groups for summary print statements (these still use unique star_index) ---
# For these, we *do* want the unique star_index because we're listing stars, not individual observations.
confirmed_variable_indices = full_star_data.loc[
    full_star_data['Simbad_Type'].isin(CONFIRMED_VARIABLE_SIMBAD_TYPES),
    'star_index'
].unique().tolist()

uncertain_unclassified_indices = full_star_data.loc[
    (full_star_data['Simbad_Type'] == 'UNCERTAIN') |
    (full_star_data['Simbad_Type'] == 'UNCLASSIFIED'),
    'star_index'
].unique().tolist()

exoplanet_indices = full_star_data.loc[
    full_star_data['Simbad_Type'] == 'EXOPLANET',
    'star_index'
].unique().tolist()

# --- 3. Print summaries ---
print(f'CONFIRMED VARIABLE STARS (based on Simbad types) - {len(confirmed_variable_indices)}: \n{confirmed_variable_indices}\n')
print(f'POSSIBLE EXOPLANET SYSTEMS - {len(exoplanet_indices)}: \n{exoplanet_indices}\n')
print(f'UNCERTAIN OR UNCLASSIFIED STARS (potentially variable, need manual check) - {len(uncertain_unclassified_indices)}: \n{uncertain_unclassified_indices}\n')

print("--- 'is_variable' column added to data. Sample: ---")
print(full_star_data[['star_index', 'Simbad_Type', 'is_variable']].head(10))

# --- NEW: Save the updated DataFrame to a new CSV file ---
print(f"\nSaving updated data with 'is_variable' column to: {output_data_with_variability_flag_path}")
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)
full_star_data.to_csv(output_data_with_variability_flag_path, index=False)
print("Save complete.")
