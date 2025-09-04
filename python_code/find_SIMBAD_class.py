from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import sys
import os
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# --- Configuration ---
input_full_data_path = sys.argv[1] + "/tracked_stars.csv"
output_classified_full_data_path = sys.argv[2] + "/tracked_stars_with_simbad.csv"
temp_output_path = sys.argv[2] + "/tracked_stars_with_simbad_temp.csv"  # Temporary file for partial results
radius_arcsec = float(sys.argv[3])
search_radius = radius_arcsec * u.arcminute
uncertainThreshold = 2.0 * u.arcsecond
unclassifiedThreshold = 10.0 * u.arcsecond

# Retry configuration
MAX_RETRIES = 3
RETRY_BACKOFF_FACTOR = 2  # Delay doubles each retry (2s, 4s, 8s)
QUERY_DELAY = 1  # Seconds between queries to avoid rate limiting

# --- End Configuration ---

try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        pass

# Configure SIMBAD with retries
custom_simbad = Simbad()
custom_simbad.add_votable_fields('otype')
custom_simbad.add_votable_fields('mesdistance')

# Add retry mechanism to SIMBAD session
session = custom_simbad._session
retries = Retry(total=MAX_RETRIES, backoff_factor=RETRY_BACKOFF_FACTOR, status_forcelist=[429, 500, 502, 503, 504])
adapter = HTTPAdapter(max_retries=retries)
session.mount('http://', adapter)
session.mount('https://', adapter)

# Load existing classifications if available
classified_stars = {}
if os.path.exists(temp_output_path):
    print(f"Loading existing classifications from {temp_output_path}...")
    temp_df = pd.read_csv(temp_output_path)
    classified_stars = dict(zip(temp_df['star_index'], zip(temp_df['Simbad_Type'], temp_df['Simbad_Name'])))

print("Loading full data...")
full_star_data = pd.read_csv(input_full_data_path)

# Identify unique stars and their first occurrence coordinates
unique_star_indices = full_star_data['star_index'].unique()
print(f"Found {len(unique_star_indices)} unique stars for classification.")

first_observations = full_star_data.drop_duplicates(subset=['star_index'])
first_observations.set_index('star_index', inplace=True)

# Initialize dictionaries to store results
unique_star_types = {k: v[0] for k, v in classified_stars.items()}
unique_simbad_names = {k: v[1] for k, v in classified_stars.items()}
all_found_otypes = []  # Track unique OTYPEs

# Skip already classified stars
remaining_stars = [idx for idx in unique_star_indices if idx not in classified_stars]
print(f"Found {len(classified_stars)} already classified stars. Processing {len(remaining_stars)} remaining stars.")

print("Starting Simbad classification for unique stars...")
for i, star_idx in enumerate(remaining_stars):
    clear_output(wait=True)
    print(f'PROGRESS (Unique Stars): {((i+1+len(classified_stars))/len(unique_star_indices))*100:.1f}% ({i+1}/{len(remaining_stars)} remaining)')

    coord1 = first_observations.loc[star_idx, 'ALPHA_J2000']
    coord2 = first_observations.loc[star_idx, 'DELTA_J2000']
    c = SkyCoord(coord1, coord2, frame='icrs', unit='deg')

    for attempt in range(MAX_RETRIES):
        try:
            result_table = custom_simbad.query_region(c, radius=search_radius)
            break  # Success, exit retry loop
        except Exception as e:
            print(f"Error querying star {star_idx} (attempt {attempt+1}/{MAX_RETRIES}): {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF_FACTOR * (2 ** attempt))
            else:
                print(f"Failed to query star {star_idx} after {MAX_RETRIES} attempts. Classifying as UNCLASSIFIED.")
                result_table = None

    # Process query results
    current_star_type = 'UNCLASSIFIED'
    current_simbad_name = 'No Simbad Name'

    try:
        if result_table is not None and len(result_table) > 0:
            if 'main_id' in result_table.colnames:
                current_simbad_name = result_table['main_id'][0]

            if 'otype' in result_table.colnames:
                simbad_otype = result_table['otype'][0]
                current_star_type = simbad_otype
                if simbad_otype not in all_found_otypes:
                    all_found_otypes.append(simbad_otype)

                if len(result_table) > 1 and (result_table['otype'][0] == 'Planet' or result_table['otype'][1] == 'Planet'):
                    current_star_type = 'EXOPLANET'

                if 'mesdistance.dist' in result_table.colnames:
                    distance_0 = result_table['mesdistance.dist'][0] * u.arcsecond

                    if len(result_table) > 1:
                        distance_1 = result_table['mesdistance.dist'][1] * u.arcsecond
                        if (distance_1 - distance_0) < uncertainThreshold:
                            current_star_type = 'UNCERTAIN'
                    
                    if distance_0 > unclassifiedThreshold:
                        current_star_type = 'UNCLASSIFIED'
            else:
                print(f"Warning: 'otype' column not found for star {star_idx}. Classified as UNCLASSIFIED.")
        else:
            print(f"No results for star {star_idx}. Classified as UNCLASSIFIED.")

    except (TypeError, IndexError) as e:
        print(f"Exception during processing star {star_idx} ({c.to_string('hmsdms')}): {e}. Classifying as UNCLASSIFIED.")
        current_star_type = 'UNCLASSIFIED'
        current_simbad_name = 'Error during Simbad query'

    unique_star_types[star_idx] = current_star_type
    unique_simbad_names[star_idx] = current_simbad_name

    # Save partial results every 100 stars
    if (i + 1) % 100 == 0:
        temp_df = pd.DataFrame({
            'star_index': list(unique_star_types.keys()),
            'Simbad_Type': list(unique_star_types.values()),
            'Simbad_Name': list(unique_simbad_names.values())
        })
        temp_df.to_csv(temp_output_path, index=False)
        print(f"Saved partial results to {temp_output_path} after {i + 1 + len(classified_stars)} stars.")

    # Delay to avoid rate limiting
    time.sleep(QUERY_DELAY)

# Save final partial results
temp_df = pd.DataFrame({
    'star_index': list(unique_star_types.keys()),
    'Simbad_Type': list(unique_star_types.values()),
    'Simbad_Name': list(unique_simbad_names.values())
})
temp_df.to_csv(temp_output_path, index=False)
print(f"Saved final partial results to {temp_output_path}")

# Map classifications to full DataFrame
print("Mapping classifications to all observations...")
full_star_data['Simbad_Type'] = full_star_data['star_index'].map(unique_star_types)
full_star_data['Simbad_Name'] = full_star_data['star_index'].map(unique_simbad_names)

# Save final output
full_star_data.to_csv(output_classified_full_data_path, index=False)
print(f"Classified full data saved to: {output_classified_full_data_path}")
print("Unique Simbad OTYPEs found:", all_found_otypes)


