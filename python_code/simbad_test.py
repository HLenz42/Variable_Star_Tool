from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import os
import sys

# Assume your input star catalog is this file
# input_star_catalog_path = "/home/hlenz/Projects/variable_star_tool/pipline_variability_detect/start_tracked_cats/extracted_star_subsets_4simbad.csv"
input_star_catalog_path = sys.argv[1] # path_2_star_tracked_cat,

# Define where you want to save the classified stars (this should be a *new* file)
# output_classified_stars_path = "/home/hlenz/Projects/variable_star_tool/pipline_variability_detect/start_tracked_cats/classified_stars_with_types.csv"
output_classified_stars_path = sys.argv[2] #path_2_star_tracked_simbad_cat

# Load your input star data
starParams = pd.read_csv(input_star_catalog_path)

radius_arcsec = float(sys.argv[3])

# You also need to define numStars, uncertainThreshold, and unclassifiedThreshold
numStars = len(starParams)
uncertainThreshold = 2.0 * u.arcsecond # Example: 2 arcseconds
unclassifiedThreshold = radius_arcsec * u.arcsecond # Example: 10 arcseconds

try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        pass

if not os.path.isfile(output_classified_stars_path):
    starTypeList = []
    # Add a list to store Simbad names
    simbadNameList = [] # NEW: List to store object names

    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('otype')
    custom_simbad.add_votable_fields('mesdistance')

    types = [] # This list stores unique OTYPEs found, useful for summary

    print("Starting Simbad classification...")
    for i in range(numStars):
        clear_output(wait=True)
        print(f'PROGRESS: {(i/(numStars-1))*100:.1f}%')

        coord1 = starParams['ALPHA_J2000'].iloc[i]
        coord2 = starParams['DELTA_J2000'].iloc[i]

        c = SkyCoord(coord1, coord2, frame='icrs', unit='deg')

        r = 2 * u.arcminute

        result_table = custom_simbad.query_region(c, radius=r)

        # --- DEBUGGING STEP (can remove after successful run, but good to keep for a bit) ---
        if i < 5:
            print(f"\n--- Debugging Star {i} ({c.to_string('hmsdms')}) ---")
            if result_table is None:
                print("Result table is None (no response from Simbad?).")
            elif len(result_table) == 0:
                print("Result table is empty (no objects found in region).")
                print(f"Querying region RA={c.ra.deg} Dec={c.dec.deg} with radius {r}")
            else:
                print(f"Columns in result_table: {result_table.colnames}")
                print(f"First row of result_table:\n{result_table[0]}")
            print("-------------------------------------------\n")
        # --- END DEBUGGING STEP ---


        try:
            if len(result_table) == 0:
                starTypeList.append('UNCLASSIFIED')
                simbadNameList.append('No Simbad Name') # NEW: No name if no object found
                continue

            # Retrieve object type (lowercase 'otype')
            starType = result_table['otype'][0]
            if starType not in types:
                types.append(starType)
            starTypeList.append(starType)

            # Retrieve object name (from 'main_id')
            simbadName = result_table['main_id'][0] # NEW: Get the main_id
            simbadNameList.append(simbadName) # NEW: Add to list

            # Access the distance using 'mesdistance.dist'
            if 'mesdistance.dist' not in result_table.colnames:
                starTypeList[-1] = 'UNCLASSIFIED'
                print(f"Warning: 'mesdistance.dist' column not found for star {i}. Cannot assess uncertainty/unclassified status for this star, marking UNCLASSIFIED.")
                continue

            # Now, proceed with distance calculations using 'mesdistance.dist'
            if len(result_table) > 1:
                distance_0 = result_table['mesdistance.dist'][0] * u.arcsecond
                distance_1 = result_table['mesdistance.dist'][1] * u.arcsecond

                if (result_table['otype'][0] == 'Planet' or result_table['otype'][1] == 'Planet'):
                    starTypeList[-1] = 'EXOPLANET'
                elif (distance_1 - distance_0) < uncertainThreshold:
                    starTypeList[-1] = 'UNCERTAIN'
                elif distance_0 > unclassifiedThreshold:
                    starTypeList[-1] = 'UNCLASSIFIED'
            elif len(result_table) == 1:
                distance_0 = result_table['mesdistance.dist'][0] * u.arcsecond
                if distance_0 > unclassifiedThreshold:
                    starTypeList[-1] = 'UNCLASSIFIED'

        except (TypeError, IndexError) as e:
            print(f"Exception during processing star {i} ({c.to_string('hmsdms')}): {e}. Classifying as UNCLASSIFIED.")
            starTypeList.append('UNCLASSIFIED')
            simbadNameList.append('Error or No Simbad Name') # NEW: Handle name in error case
            continue

    starParams['Type'] = starTypeList
    starParams['Simbad_Name'] = simbadNameList # NEW: Add the new column to DataFrame

    print("\nSimbad classification complete.")
    print("Unique types found:", types)
    print("First few classified stars:")
    print(starParams.head())

    starParams.to_csv(output_classified_stars_path, index=False)
    print(f"Classified stars saved to: {output_classified_stars_path}")

else:
    print('Already completed. Skip (classified output file already exists).')