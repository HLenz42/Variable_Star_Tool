import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
import glob
import os
from tqdm import tqdm
import sys


def get_julian_date(fits_path):
    """Extract DATE-OBS from FITS header and convert to Julian date."""
    try:
        with fits.open(fits_path) as hdul:
            header = hdul[0].header
            date_obs = header.get('DATE-OBS')
            if date_obs:
                t = Time(date_obs, format='isot', scale='utc')
                return t.jd
            else:
                return np.nan
    except Exception as e:
        print(f"Error reading FITS file {fits_path}: {e}")
        return np.nan


def match_catalogues(csv_files, search_radius_arcsec, output_path, fits_folder):
    # Load all catalogues
    catalogues = [pd.read_csv(f) for f in csv_files]
    image_ids = [os.path.basename(f).replace('.csv', '') for f in csv_files]

    # Get Julian dates for each image
    julian_dates = {}
    for image_id in image_ids:
        # Construct FITS filename: remove '_stars' and add '.fit'
        fits_name = image_id.replace('_stars', '') + '.fit'
        
        fits_path = os.path.join(fits_folder, fits_name)
        julian_dates[image_id] = get_julian_date(fits_path)
        fits_name = image_id.replace('_stars', '') + '.fits'
    # Assume each CSV has 'RA', 'Dec' columns in degrees
    skycoords = [SkyCoord(ra=cat['ALPHA_J2000'].values * u.deg, dec=cat['DELTA_J2000'].values * u.deg) for cat in catalogues]

    # Get the column names from the first catalogue (assuming all catalogues have the same columns)
    columns = catalogues[0].columns.tolist()
    output_columns = ['image_id', 'julian_date', 'star_index'] + columns  # Add julian_date to output columns

    # Start with first catalogue as the reference
    ref_cat = catalogues[0].copy()
    ref_coords = skycoords[0]
    matched_rows = []
    star_indices = []  # Track star indices for matched stars
    n_matched_images = []  # Track number of images each star was matched in

    print(f"[⧗] Matching stars from {len(csv_files)} catalogues within {search_radius_arcsec} arcsec...")

    for i in tqdm(range(len(ref_cat))):
        ref_star = ref_coords[i]
        matched_indices = [i]  # Store index for first catalogue
        n_matched = 1  # Count first image match

        # Try to match this star in each subsequent catalogue
        for j in range(1, len(catalogues)):
            sep = ref_star.separation(skycoords[j])
            match_idx = np.where(sep.arcsec < search_radius_arcsec)[0]

            if len(match_idx) == 0:
                matched_indices.append(None)  # No match found
            else:
                closest = match_idx[np.argmin(sep.arcsec[match_idx])]
                matched_indices.append(closest)
                n_matched += 1

        # Include stars matched in at least 80% of images (adjust threshold as needed)
        if n_matched >= len(catalogues) * (min_track_percent/100):
            star_indices.append(matched_indices)
            n_matched_images.append(n_matched)

    # Create rows by image
    for img_idx in range(len(image_ids)):
        for star_idx, matched_indices in enumerate(star_indices):
            row_data = {
                'image_id': image_ids[img_idx],    # LEADS TO ORIGINAL IMAGE
                'julian_date': julian_dates[image_ids[img_idx]],  # Add Julian date
                'star_index': star_idx  # Sequential index for matched stars
            }
            if matched_indices[img_idx] is not None:
                # Add data if index exists
                for key, value in catalogues[img_idx].iloc[matched_indices[img_idx]].to_dict().items():
                    row_data[key] = value
            # Else, row_data remains with only image_id, julian_date, and star_index, other fields are NaN (implicitly)
            matched_rows.append(row_data)

    print(f"[✓] Matched {len(star_indices)} stars across at least 80% of the images.")

    # Save to CSV — use pandas.json_normalize to flatten
    df_flat = pd.json_normalize(matched_rows)
    # Ensure all columns are present, fill missing with NaN
    for col in output_columns:
        if col not in df_flat.columns:
            df_flat[col] = np.nan
    df_flat = df_flat[output_columns]  # Reorder columns
    df_flat.to_csv(output_path, index=False)
    print(f"Saved tracked catalogue to: {output_path}")

    # === Summary ===
    print("\n========== Star Tracking Summary ==========\n")
    print(f"Tracked {len(star_indices)} stars across at least {min_track_percent}% of the images")
    print(f"Number of images used: {len(catalogues)}")
    print(f"Matching radius used: {search_radius_arcsec} arcseconds")
    print(f"Data structure: Data grouped by image, with {len(star_indices)} rows per image block, "
          f"blank rows (NaN for data fields) for unmatched stars")
    print("==========================================\n")


if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_csv = sys.argv[2]
    radius_arcsec = float(sys.argv[3])
    fits_folder = sys.argv[4] 
    min_track_percent = float(sys.argv[5])

    # Fix: Ensure output_csv is a file, not a directory
    output_file = os.path.join(output_csv, "tracked_stars.csv")  # Append filename
    os.makedirs(output_csv, exist_ok=True)  # Ensure output directory exists

    all_csvs = sorted(glob.glob(os.path.join(input_folder, "*.csv")))

    match_catalogues(all_csvs, search_radius_arcsec=radius_arcsec, output_path=output_file, fits_folder=fits_folder)