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
                print(f"Warning: DATE-OBS missing in {fits_path}")
                return np.nan
    except Exception as e:
        print(f"Error reading FITS file {fits_path}: {e}")
        return np.nan

def estimate_median_error(catalog):
    """Estimate median positional uncertainty from ERRAWIN_IMAGE and ERRBWIN_IMAGE."""
    if 'ERRAWIN_IMAGE' in catalog.columns and 'ERRBWIN_IMAGE' in catalog.columns:
        errors = np.sqrt(catalog['ERRAWIN_IMAGE']**2 + catalog['ERRBWIN_IMAGE']**2)
        return np.nanmedian(errors) * 3600  # Convert to arcsec
    return None

def match_catalogues(csv_files, search_radius_arcsec, output_path, fits_folder, min_match_fraction):
    """Match stars across catalogues with magnitude and error weighting."""
    # Load all catalogues
    catalogues = [pd.read_csv(f) for f in csv_files]
    image_ids = [os.path.basename(f).replace('.csv', '') for f in csv_files]

    # Get Julian dates for each image
    julian_dates = {}
    for image_id in image_ids:
        fits_name = image_id.replace('_stars', '') + '.fit'
        fits_path = os.path.join(fits_folder, fits_name)
        if not os.path.exists(fits_path):
            fits_name = image_id.replace('_stars', '') + '.fits'
            fits_path = os.path.join(fits_folder, fits_name)
        julian_dates[image_id] = get_julian_date(fits_path)

    # Create SkyCoord objects for each catalogue
    skycoords = [SkyCoord(ra=cat['ALPHA_J2000'].values * u.deg, dec=cat['DELTA_J2000'].values * u.deg) for cat in catalogues]

    # Estimate median positional errors for adaptive radius
    median_errors = [estimate_median_error(cat) for cat in catalogues]
    adaptive_radius = max([r for r in median_errors if r is not None], default=search_radius_arcsec)
    print(f"[⧗] Using adaptive search radius: {adaptive_radius:.2f} arcsec (user-specified: {search_radius_arcsec:.2f})")

    # Get column names from the first catalogue
    columns = catalogues[0].columns.tolist()
    output_columns = ['image_id', 'julian_date', 'star_index'] + columns

    # Reference catalogue (first image)
    ref_cat = catalogues[0].copy()
    ref_coords = skycoords[0]
    matched_rows = []
    star_indices = []
    n_matched_images = []
    separations = []  # Track separations for diagnostics
    mag_diffs = []   # Track magnitude differences

    print(f"[⧗] Matching {len(ref_cat)} stars from {len(csv_files)} catalogues...")

    for i in tqdm(range(len(ref_cat)), desc="Matching stars"):
        ref_star = ref_coords[i]
        ref_mag = ref_cat['MAG_AUTO'].iloc[i] if 'MAG_AUTO' in ref_cat.columns else np.nan
        matched_indices = [i]  # Store index for first catalogue
        n_matched = 1

        for j in range(1, len(catalogues)):
            sep = ref_star.separation(skycoords[j])
            match_idx = np.where(sep.arcsec < adaptive_radius)[0]

            if len(match_idx) == 0:
                matched_indices.append(None)
            else:
                # If multiple matches, select based on magnitude and positional error
                if len(match_idx) > 1 and 'MAG_AUTO' in catalogues[j].columns:
                    mags = catalogues[j]['MAG_AUTO'].iloc[match_idx].values
                    mag_diff = np.abs(mags - ref_mag)
                    errors = np.sqrt(catalogues[j]['ERRAWIN_IMAGE'].iloc[match_idx]**2 +
                                     catalogues[j]['ERRBWIN_IMAGE'].iloc[match_idx]**2) * 3600
                    # Weighted score: combine positional and magnitude differences
                    weights = sep.arcsec[match_idx] / adaptive_radius + 0.5 * mag_diff / np.nanstd(mags)
                    if errors is not None:
                        weights += errors / adaptive_radius  # Penalize high uncertainty
                    closest = match_idx[np.argmin(weights)]
                else:
                    closest = match_idx[np.argmin(sep.arcsec[match_idx])]
                matched_indices.append(closest)
                n_matched += 1
                separations.append(sep.arcsec[closest])
                if 'MAG_AUTO' in catalogues[j].columns:
                    mag_diffs.append(np.abs(catalogues[j]['MAG_AUTO'].iloc[closest] - ref_mag))

        # Include stars matched in at least min_match_fraction of images
        if n_matched >= len(catalogues) * min_match_fraction:
            star_indices.append(matched_indices)
            n_matched_images.append(n_matched)

    # Create output rows
    for img_idx in range(len(image_ids)):
        for star_idx, matched_indices in enumerate(star_indices):
            row_data = {
                'image_id': image_ids[img_idx],
                'julian_date': julian_dates[image_ids[img_idx]],
                'star_index': star_idx
            }
            if matched_indices[img_idx] is not None:
                row_data.update(catalogues[img_idx].iloc[matched_indices[img_idx]].to_dict())
            matched_rows.append(row_data)

    print(f"[✓] Matched {len(star_indices)} stars across at least {min_match_fraction*100:.0f}% of the images.")

    # Save to CSV
    df_flat = pd.json_normalize(matched_rows)
    for col in output_columns:
        if col not in df_flat.columns:
            df_flat[col] = np.nan
    df_flat = df_flat[output_columns]
    df_flat.to_csv(output_path, index=False)
    print(f"Saved tracked catalogue to: {output_path}")

    # Diagnostics
    print("\n========== Star Tracking Summary ==========\n")
    print(f"Tracked {len(star_indices)} stars across at least {min_match_fraction*100:.0f}% of the images")
    print(f"Number of images: {len(catalogues)}")
    print(f"Effective search radius: {adaptive_radius:.2f} arcsec")
    if separations:
        print(f"Median separation: {np.nanmedian(separations):.2f} arcsec")
    if mag_diffs:
        print(f"Median magnitude difference: {np.nanmedian(mag_diffs):.2f} mag")
    print("==========================================\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py input_folder output_folder radius_arcsec fits_folder")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_csv = sys.argv[2]
    radius_arcsec = float(sys.argv[3])
    fits_folder = sys.argv[4]
    min_match_fraction float(sys.argv[5]) / 100

    output_file = os.path.join(output_csv, "tracked_stars.csv")
    os.makedirs(output_csv, exist_ok=True)

    all_csvs = sorted(glob.glob(os.path.join(input_folder, "*.csv")))
    if not all_csvs:
        print(f"No CSV files found in {input_folder}")
        sys.exit(1)

    match_catalogues(all_csvs, radius_arcsec, output_file, fits_folder, min_match_fraction)