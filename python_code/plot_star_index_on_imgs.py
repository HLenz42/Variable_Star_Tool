import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.patches import Circle
import numpy as np
import re
from matplotlib import rcParams
import seaborn as sns
from matplotlib.gridspec import GridSpec

# Set matplotlib parameters
rcParams["figure.figsize"] = (10,10)
rcParams["font.size"] = 13

tracked_csv = "tracked_stars.csv"

sns.set(style="whitegrid")

def plot_quadrant_with_labeled_stars(fits_path, stars_df, output_path, quadrant_name, target, date):
    """
    Plot a quadrant FITS image with labeled circles around stars, cycling through colors.
    Labels show star_index, and color bar is removed since tracking percentage is not used.
    """
    # Load FITS file
    try:
        hdu = fits.open(fits_path)
        img_data = hdu[0].data
        img_wcs = WCS(hdu[0].header)
    except Exception as e:
        print(f"[!] Error loading FITS file {fits_path}: {e}")
        return

    # Create figure
    fig = plt.figure(figsize=(10,10))  # No color bar, so standard size
    ax = fig.add_subplot(111, projection=img_wcs)
    im = ax.imshow(img_data, origin='lower', cmap='gray', vmin=np.percentile(img_data, 5), vmax=np.percentile(img_data, 99))

    # Define a list of distinct colors to cycle through
    colors = ['red', 'cornflowerblue', 'green', 'yellow', 'cyan', 'magenta', 'lime']
    num_colors = len(colors)

    # Plot circles and labels for each star
    for idx, (_, row) in enumerate(stars_df.iterrows()):
        if pd.isna(row['ALPHA_J2000']) or pd.isna(row['DELTA_J2000']):
            continue  # Skip if coordinates are NaN
        sky_coord = SkyCoord(ra=row['ALPHA_J2000'] * u.deg, dec=row['DELTA_J2000'] * u.deg)
        pixel_x, pixel_y = img_wcs.world_to_pixel(sky_coord)

        # Cycle through colors based on index
        color = colors[idx % num_colors]

        # Draw circle
        circle = Circle((pixel_x, pixel_y), radius=8, edgecolor=color, facecolor='none', linewidth=1)
        ax.add_patch(circle)

        # Add label (star_index) slightly offset from the circle
        star_index = int(row['star_index'])
        ax.text(pixel_x + 7, pixel_y + 7, f'{star_index}', fontsize=4, color=color, ha='left', va='bottom')

    ax.set_title(f"{quadrant_name}")
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    save_path = os.path.join(output_path, f"star_index_{quadrant_name}_{target}_{date}.png")
    plt.savefig(save_path, bbox_inches='tight',dpi=400)
    plt.close()

def summarize_star_counts(df):
    """
    Print a summary of the number of stars per quadrant.
    """
    quadrant_counts = df.groupby('QUADRANT_FILE')['star_index'].nunique()
    print("\n===== Star Count Summary =====")
    print(f"Total unique stars: {df['star_index'].nunique()}")
    for quadrant, count in quadrant_counts.items():
        print(f"Quadrant {quadrant}: {count} stars")
    print("================================\n")

def extract_image_id(quadrant_file):
    """
    Extract the image identifier from QUADRANT_FILE by removing the quadrant suffix (_q#).
    Returns None if the format doesn't match.
    """
    if not isinstance(quadrant_file, str):
        return None
    # Match _q followed by 0-15, optionally followed by _wcs_se.ldac or similar
    match = re.match(r'(.+?)_q(\d{1,2})(?:_wcs_se\.ldac|_se\.ldac|_wcs)?$', quadrant_file)
    if match:
        image_id, quadrant_num = match.groups()
        if int(quadrant_num) in range(16):  # Ensure quadrant is 0-15
            return image_id
    return None

def main(csv_path, fits_folder, output_folder, target, date):
    """
    Main function to generate quadrant plots for the first image's 16 quadrants, with labeled stars.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Load CSV
    df = pd.read_csv(os.path.join(csv_path, tracked_csv))

    # Print star count summary
    summarize_star_counts(df)

    # Extract image identifiers from QUADRANT_FILE
    df['image_id'] = df['QUADRANT_FILE'].apply(extract_image_id)
    valid_image_ids = df['image_id'].dropna().unique()
    
    if not valid_image_ids.size:
        print("[!] No valid image identifiers found in QUADRANT_FILE. Check format (e.g., <image>_q0_wcs_se.ldac).")
        print(f"Available QUADRANT_FILE values: {df['QUADRANT_FILE'].dropna().unique()}")
        sys.exit(1)

    # Select the first image (alphabetically sorted for consistency)
    first_image_id = sorted(valid_image_ids)[0]
    print(f">>> Selected first image: {first_image_id}")

    # Define the 16 quadrants for the first image
    valid_quadrants = [
        f"{first_image_id}_q{i}_wcs_se.ldac" for i in range(16)
    ]

    # Filter DataFrame to include only these quadrants
    df_filtered = df[df['QUADRANT_FILE'].isin(valid_quadrants)]

    if df_filtered.empty:
        print(f"[!] No stars found for quadrants of {first_image_id}. Check QUADRANT_FILE values.")
        print(f"Expected quadrants: {valid_quadrants}")
        print(f"Available quadrants: {df['QUADRANT_FILE'].dropna().unique()}")
        sys.exit(1)

    # Process each of the 16 quadrants
    for quadrant_file in valid_quadrants:
        quadrant_name = quadrant_file.replace("_se.ldac", "").replace("_wcs_se.ldac", "")
        fits_name = quadrant_name + ".fit"
        fits_path = os.path.join(fits_folder, fits_name)

        if not os.path.exists(fits_path):
            print(f"[!] FITS file not found: {fits_path}")
            continue

        stars_in_quadrant = df_filtered[df_filtered['QUADRANT_FILE'] == quadrant_file]

        if stars_in_quadrant.empty:
            print(f"[!] No stars found in quadrant {quadrant_name}. Skipping.")
            continue

        print(f">>> Plotting {quadrant_name} with {len(stars_in_quadrant)} stars...")
        plot_quadrant_with_labeled_stars(fits_path, stars_in_quadrant, output_dir, quadrant_name, target, date)

    print(f"[✓] All valid quadrant plots generated for {first_image_id}.")

if __name__ == "__main__":
    csv_dir = sys.argv[1]          # Folder containing 'tracked_stars.csv'
    fits_dir = sys.argv[2]         # Folder with WCS-calibrated quadrant FITS files
    output_dir = sys.argv[3]       # Folder to save plots
    target = sys.argv[4]           # Target being observed
    date = sys.argv[5]             # Date observed

    main(csv_path=csv_dir, fits_folder=fits_dir, output_folder=output_dir, target=target, date=date)