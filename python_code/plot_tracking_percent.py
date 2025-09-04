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
rcParams["figure.figsize"] = (9, 9)
rcParams["font.size"] = 13

tracked_csv = "tracked_stars.csv"

sns.set(style="whitegrid")

def plot_quadrant_with_stars(fits_path, stars_df, output_path, quadrant_name, cmap, target, date):
    """
    Plot a quadrant FITS image with circles around stars, colored by tracking percentage.
    Color bar height matches the image height.
    """
    # Load FITS file
    try:
        hdu = fits.open(fits_path)
        img_data = hdu[0].data
        img_wcs = WCS(hdu[0].header)
    except Exception as e:
        print(f"[!] Error loading FITS file {fits_path}: {e}")
        return

    # Create figure with GridSpec for precise layout
    fig = plt.figure(figsize=(18, 17))  # Adjusted width to accommodate color bar
    gs = GridSpec(1, 2, width_ratios=[16, 1], wspace=0.05)  # Image:Colorbar ratio

    # Create axes for image with WCS projection
    ax = fig.add_subplot(gs[0], projection=img_wcs)
    im = ax.imshow(img_data, origin='lower', cmap='gray', vmin=np.percentile(img_data, 5), vmax=np.percentile(img_data, 99))

    # Plot circles for each star
    for _, row in stars_df.iterrows():
        if pd.isna(row['ALPHA_J2000']) or pd.isna(row['DELTA_J2000']):
            continue  # Skip if coordinates are NaN
        sky_coord = SkyCoord(ra=row['ALPHA_J2000'] * u.deg, dec=row['DELTA_J2000'] * u.deg)
        pixel_x, pixel_y = img_wcs.world_to_pixel(sky_coord)

        # Get tracking percentage and map to color
        tracking_percentage = row['tracking_percentage']
        color = cmap(tracking_percentage / 100)  # Normalize to [0, 1] for colormap

        # Draw circle
        circle = Circle((pixel_x, pixel_y), radius=8, edgecolor=color, facecolor='none', linewidth=1)
        ax.add_patch(circle)

    # Create axes for color bar, aligned with image height
    cax = fig.add_subplot(gs[1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(80, 100))
    cbar = fig.colorbar(sm, cax=cax, label='Tracking Percentage (%)')
    cbar.set_ticks(np.linspace(0, 100, 21))  # Show 0%, 20%, 40%, 60%, 80%, 100%

    ax.set_title(f"{quadrant_name}")
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    # plt.tight_layout()
    save_path = os.path.join(output_path, f"{quadrant_name}_{target}_{date}.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()

def summarize_tracking_percentages(df):
    """
    Print a summary of tracking percentages and NaN diagnostics.
    """
    tracking_percentages = df.groupby('star_index')['tracking_percentage'].first()
    print("\n===== Tracking Percentage Summary =====")
    print(f"Total stars: {len(tracking_percentages)}")
    print(f"Mean tracking percentage: {tracking_percentages.mean():.1f}%")
    print(f"Median tracking percentage: {tracking_percentages.median():.1f}%")
    print(f"Min tracking percentage: {tracking_percentages.min():.1f}%")
    print(f"Max tracking percentage: {tracking_percentages.max():.1f}%")
    print("\nPercentage distribution:")
    bins = np.histogram(tracking_percentages, bins=[0, 20, 40, 60, 80, 100])[0]
    for i, count in enumerate(bins):
        print(f"{i*20}-{(i+1)*20}%: {count} stars")
    print("\nNaN Diagnostics:")
    nan_quadrants = df['QUADRANT_FILE'].isna().sum()
    print(f"Rows with NaN QUADRANT_FILE: {nan_quadrants} ({100 * nan_quadrants / len(df):.1f}%)")
    print("=======================================\n")

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
    Main function to generate quadrant plots for the first image's 16 quadrants, with stars colored by tracking percentage.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Load CSV
    df = pd.read_csv(os.path.join(csv_path, tracked_csv))

    # Calculate total number of unique images (based on julian_date)
    total_images = len(df['julian_date'].unique())  # or df['image_id'].unique()

    # Calculate tracking percentage for each star
    non_nan_counts = df.groupby('star_index').apply(
        lambda x: x[['ALPHA_J2000', 'DELTA_J2000']].notna().all(axis=1).sum(),
        include_groups=False
    )
    df = df.merge(
        non_nan_counts.rename('tracking_count').to_frame(),
        left_on='star_index',
        right_index=True
    )
    df['tracking_percentage'] = (df['tracking_count'] / total_images) * 100

    # Print tracking percentage summary
    summarize_tracking_percentages(df)

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

    # Use viridis colormap for tracking percentage
    # cmap = plt.cm.viridis
    cmap = plt.cm.hsv_r

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
        plot_quadrant_with_stars(fits_path, stars_in_quadrant, output_dir, quadrant_name, cmap, target, date)

    print(f"[✓] All valid quadrant plots generated for {first_image_id}.")

if __name__ == "__main__":
    
    csv_dir = sys.argv[1]          # Folder containing 'tracked_stars_with_simbad.csv'
    fits_dir = sys.argv[2]         # Folder with WCS-calibrated quadrant FITS files
    output_dir = sys.argv[3]       # Folder to save plots
    target = sys.argv[4]           # Target being observed
    date = sys.argv[5]             # Date observed

    main(csv_path=csv_dir, fits_folder=fits_dir, output_folder=output_dir, target=target, date=date)