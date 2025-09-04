import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.patches import Circle
from collections import defaultdict
import numpy as np

from matplotlib import rcParams
rcParams["figure.figsize"] = (8, 8)
rcParams["font.size"] = 10

import seaborn as sns
sns.set(style="whitegrid")

def plot_quadrant_with_stars(fits_path, stars_df, output_path, quadrant_name, color_map, target, date):
    hdu = fits.open(fits_path)
    img_data = hdu[0].data
    img_wcs = WCS(hdu[0].header)

    fig = plt.figure()
    ax = plt.subplot(projection=img_wcs)
    ax.imshow(img_data, origin='lower', cmap='gray', vmin=np.percentile(img_data, 5), vmax=np.percentile(img_data, 99))

    legend_labels = {}

    for _, row in stars_df.iterrows():
        sky_coord = SkyCoord(ra=row['ALPHA_J2000'] * u.deg, dec=row['DELTA_J2000'] * u.deg)
        pixel_x, pixel_y = img_wcs.world_to_pixel(sky_coord)

        star_type = row['Simbad_Type']
        color = color_map.get(star_type, 'gray')

        circle = Circle((pixel_x, pixel_y), radius=10, edgecolor=color, facecolor='none', linewidth=1.5)
        ax.add_patch(circle)

        if star_type not in legend_labels:
            legend_labels[star_type] = color

    legend_patches = [Circle((0, 0), radius=5, edgecolor=color, facecolor='none', label=label) for label, color in legend_labels.items()]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=8, title="Simbad Types")

    ax.set_title(f"{quadrant_name}")
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    plt.tight_layout()
    save_path = os.path.join(output_path, f"{quadrant_name}_{target}_{date}.png")
    plt.savefig(save_path)
    plt.close()

def summarize_types(df):
    counts = df['Simbad_Type'].value_counts()
    total = counts.sum()
    print("\n===== SIMBAD Type Summary =====")
    for obj_type, count in counts.items():
        percent = 100 * count / total
        print(f"{obj_type:<20}: {count:4} ({percent:.1f}%)")
    print("===============================\n")

def main(csv_path, fits_folder, output_folder, target, date):
    os.makedirs(output_folder, exist_ok=True)

    df = pd.read_csv(os.path.join(csv_path, "tracked_stars_with_simbad.csv"))

    # Print SIMBAD type summary
    summarize_types(df)

    # Build color map
    unique_types = df['Simbad_Type'].unique()
    palette = sns.color_palette("husl", len(unique_types))
    color_map = dict(zip(unique_types, palette))

    # Process each quadrant
    for quadrant_file in df['QUADRANT_FILE'].unique():
        quadrant_name = quadrant_file.replace("_se.ldac", "")
        fits_name = quadrant_name + ".fit"
        fits_path = os.path.join(fits_folder, fits_name)

        if not os.path.exists(fits_path):
            print(f">>> Warning: FITS file not found: {fits_path}")
            continue

        stars_in_quadrant = df[df['QUADRANT_FILE'] == quadrant_file]

        print(f">>> Plotting {quadrant_name} with {len(stars_in_quadrant)} stars...")
        plot_quadrant_with_stars(fits_path, stars_in_quadrant, output_folder, quadrant_name, color_map, target, date)

    print(">>> All quadrant plots generated.")

if __name__ == "__main__":
    csv_dir = sys.argv[1]          # Folder containing 'tracked_stars_with_simbad.csv'
    fits_dir = sys.argv[2]         # Folder with WCS-calibrated quadrant FITS files
    output_dir = sys.argv[3]       # Folder to save plots
    target = sys.argv[4]           # Target being observed
    date = sys.argv[5]             # date observed 

    main(csv_dir, fits_dir, output_dir, target, date)
