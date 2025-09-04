import os
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import Circle
from tqdm import tqdm

# --- Command-line arguments ---
data_dir = sys.argv[1]
output_dir = sys.argv[2]
target_id = sys.argv[3]
date_code = sys.argv[4]
star_index = int(sys.argv[5])
fits_folder = sys.argv[6]
padding = 150  # Size of cutout in pixels
circle_radius = 4  # Radius of circle in pixels
fps = 4  # Frames per second for video

# --- Validate input directories ---
if not os.path.exists(data_dir) or not os.path.exists(fits_folder):
    raise FileNotFoundError("❌ Data or FITS folder does not exist")

# --- Load and filter tracked stars data ---
csv_path = os.path.join(data_dir, f"tracked_stars.csv")
df = pd.read_csv(csv_path)
df['star_index'] = df['star_index'].astype(int)
df = df[df['star_index'] == star_index].sort_values("julian_date")

# --- Filter out rows with NaN RA/Dec ---
df = df.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])
if df.empty:
    raise RuntimeError("❌ No valid coordinates for target star.")

# image_ids = df["image_id"].values
image_ids = df["QUADRANT_FILE"].values
times = df["julian_date"].values
ra_values = df["ALPHA_J2000"].values
dec_values = df["DELTA_J2000"].values

# --- Compute mean RA/Dec as fixed cutout center ---
mean_coord = SkyCoord(ra=np.mean(ra_values) * u.deg, dec=np.mean(dec_values) * u.deg)
print(f"[⧗] Using mean coordinates for fixed cutout: RA={mean_coord.ra.deg:.6f}, Dec={mean_coord.dec.deg:.6f}")

# --- Preload cutouts ---
cutouts = []
valid_times = []
star_positions = []  # Store (x, y) pixel positions of target star
for img_id, jd, ra, dec in tqdm(zip(image_ids, times, ra_values, dec_values), total=len(image_ids), desc="Processing frames"):
    target_coord = SkyCoord(ra * u.deg, dec * u.deg)
    # base_id = img_id.replace("_stars", "")
    base_id = img_id.replace("_se.ldac", "")
    found = False
    # fits_path = os.path.join(fits_folder, f"{base_id}_wcs.fit")
    fits_path = os.path.join(fits_folder, f"{base_id}.fits")
    if not os.path.exists(fits_path):
        print(f"⚠️ Warning: FITS file not found for {fits_path}")
        continue
    try:
        with fits.open(fits_path) as hdul:
            data = hdul[0].data
            wcs = WCS(hdul[0].header)
            # Create cutout centered on mean_coord (fixed across all frames)
            cut = Cutout2D(data, position=mean_coord, size=(padding, padding), wcs=wcs)
            # Compute target star’s pixel position within the cutout
            x_star, y_star = cut.wcs.world_to_pixel(target_coord)
            # Check if the cutout is valid (data exists and star is within bounds)
            if not (0 <= x_star < padding and 0 <= y_star < padding):
                print(f"⚠️ Warning: Star out of cutout bounds for image_id {img_id}")
                continue
            cutouts.append((cut.data, x_star, y_star))
            valid_times.append(jd)
            star_positions.append((x_star, y_star))
            found = True
    except Exception as e:
        print(f"⚠️ Warning: Failed to process {fits_path}: {e}")
        continue
    if not found:
        print(f"⚠️ Warning: No valid FITS file for image_id {img_id}")

# --- Check if any valid cutouts were found ---
if not cutouts:
    raise RuntimeError("❌ No valid frames to create timelapse.")

# --- Compute global contrast ---
all_data = [cut[0] for cut in cutouts]
vmin = np.percentile(np.concatenate(all_data), 5) if all_data else 0
vmax = np.percentile(np.concatenate(all_data), 99.5) if all_data else 1

# --- Set up figure ---
fig, ax = plt.subplots(figsize=(6, 6))
img = ax.imshow(np.zeros((padding, padding)), origin="lower", cmap="gray")
time_text = ax.text(0.5, 1.01, '', transform=ax.transAxes, ha='center')
circle = None

# --- Update function ---
def update(i):
    global circle
    data, x, y = cutouts[i]
    img.set_data(data)
    img.set_clim(vmin, vmax)
    time_text.set_text(f"JD: {valid_times[i]:.5f}")
    if circle:
        circle.remove()
    if x is not None and y is not None:
        circle = Circle((x, y), radius=circle_radius, edgecolor='red', facecolor='none', linewidth=1.5)
        ax.add_patch(circle)
    return [img, time_text, circle] if circle else [img, time_text]

# --- Save animation ---
os.makedirs(output_dir, exist_ok=True)
out_path = os.path.join(output_dir, f"fixed_field_timelapse_{target_id}_{date_code}_star{star_index}.mp4")
ani = FuncAnimation(fig, update, frames=len(cutouts), blit=True)
ani.save(out_path, writer=FFMpegWriter(fps=fps))

print(f"✅ Saved video to {out_path}")