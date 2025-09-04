# star_track_video_fixed_wcsaxes.py
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
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

# --- Command-line arguments ---
data_dir = sys.argv[1]
output_dir = sys.argv[2]
target_id = sys.argv[3]
date_code = sys.argv[4]
star_index = int(sys.argv[5])
fits_folder = sys.argv[6]
padding = 150
vignet_shape = (20, 20)
circle_radius = 4
fps = 4

if not os.path.exists(data_dir) or not os.path.exists(fits_folder):
    raise FileNotFoundError(">>> Data or FITS folder does not exist")

# --- Input file paths ---
csv_path = os.path.join(data_dir, f"tracked_stars.csv")
df = pd.read_csv(csv_path)
df['star_index'] = df['star_index'].astype(int)
df = df[df['star_index'] == star_index].sort_values("julian_date")

# --- Filter out rows with NaN RA/Dec ---
df = df.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])
if df.empty:
    raise RuntimeError(">>> No valid coordinates for target star.")

image_ids = df["image_id"].values
times = df["julian_date"].values
a_image = df["A_IMAGE"].values
b_image = df["B_IMAGE"].values

# --- Function to extract VIGNET from cutout ---
def extract_vignet(cutout_data, x_c, y_c, vignet_shape=(20, 20)):
    height, width = vignet_shape
    half_h, half_w = height // 2, width // 2
    y_min = int(y_c - half_h)
    y_max = int(y_c + half_h)
    x_min = int(x_c - half_w)
    x_max = int(x_c + half_w)

    if (y_min < 0 or x_min < 0 or
        y_max > cutout_data.shape[0] or
        x_max > cutout_data.shape[1]):
        return None
    return cutout_data[y_min:y_max, x_min:x_max]

# --- Function to compute radial profile ---
def compute_radial_profile(vignet, center=None):
    if vignet is None:
        return None, None
    height, width = vignet.shape
    if center is None:
        center = ((width - 1) / 2, (height - 1) / 2)
    y, x = np.indices((height, width))
    distances = np.sqrt((x - center[0])**2 + (y - center[1])**2).flatten()
    intensities = vignet.flatten()
    return distances, intensities

# --- Preload cutouts and VIGNETs ---
cutouts = []
vignets = []
valid_times = []
metadata = []
for img_id, jd in tqdm(zip(image_ids, times), total=len(image_ids), desc="Processing frames"):
    row = df[df['image_id'] == img_id].iloc[0]
    target_coord = SkyCoord(row["ALPHA_J2000"] * u.deg, row["DELTA_J2000"] * u.deg)
    base_id = img_id.replace("_stars", "")
    found = False
    for q in range(16):
        fits_path = os.path.join(fits_folder, f"{base_id}_q{q}_wcs.fit")
        if not os.path.exists(fits_path):
            continue
        try:
            with fits.open(fits_path) as hdul:
                data = hdul[0].data
                wcs = WCS(hdul[0].header)
                x, y = wcs.world_to_pixel(target_coord)
                if not (0 <= x < data.shape[1] and 0 <= y < data.shape[0]):
                    continue
                cut = Cutout2D(data, position=(x, y), size=(padding, padding), wcs=wcs)
                x_c, y_c = cut.input_position_cutout
                vignet = extract_vignet(cut.data, x_c, y_c, vignet_shape)
                if vignet is None or vignet.shape != vignet_shape:
                    continue
                cutouts.append((cut.data, cut.wcs))
                vignets.append(vignet)
                valid_times.append(jd)
                metadata.append({
                    'star_index': star_index,
                    'mag_auto': row.get('MAG_AUTO', np.nan),
                    'mag_aper': row.get('MAG_APER', np.nan),
                    'ra': row['ALPHA_J2000'],
                    'dec': row['DELTA_J2000']
                })
                found = True
                break
        except Exception as e:
            print(f">>> Warning: Failed to process {fits_path}: {e}")
            continue
    if not found:
        print(f">>> Warning: No valid FITS file or star out of bounds for image_id {img_id}")

if not cutouts or not vignets:
    raise RuntimeError(">>> No valid frames to create timelapse.")

# --- Global contrast ---
all_data = [cut[0] for cut in cutouts]
vmin = np.percentile(np.concatenate(all_data), 5) if all_data else 0
vmax = np.percentile(np.concatenate(all_data), 99.5) if all_data else 1

all_vignets = np.array(vignets)
vignet_vmin = np.percentile(all_vignets, 5) if all_vignets.size else 0
vignet_vmax = np.percentile(all_vignets, 99.5) if all_vignets.size else 1

# --- Figure ---
fig = plt.figure(figsize=(12, 12))
ax1 = fig.add_subplot(221, projection=cutouts[0][1])  # WCSAxes for RA/Dec
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223, projection='3d')
ax4 = fig.add_subplot(224)

img2 = ax2.imshow(np.zeros(vignet_shape), cmap='viridis', origin='lower')
time_text = ax1.text(0.5, 1.05, '', transform=ax1.transAxes, ha='center')
fig.suptitle("", fontsize=14)

# Labels
ax1.set_xlabel("RA (J2000)")
ax1.set_ylabel("Dec (J2000)")
ax1.set_title('Star Cutout (WCS)')
ax2.set_title('VIGNET Cutout (2D)')
ax3.set_title('VIGNET Surface (3D)')
ax4.set_xlabel('Radius (pixels)')
ax4.set_ylabel('Counts (ADU)')
ax4.set_title('Radial Profile')
ax4.grid(True, linestyle='--', alpha=0.5)
plt.colorbar(img2, ax=ax2, label='Counts (ADU)')

# --- Update function ---
def update(i):
    cut_data, cut_wcs = cutouts[i]
    row = df[df['image_id'] == image_ids[i]].iloc[0]

    # Main cutout (RA/Dec)
    ax1.clear()
    ax1.imshow(cut_data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    ax1.set_title("Star Cutout (WCS)")
    ax1.set_xlabel("RA (J2000)")
    ax1.set_ylabel("Dec (J2000)")

    target_coord = SkyCoord(row["ALPHA_J2000"] * u.deg, row["DELTA_J2000"] * u.deg)
    tx, ty = cut_wcs.world_to_pixel(target_coord)
    ax1.add_patch(Circle((tx, ty), radius=circle_radius, edgecolor='red', facecolor='none', linewidth=1.5))

    # VIGNET 2D
    vignet = vignets[i]
    img2.set_data(vignet)
    img2.set_clim(vignet_vmin, vignet_vmax)

    # VIGNET 3D
    ax3.clear()
    x_grid, y_grid = np.arange(vignet_shape[1]), np.arange(vignet_shape[0])
    X, Y = np.meshgrid(x_grid, y_grid)
    ax3.plot_surface(X, Y, vignet, cmap='viridis', vmin=vignet_vmin, vmax=vignet_vmax)
    ax3.set_title("VIGNET Surface (3D)")

    # Radial profile with per-frame aperture lines
    distances, intensities = compute_radial_profile(vignet)
    ax4.clear()
    ax4.scatter(distances, intensities, color='blue', s=10, label='Pixel Intensities')
    if len(distances) > 4:
        poly_coeffs = np.polyfit(distances, intensities, deg=4)
        poly_func = np.poly1d(poly_coeffs)
        r_fine = np.linspace(0, np.max(distances), 100)
        ax4.plot(r_fine, poly_func(r_fine), color="black")

    ax4.axvline(x=3, color='red', linestyle='--', label='Fixed Aperture')
    ax4.axvline(x=(a_image[i]*6 + b_image[i]*6)/2, color='green', linestyle='--', label='Auto Aperture AVG')
    ax4.axvline(x=a_image[i]*6, color='orange', linestyle='--', label='Major Axis')
    ax4.axvline(x=b_image[i]*6, color='purple', linestyle='--', label='Minor Axis')

    ax4.set_xlabel("Radius (pixels)")
    ax4.set_ylabel("Counts (ADU)")
    ax4.set_title("Radial Profile")
    ax4.grid(True, linestyle="--", alpha=0.5)
    ax4.legend()

    # Time text
    time_text.set_text(f"JD: {valid_times[i]:.5f}")

    # Suptitle
    meta = metadata[i]
    fig.suptitle(
        f"Star {meta['star_index']}: RA={meta['ra']:.6f}, Dec={meta['dec']:.6f}, "
        f"MAG_AUTO={meta['mag_auto']:.2f}, MAG_APER={meta['mag_aper']:.2f}",
        fontsize=14
    )

    return [img2]

# --- Save animation ---
os.makedirs(output_dir, exist_ok=True)
out_path = os.path.join(output_dir, f"centered_timelapse_with_vignet_{target_id}_{date_code}_star{star_index}.mp4")
ani = FuncAnimation(fig, update, frames=len(cutouts), blit=False)
ani.save(out_path, writer=FFMpegWriter(fps=fps))

print(f">>> Saved video to {out_path}")



# # simple_star_video_centered.py
# import os
# import sys
# import numpy as np
# import pandas as pd
# from astropy.io import fits
# from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
# from astropy.nddata import Cutout2D
# import astropy.units as u
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation, FFMpegWriter
# from matplotlib.patches import Circle
# from tqdm import tqdm

# # --- Command-line arguments ---
# data_dir = sys.argv[1]
# output_dir = sys.argv[2]
# target_id = sys.argv[3]
# date_code = sys.argv[4]
# star_index = int(sys.argv[5])
# fits_folder = sys.argv[6]
# padding = 150 # 70
# circle_radius = 4
# fps = 4

# if not os.path.exists(data_dir) or not os.path.exists(fits_folder):
#     raise FileNotFoundError("❌ Data or FITS folder does not exist")

# # --- Input file paths ---
# csv_path = os.path.join(data_dir, f"tracked_stars.csv")
# df = pd.read_csv(csv_path)
# df['star_index'] = df['star_index'].astype(int)
# df = df[df['star_index'] == star_index].sort_values("julian_date")

# # --- Filter out rows with NaN RA/Dec ---
# df = df.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])
# if df.empty:
#     raise RuntimeError("❌ No valid coordinates for target star.")

# image_ids = df["image_id"].values
# times = df["julian_date"].values

# # --- Set up figure ---
# fig, ax = plt.subplots(figsize=(6, 6))
# img = ax.imshow(np.zeros((padding, padding)), origin="lower", cmap="gray")
# time_text = ax.text(0.5, 1.01, '', transform=ax.transAxes, ha='center')
# circle = None

# # --- Preload cutouts ---
# cutouts = []
# valid_times = []
# for img_id, jd in tqdm(zip(image_ids, times), total=len(image_ids), desc="Processing frames"):
#     row = df[df['image_id'] == img_id].iloc[0]
#     target_coord = SkyCoord(row["ALPHA_J2000"] * u.deg, row["DELTA_J2000"] * u.deg)
#     base_id = img_id.replace("_stars", "")
#     found = False
#     for q in range(16):
#         # fits_path = os.path.join(fits_folder, f"{base_id}_wcs.fit")
#         fits_path = os.path.join(fits_folder, f"{base_id}_q{q}_wcs.fits")
#         if not os.path.exists(fits_path):
#             continue
#         try:
#             with fits.open(fits_path) as hdul:
#                 data = hdul[0].data
#                 wcs = WCS(hdul[0].header)
#                 x, y = wcs.world_to_pixel(target_coord)
#                 if not (0 <= x < data.shape[1] and 0 <= y < data.shape[0]):
#                     continue
#                 cut = Cutout2D(data, position=(x, y), size=(padding, padding), wcs=wcs)
#                 x_c, y_c = cut.input_position_cutout
#                 cutouts.append((cut.data, x_c, y_c))
#                 valid_times.append(jd)
#                 found = True
#                 break
#         except Exception as e:
#             print(f"⚠️ Warning: Failed to process {fits_path}: {e}")
#             continue
#     if not found:
#         print(f"⚠️ Warning: No valid FITS file or star out of bounds for image_id {img_id}")

# # --- Check if any valid cutouts were found ---
# if not cutouts:
#     raise RuntimeError("❌ No valid frames to create timelapse.")

# # --- Compute global contrast ---
# all_data = [cut[0] for cut in cutouts]
# vmin = np.percentile(np.concatenate(all_data), 5) if all_data else 0
# vmax = np.percentile(np.concatenate(all_data), 99.5) if all_data else 1

# # --- Update function ---
# def update(i):
#     global circle
#     data, x, y = cutouts[i]
#     img.set_data(data)
#     img.set_clim(vmin, vmax)
#     time_text.set_text(f"JD: {valid_times[i]:.5f}")
#     if circle:
#         circle.remove()
#     if x is not None and y is not None:
#         circle = Circle((x, y), radius=circle_radius, edgecolor='red', facecolor='none', linewidth=1.5)
#         ax.add_patch(circle)
#     return [img, time_text, circle] if circle else [img, time_text]

# # --- Save animation ---
# os.makedirs(output_dir, exist_ok=True)
# out_path = os.path.join(output_dir, f"centered_timelapse_{target_id}_{date_code}_star{star_index}.mp4")
# ani = FuncAnimation(fig, update, frames=len(cutouts), blit=True)
# ani.save(out_path, writer=FFMpegWriter(fps=fps))

# print(f"✅ Saved video to {out_path}")