# star_video_revamped.py – generates 3x2 panel animations with FITS cutout, VIGNET, radial profile, and differential photometry

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u

# --- Command-line arguments ---
data_dir = sys.argv[1]
output_dir = sys.argv[2]
target_id = sys.argv[3]
date_code = sys.argv[4]
star_index = int(sys.argv[5])
fits_folder = sys.argv[6]
w = sys.argv[7]

# --- Filenames ---
csv_path = os.path.join(data_dir, f"tracked_stars_with_differential_mags_{target_id}_{date_code}.csv")
json_path = os.path.join(data_dir, f"comparison_star_map_{target_id}_{date_code}.json")
output_path = os.path.join(output_dir, f"star_timelapse_{target_id}_star{star_index}_{w}.mp4")
os.makedirs(output_dir, exist_ok=True)

# --- Load CSV and JSON ---
df = pd.read_csv(csv_path)
df['star_index'] = df['star_index'].astype(int)
with open(json_path, 'r') as f:
    comp_map = json.load(f)
comparison_stars_indices = comp_map.get(str(star_index), [])

target_df = df[df['star_index'] == star_index].sort_values("julian_date").copy()
if target_df.empty:
    raise RuntimeError("No valid rows for target star.")

# --- Parse VIGNET string into 20x20 array ---
def parse_vignet(vignet_str):
    """Parse VIGNET column string into a 20x20 numpy array (round brackets, ; separators)."""
    if pd.isna(vignet_str) or not isinstance(vignet_str, str):
        return None
    clean = vignet_str.strip("()")
    row_strs = clean.split(");(")
    rows = []
    for r in row_strs:
        r_clean = r.replace("(", "").replace(")", "")
        if not r_clean.strip():
            continue
        row_vals = [float(x) for x in r_clean.split(";") if x.strip()]
        rows.append(row_vals)
    try:
        arr = np.array(rows)
        if arr.shape == (20, 20):
            return arr
        else:
            print(f"Warning: Parsed VIGNET shape {arr.shape} (expected 20x20)")
            return None
    except Exception as e:
        print(f"Error parsing VIGNET: {e}")
        return None

def compute_radial_profile(vignet, center=None):
    if vignet is None:
        return None, None
    h, w = vignet.shape
    if center is None:
        center = ((w - 1) / 2, (h - 1) / 2)
    y, x = np.indices((h, w))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2).flatten()
    vals = vignet.flatten()
    return r, vals

# --- Coordinates ---
row0 = target_df.dropna(subset=["ALPHA_J2000", "DELTA_J2000"]).iloc[0]
target_coord = SkyCoord(row0['ALPHA_J2000'] * u.deg, row0['DELTA_J2000'] * u.deg)
ra_deg = row0['ALPHA_J2000']
dec_deg = row0['DELTA_J2000']

# --- Arrays ---
times = target_df["julian_date"].values
image_ids = target_df["image_id"].values
mag_auto = target_df["MAG_AUTO"].values
mag_aper = target_df["MAG_APER"].values
mag_auto_err = target_df["MAGERR_AUTO"].values
mag_aper_err = target_df["MAGERR_APER"].values
diff_mag = target_df["differential_mag"].values
diff_mag_err = target_df["differential_mag_err"].values
weighted_diff_mag = target_df["weighted_differential_mag"].values
weighted_diff_mag_err = target_df["weighted_differential_mag_err"].values
vignets = [parse_vignet(v) for v in target_df["VIGNET"].values]

# --- Preload FITS cutouts (for WCS panel) ---
padding = 100
preload_cutouts = []
for img_id in image_ids:
    base_id = img_id.replace("_stars", "")
    found = False
    for q in range(16):
        fpath = os.path.join(fits_folder, f"{base_id}_q{q}_wcs.fit")
        if not os.path.exists(fpath):
            continue
        try:
            with fits.open(fpath) as hdul:
                data = hdul[0].data
                wcs = WCS(hdul[0].header)
                x, y = wcs.world_to_pixel(target_coord)
                if not (0 <= x < data.shape[1] and 0 <= y < data.shape[0]):
                    continue
                cut = Cutout2D(data, position=(x, y), size=(padding, padding), wcs=wcs)
                preload_cutouts.append((cut.data, cut.wcs))
                found = True
                break
        except Exception:
            continue
    if not found:
        preload_cutouts.append((None, None))

# --- Setup figure (3x2 grid) ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
(ax_img, ax_v2d, ax_v3d, ax_radial, ax_lc, ax_err) = axes.flatten()
ax_v3d.remove()
ax_v3d = fig.add_subplot(2, 3, 3, projection='3d')
fig.subplots_adjust(top=0.8)

# --- Update function ---
def update(i):
    cut_data, cut_wcs = preload_cutouts[i]
    vignet = vignets[i]

    # FITS cutout with WCS + target + comps
    ax_img.clear()
    if cut_data is not None:
        # Percentile scaling
        vmin, vmax = np.percentile(cut_data, [5, 99])
        ax_img.imshow(cut_data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ax_img.set_title("Image Cutout (WCS)")
        ax_img.set_xlabel("RA (J2000)")
        ax_img.set_ylabel("Dec (J2000)")
        if cut_wcs:
            tx, ty = cut_wcs.world_to_pixel(target_coord)
            ax_img.add_patch(Circle((tx, ty), 5, edgecolor='red', facecolor='none'))
            frame_df = df[df['image_id'] == image_ids[i]]
            for comp_idx in comparison_stars_indices:
                comp_row = frame_df[frame_df['star_index'] == comp_idx]
                if not comp_row.empty:
                    cx, cy = cut_wcs.world_to_pixel(
                        SkyCoord(comp_row.iloc[0]['ALPHA_J2000'] * u.deg,
                                 comp_row.iloc[0]['DELTA_J2000'] * u.deg))
                    ax_img.add_patch(Circle((cx, cy), 5, edgecolor='cyan', facecolor='none'))


    # Vignet 2D
    ax_v2d.clear()
    if vignet is not None:
        ax_v2d.imshow(vignet, origin="lower", cmap="viridis")
    ax_v2d.set_title("VIGNET 2D")
    ax_v2d.set_xlabel("X (pixels)")
    ax_v2d.set_ylabel("Y (pixels)")

    # Vignet 3D
    ax_v3d.clear()
    if vignet is not None:
        X, Y = np.meshgrid(np.arange(vignet.shape[1]), np.arange(vignet.shape[0]))
        ax_v3d.plot_surface(X, Y, vignet, cmap="viridis")
    ax_v3d.set_title("VIGNET 3D")
    ax_v3d.set_xlabel("X (pixels)")
    ax_v3d.set_ylabel("Y (pixels)")
    ax_v3d.set_zlabel("Counts (ADU)")

    # Radial profile
    ax_radial.clear()
    if vignet is not None:
        r, vals = compute_radial_profile(vignet)
        ax_radial.scatter(r, vals, s=5, color="blue")
        if len(r) > 10:
            coeffs = np.polyfit(r, vals, deg=4)
            rr = np.linspace(0, r.max(), 100)
            ax_radial.plot(rr, np.poly1d(coeffs)(rr), color="black")
    ax_radial.set_title("Radial Profile")
    ax_radial.set_xlabel("Radius (pixels)")
    ax_radial.set_ylabel("Counts (ADU)")

    # Differential photometry (error bars + marker)
    ax_lc.clear()
    ax_lc.errorbar(times, diff_mag, yerr=diff_mag_err, fmt='o', color="blue", markersize=2)
    if np.isfinite(diff_mag[i]):
        ax_lc.plot(times[i], diff_mag[i], 'ro')
    ax_lc.invert_yaxis()
    ax_lc.set_title("Differential Light Curve")
    ax_lc.set_xlabel("Julian Date")
    ax_lc.set_ylabel("Differential Magnitude")

    # Error vs time
    ax_err.clear()
    ax_err.plot(times, mag_auto_err, label="MAGERR_AUTO")
    ax_err.plot(times, mag_aper_err, label="MAGERR_APER")
    ax_err.plot(times, diff_mag_err, label="Diff Mag Err")
    ax_err.plot(times, weighted_diff_mag_err, label="Weighted Diff Mag Err")
    ax_err.axvline(times[i], color="red", linestyle="--", linewidth=1)
    ax_err.legend()
    ax_err.set_title("Error vs Time")
    ax_err.set_xlabel("Julian Date")
    ax_err.set_ylabel("Measurement / Calculation Error")

    # Time calculations
    total_hours = (times[-1] - times[0]) * 24
    elapsed_hours = (times[i] - times[0]) * 24

    # Update suptitle
    fig.suptitle(
        f"Target {target_id} ({date_code}) | RA={ra_deg:.5f} Dec={dec_deg:.5f}\n"
        f"Observation Length: {elapsed_hours:.2f}/{total_hours:.2f} Hr | Frame {i+1}/{len(image_ids)}\n"
        f"MAG_APER={mag_aper[i]:.2f}, "
        f"MAG_AUTO={mag_auto[i]:.2f}, "
        f"Diff Mag={diff_mag[i]:.2f}, "
        f"WDiff Mag={weighted_diff_mag[i]:.2f}"
    )

# --- Animate ---
ani = FuncAnimation(fig, update, frames=len(image_ids), blit=False)
ani.save(output_path, writer=FFMpegWriter(fps=5))
print(f">>> Saved to {output_path}")







# # star_video.py – generates .mp4 animations for star light curves and image cutouts

# import os
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle
# from matplotlib.animation import FuncAnimation, FFMpegWriter
# from astropy.io import fits
# from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
# import astropy.units as u
# import json
# import sys
# from astropy.nddata import Cutout2D
# from astropy.visualization.wcsaxes import WCSAxes

# # --- Command-line arguments ---
# data_dir = sys.argv[1]
# output_dir = sys.argv[2]
# target_id = sys.argv[3]
# date_code = sys.argv[4]
# star_index = int(sys.argv[5])
# fits_folder = sys.argv[6]
# w = sys.argv[9]

# # --- Filenames ---
# csv_path = os.path.join(data_dir, f"tracked_stars_with_differential_mags_{target_id}_{date_code}.csv")
# json_path = os.path.join(data_dir, f"comparison_star_map_{target_id}_{date_code}.json")
# output_path = os.path.join(output_dir, f"star_timelapse_{target_id}_star{star_index}_{w}.mp4")

# os.makedirs(output_dir, exist_ok=True)

# # --- Load CSV and JSON ---
# df = pd.read_csv(csv_path)
# df['star_index'] = df['star_index'].astype(int)

# with open(json_path, 'r') as f:
#     comp_map = json.load(f)

# comparison_stars_indices = comp_map.get(str(star_index), [])

# # Use the full target_df without dropping NaNs initially for plotting gaps
# target_df_all = df[df['star_index'] == star_index].sort_values("julian_date").copy()

# # Determine a constant target_skycoord for WCS transformations
# initial_target_coords_row = target_df_all.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])
# if initial_target_coords_row.empty:
#     raise RuntimeError(">>> Target star has no valid celestial coordinates (ALPHA_J2000, DELTA_J2000) in the CSV.")
# target_skycoord_constant = SkyCoord(initial_target_coords_row.iloc[0]['ALPHA_J2000'] * u.deg, 
#                                     initial_target_coords_row.iloc[0]['DELTA_J2000'] * u.deg)

# # --- Extract all data for light curve plotting, including NaNs ---
# times = target_df_all["julian_date"].values
# diff_mags = target_df_all[sys.argv[7]].values
# diff_errs = target_df_all[sys.argv[8]].values
# image_ids = target_df_all["image_id"].values

# # --- Preload FITS cutouts ---
# preload_frame_data = [] # Stores (cut_data, cut_wcs, image_id, vmin, vmax, cutout_dim, target_mag_is_nan)
# padding_pixels = 20 # Padding around the bounding box of stars for cutout

# print("Loading FITS cutouts...")
# initial_last_valid_info = None # Stores (cut_data, cut_wcs, vmin, vmax, cutout_dim)
# for i in range(len(image_ids)):
#     image_id_current = image_ids[i]
#     target_mag_is_nan_for_frame = np.isnan(diff_mags[i])

#     base_id = image_id_current.replace("_stars","")

#     current_frame_fits_info = None

#     current_image_df_subset = df[df['image_id'] == image_id_current]
#     star_coords_for_cutout = [target_skycoord_constant] 
#     for comp_idx in comparison_stars_indices:
#         comp_row = current_image_df_subset[current_image_df_subset['star_index'] == comp_idx]
#         if not comp_row.empty and np.isfinite(comp_row.iloc[0]["ALPHA_J2000"]) and np.isfinite(comp_row.iloc[0]["DELTA_J2000"]):
#             star_coords_for_cutout.append(SkyCoord(comp_row.iloc[0]["ALPHA_J2000"]*u.deg, comp_row.iloc[0]["DELTA_J2000"]*u.deg))

#     if not star_coords_for_cutout:
#         preload_frame_data.append((None, None, image_id_current, None, None, None, target_mag_is_nan_for_frame))
#         print(f"  No valid star coordinates for cutout calculation for image_id: {image_id_current}. Marking FITS data as invalid.")
#         continue

#     found_valid_quadrant_for_frame = False
#     for q in range(16):
#         fits_file = os.path.join(fits_folder, f"{base_id}_q{q}_wcs.fit")
#         if not os.path.exists(fits_file):
#             continue
#         try:
#             with fits.open(fits_file) as hdul:
#                 data = hdul[0].data
#                 wcs = WCS(hdul[0].header)

#                 star_pixel_coords_in_quadrant = []
#                 all_stars_in_this_quadrant = True
#                 for star_coord in star_coords_for_cutout:
#                     pix_x, pix_y = wcs.world_to_pixel(star_coord)
#                     if not (-5 <= pix_x < data.shape[1] + 5 and -5 <= pix_y < data.shape[0] + 5):
#                         all_stars_in_this_quadrant = False
#                         break
#                     star_pixel_coords_in_quadrant.append((pix_x, pix_y))
                
#                 if not all_stars_in_this_quadrant:
#                     continue

#                 min_x = np.floor(min([p[0] for p in star_pixel_coords_in_quadrant]) - padding_pixels).astype(int)
#                 max_x = np.ceil(max([p[0] for p in star_pixel_coords_in_quadrant]) + padding_pixels).astype(int)
#                 min_y = np.floor(min([p[1] for p in star_pixel_coords_in_quadrant]) - padding_pixels).astype(int)
#                 max_y = np.ceil(max([p[1] for p in star_pixel_coords_in_quadrant]) + padding_pixels).astype(int)

#                 width = max_x - min_x
#                 height = max_y - min_y
#                 cutout_dim = max(width, height)
                
#                 cutout_center_x_orig = min_x + cutout_dim / 2.0
#                 cutout_center_y_orig = min_y + cutout_dim / 2.0

#                 cutout = Cutout2D(data, position=(cutout_center_x_orig, cutout_center_y_orig),
#                                   size=(cutout_dim, cutout_dim), wcs=wcs)
                
#                 cut_data = cutout.data
#                 cut_wcs = cutout.wcs

#                 finite_cut_data = cut_data[np.isfinite(cut_data)]
#                 if finite_cut_data.size == 0 or np.all(finite_cut_data == 0) or (finite_cut_data.size > 0 and np.std(finite_cut_data) < 1e-6):
#                     continue
                
#                 vmin_frame = np.nanpercentile(finite_cut_data, 10)
#                 vmax_frame = np.nanpercentile(finite_cut_data, 99.5)

#                 current_frame_fits_info = (cut_data, cut_wcs, image_id_current, vmin_frame, vmax_frame, cutout_dim)
#                 found_valid_quadrant_for_frame = True
#                 print(f"  Loaded quadrant {q} for image_id: {image_id_current}")
#                 if initial_last_valid_info is None:
#                     initial_last_valid_info = (cut_data, cut_wcs, vmin_frame, vmax_frame, cutout_dim)
#                 break 
#         except Exception as e:
#             pass

#     if current_frame_fits_info:
#         preload_frame_data.append(current_frame_fits_info + (target_mag_is_nan_for_frame,))
#     else:
#         preload_frame_data.append((None, None, image_id_current, None, None, None, target_mag_is_nan_for_frame))
#         print(f"  No valid FITS data found for image_id: {image_id_current}. This frame will display previous valid image.")


# if initial_last_valid_info is None:
#     raise RuntimeError(">>> No valid FITS data found for any frame. Cannot start animation. Check FITS files and star coordinates.")

# initial_cut_data, initial_cut_wcs, initial_vmin, initial_vmax, initial_cutout_dim = initial_last_valid_info

# # --- Calculate Observation Duration ---
# if len(times) > 1 and np.isfinite(times[0]) and np.isfinite(times[-1]):
#     observation_start_jd = times[0]
#     total_observation_duration_days = times[-1] - observation_start_jd
#     total_observation_duration_hours = total_observation_duration_days * 24
# else:
#     observation_start_jd = times[0] if len(times) > 0 and np.isfinite(times[0]) else 0
#     total_observation_duration_hours = 0


# # --- Set up figure ---
# # Increased figure height slightly to make more room for text below plots
# fig = plt.figure(figsize=(14, 8)) #constrained_layout=True)

# # Add an overall title to the figure
# fig.suptitle(f"Light Curve and Image Cutout for Target Star {target_id} ({date_code}) \n Observation Length: {total_observation_duration_hours:.2f} Hr", fontsize=16)

# # Modify gridspec to allocate more space at the bottom for the text
# # 'bottom' parameter defines the lower boundary of the grid in figure coordinates (0 to 1)
# gs = fig.add_gridspec(1, 2, bottom=0.15) # Increased from default to give space for fig.text

# # Image subplot setup using WCSAxes
# img_ax = fig.add_subplot(gs[0, 0], projection=initial_cut_wcs)
# img = img_ax.imshow(initial_cut_data, origin="lower", cmap="gray", vmin=initial_vmin, vmax=initial_vmax)

# img_ax.set_title("Image Cutout")
# # img_ax.grid(True, linestyle='--', alpha=0.9, color='white')

# # Configure the WCSAxes coordinate system initially
# ra_axis = img_ax.coords[0]
# dec_axis = img_ax.coords[1]

# ra_axis.set_axislabel("Right Ascension (J2000)")
# dec_axis.set_axislabel("Declination (J2000)")

# ra_axis.set_major_formatter('hh:mm:ss')
# dec_axis.set_major_formatter('dd:mm:ss')

# ra_axis.set_ticklabel(exclude_overlapping=True)
# dec_axis.set_ticklabel(exclude_overlapping=True)

# img_ax.coords.frame.set_linewidth(1.5)
# img_ax.coords[0].set_ticks(color='white')
# img_ax.coords[1].set_ticks(color='white')


# # --- Add warning boxes ---
# no_fits_indicator_text = img_ax.text(0.95, 0.95, 'NO IMAGE DATA (Prev. Frame)',
#                                  horizontalalignment='right', verticalalignment='top',
#                                  transform=img_ax.transAxes,
#                                  bbox=dict(boxstyle='round,pad=0.5', fc='red', alpha=0.7, ec='none'),
#                                  color='white', fontsize=10, zorder=10)
# no_fits_indicator_text.set_visible(False)

# nan_mag_indicator_text = img_ax.text(0.05, 0.95, 'TARGET MAG. NaN',
#                                       horizontalalignment='left', verticalalignment='top',
#                                       transform=img_ax.transAxes,
#                                       bbox=dict(boxstyle='round,pad=0.5', fc='orange', alpha=0.7, ec='none'),
#                                       color='white', fontsize=10, zorder=10)
# nan_mag_indicator_text.set_visible(False)

# # Light curve subplot setup
# lc_ax = fig.add_subplot(gs[0, 1])
# lc_ax.plot(times, diff_mags, 'o-', color='deepskyblue', markersize=2)
# smooth = pd.Series(diff_mags).rolling(window=5, center=True, min_periods=1).median()
# lc_ax.plot(times, smooth, color='black', label='5-Frame Median')
# lc_ax.fill_between(times, diff_mags - diff_errs, diff_mags + diff_errs, color='gray', alpha=0.3)
# frame_marker, = lc_ax.plot([], [], 'ro', markersize=6)

# lc_ax.set_xlabel("Julian Date")
# lc_ax.set_ylabel("Differential Magnitude")
# lc_ax.invert_yaxis()
# lc_ax.set_title("Differential Light Curve")
# # lc_ax.grid(True, linestyle='--', alpha=0.9)

# # --- Set fixed y-axis limits for the light curve to prevent jumping, handling NaNs
# valid_lc_data_mask = np.isfinite(diff_mags) & np.isfinite(diff_errs)
# finite_diff_mags_filtered = diff_mags[valid_lc_data_mask]
# finite_diff_errs_filtered = diff_errs[valid_lc_data_mask]

# if finite_diff_mags_filtered.size > 0:
#     lc_y_min = np.nanmin(finite_diff_mags_filtered - finite_diff_errs_filtered)
#     lc_y_max = np.nanmax(finite_diff_mags_filtered + finite_diff_errs_filtered)
#     padding_lc = (lc_y_max - lc_y_min) * 0.1
#     lc_ax.set_ylim(lc_y_max + padding_lc, lc_y_min - padding_lc)
# else:
#     lc_ax.set_ylim(1.0, 0.0)

# finite_times = times[np.isfinite(times)]
# if finite_times.size > 0:
#     lc_ax.set_xlim(np.nanmin(finite_times), np.nanmax(finite_times))
# else:
#     lc_ax.set_xlim(0, 1)

# # Add text elements below the plots
# # Adjusted 'y' coordinates to fit within the new bottom margin
# text_image_info = fig.text(0.5, 0.08, '', ha='center', va='bottom', fontsize=12) # Moved up
# text_time_info = fig.text(0.5, 0.04, '', ha='center', va='bottom', fontsize=12) # Moved up

# # --- Global variable to store the last valid cutout information ---
# current_display_fits_info = (initial_cut_data, initial_cut_wcs, initial_cutout_dim, initial_vmin, initial_vmax)

# # --- Update function ---
# def update(frame):
#     global current_display_fits_info # Declare intent to modify global variable
    
#     current_frame_data_tuple = preload_frame_data[frame]
    
#     cut_data_for_frame = current_frame_data_tuple[0]
#     cut_wcs_for_frame = current_frame_data_tuple[1]
#     image_id_current_frame = current_frame_data_tuple[2]
#     vmin_frame = current_frame_data_tuple[3]
#     vmax_frame = current_frame_data_tuple[4]
#     cutout_dim_for_frame = current_frame_data_tuple[5]
#     target_mag_is_nan_for_frame = current_frame_data_tuple[6]

#     current_frame_has_fits_data = (cut_data_for_frame is not None)

#     display_cut_data = None
#     display_cut_wcs = None
#     display_cutout_dim = None
#     display_vmin = None
#     display_vmax = None

#     if current_frame_has_fits_data:
#         display_cut_data = cut_data_for_frame
#         display_cut_wcs = cut_wcs_for_frame
#         display_cutout_dim = cutout_dim_for_frame
#         display_vmin = vmin_frame
#         display_vmax = vmax_frame

#         current_display_fits_info = (display_cut_data, display_cut_wcs, display_cutout_dim, display_vmin, display_vmax)
#         no_fits_indicator_text.set_visible(False)
#     else:
#         display_cut_data, display_cut_wcs, display_cutout_dim, display_vmin, display_vmax = current_display_fits_info
#         no_fits_indicator_text.set_visible(True)
#         no_fits_indicator_text.set_text('NO IMAGE DATA (Prev. Frame)')

#     # Set image data and color limits
#     img.set_data(display_cut_data)
#     img.set_clim(vmin=display_vmin, vmax=display_vmax)
    
#     # Update the WCSAxes instance's WCS using reset_wcs.
#     # After resetting, re-apply the desired axis labels, formatters, and tick properties.
#     img_ax.reset_wcs(wcs=display_cut_wcs) 

#     # Re-apply axis labels and formatters (they might be reset by reset_wcs)
#     img_ax.coords[0].set_axislabel("Right Ascension (J2000)")
#     img_ax.coords[1].set_axislabel("Declination (J2000)")
#     img_ax.coords[0].set_major_formatter('hh:mm:ss')
#     img_ax.coords[1].set_major_formatter('dd:mm:ss')

#     # Ensure ticks are shown and labels don't overlap (might also be affected by reset_wcs)
#     img_ax.coords[0].set_ticklabel(exclude_overlapping=True)
#     img_ax.coords[1].set_ticklabel(exclude_overlapping=True)
#     img_ax.coords[0].set_ticks(color='white')
#     img_ax.coords[1].set_ticks(color='white')
#     img_ax.coords.frame.set_linewidth(1.5) # Reapply frame linewidth
#     # img_ax.coords.grid(True, linestyle='--', alpha=0.9, color='white') # Reapply grid

#     # Update the image's extent to match the new cutout dimensions
#     img.set_extent([-0.5, display_cutout_dim - 0.5, -0.5, display_cutout_dim - 0.5])
#     # Also reset the view limits of the WCSAxes to match the new cutout.
#     img_ax.set_xlim(-0.5, display_cutout_dim - 0.5)
#     img_ax.set_ylim(-0.5, display_cutout_dim - 0.5)

#     # --- Update warning boxes visibility ---
#     nan_mag_indicator_text.set_visible(target_mag_is_nan_for_frame)
#     if target_mag_is_nan_for_frame:
#         nan_mag_indicator_text.set_text("TARGET MAG. NaN")


#     # --- Plot Circles ---
#     # Remove ALL old patches before adding new ones
#     for patch in list(img_ax.patches):
#         patch.remove()

#     target_pixel_in_cutout = display_cut_wcs.world_to_pixel(target_skycoord_constant)
#     target_circle_current = Circle((target_pixel_in_cutout[0], target_pixel_in_cutout[1]), 5, edgecolor='red', facecolor='none', linewidth=1.5)
#     img_ax.add_patch(target_circle_current)

#     current_image_df_subset_for_comp_coords = df[df['image_id'] == image_ids[frame]]
#     for comp_idx in comparison_stars_indices:
#         comp_row = current_image_df_subset_for_comp_coords[current_image_df_subset_for_comp_coords['star_index'] == comp_idx]
#         if not comp_row.empty and np.isfinite(comp_row.iloc[0]["ALPHA_J2000"]) and np.isfinite(comp_row.iloc[0]["DELTA_J2000"]):
#             comp_ra = comp_row.iloc[0]["ALPHA_J2000"]
#             comp_dec = comp_row.iloc[0]["DELTA_J2000"]
#             try:
#                 comp_coord = SkyCoord(comp_ra*u.deg, comp_dec*u.deg)
#                 comp_pixel_cutout = display_cut_wcs.world_to_pixel(comp_coord)
                
#                 if (-padding_pixels <= comp_pixel_cutout[0] < display_cutout_dim + padding_pixels and
#                     -padding_pixels <= comp_pixel_cutout[1] < display_cutout_dim + padding_pixels):
#                     circle = Circle((comp_pixel_cutout[0], comp_pixel_cutout[1]), 5, edgecolor='cyan', facecolor='none', linewidth=1.2)
#                     img_ax.add_patch(circle)
#             except Exception as e:
#                 pass

#     # --- Update Light Curve Dot ---
#     if np.isfinite(diff_mags[frame]) and np.isfinite(times[frame]):
#         frame_marker.set_data([times[frame]], [diff_mags[frame]])
#     else:
#         frame_marker.set_data([], [])

#     # --- Update Text Information ---
#     current_image_number = frame + 1
#     text_image_info.set_text(f"Image {current_image_number}/{len(image_ids)}")
#     current_time_jd = times[frame]
#     time_into_observation_days = current_time_jd - observation_start_jd
#     time_into_observation_hours = time_into_observation_days * 24
#     text_time_info.set_text(f"Time into Observation: {time_into_observation_hours:.2f} Hours")
    
#     # Return artists that might have changed. With blit=False, Matplotlib redraws
#     # the entire figure anyway, so explicitly listing all WCSAxes internal artists
#     # is not necessary and was causing the error.
#     return [img, frame_marker, text_image_info, text_time_info, no_fits_indicator_text, nan_mag_indicator_text] + list(img_ax.patches)

# # --- Animate and save ---
# print(f"Saving animation to {output_path}...")
# ani = FuncAnimation(fig, update, frames=len(preload_frame_data), blit=False)
# try:
#     ani.save(output_path, writer=FFMpegWriter(fps=5))
#     print(">>> Animation saved successfully!")
# except Exception as e:
#     print(f">>> Error saving animation: {e}")
#     print("Ensure FFmpeg is installed and accessible in your system's PATH.")