import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.signal import savgol_filter
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
from matplotlib.patches import Circle
import seaborn as sns
import json
import sys
import os

# Set matplotlib style
from matplotlib import rcParams
rcParams["figure.figsize"] = (18, 12)
rcParams["font.size"] = 10
sns.set(style="whitegrid")

# Parse command line arguments
input_data_folder_arg = sys.argv[1]
output_plots_dir_arg = sys.argv[2]
target_id_arg = sys.argv[3]
date_string_arg = sys.argv[4]
star_index_to_plot_arg = int(sys.argv[5])
fits_folder_arg = sys.argv[6]
DIFF_MAG_COLUMN = sys.argv[7]
DIFF_MAG_ERROR_COLUMN = sys.argv[8]
w = sys.argv[9]

# File paths
input_data_path = os.path.join(
    input_data_folder_arg,
    f"tracked_stars_with_differential_mags_{target_id_arg}_{date_string_arg}.csv"
)
comparison_stars_json_path = os.path.join(
    input_data_folder_arg,
    f"comparison_star_map_{target_id_arg}_{date_string_arg}.json"
)

# Prompt for magnitude type
while True:
    print(">>> Use MAG_AUTO or MAG_APER?")
    print("(1) MAG_AUTO")
    print("(2) MAG_APER\n")
    which_mag = int(input("Answer: "))
    print()
    if which_mag == 1:
        MAG_COLUMN = 'MAG_AUTO'
        MAG_ERROR_COLUMN = 'MAGERR_AUTO'
        break
    if which_mag == 2:
        MAG_COLUMN = 'MAG_APER'
        MAG_ERROR_COLUMN = 'MAGERR_APER'
        break
    else:
        print("Error: Please select either 1 or 2\n")

print(">>> Magnitude option selection")
print(f"    MAG_COLUMN: {MAG_COLUMN}")
print(f"    MAG_ERROR_COLUMN: {MAG_ERROR_COLUMN}\n")

# Column definitions
TIME_COLUMN = 'julian_date'
STAR_NAME_COLUMN = 'Simbad_Name'
RA_COLUMN = 'ALPHA_J2000'
DEC_COLUMN = 'DELTA_J2000'
FITS_IDENTIFIER_COLUMN = 'QUADRANT_FILE'

# Plotting parameters
SG_WINDOW_LENGTH_PLOT = 15
SG_POLYNOMIAL_ORDER_PLOT = 2
Y_TICK_INCREMENT = 0.5

# Load CSV data
full_star_data = pd.read_csv(input_data_path)
full_star_data['star_index'] = full_star_data['star_index'].astype(int)

# Load comparison star mapping
target_comparison_stars = []
if os.path.exists(comparison_stars_json_path):
    with open(comparison_stars_json_path, 'r') as f:
        raw_comparison_map_indices = json.load(f)
    converted_comparison_map = {int(k): list(map(int, v)) for k, v in raw_comparison_map_indices.items()}
    target_comparison_stars = converted_comparison_map.get(star_index_to_plot_arg, [])

# Helper: clean light curve data
def get_star_data_for_plotting(star_idx, df, time_col, mag_col, err_col):
    star_data = df[df['star_index'] == star_idx].sort_values(by=time_col)
    valid = star_data[mag_col].notna() & star_data[time_col].notna() & star_data[err_col].notna()
    return (star_data[time_col][valid], star_data[mag_col][valid], star_data[err_col][valid], star_data.iloc[0])

# Helper: radial profile
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

# Main plotting
def plot_star_profile():
    star_info = full_star_data[full_star_data['star_index'] == star_index_to_plot_arg].iloc[0]
    simbad_name = star_info[STAR_NAME_COLUMN]
    simbad_type = star_info['Simbad_Type']
    is_variable = star_info['is_variable']
    quadrant_file = star_info[FITS_IDENTIFIER_COLUMN]
    fits_file = os.path.join(fits_folder_arg, quadrant_file.replace('_se.ldac', '.fit'))

    # Data
    raw_times, raw_mags, raw_errs, _ = get_star_data_for_plotting(
        star_index_to_plot_arg, full_star_data, TIME_COLUMN, MAG_COLUMN, MAG_ERROR_COLUMN
    )
    diff_times, diff_mags, diff_errs, _ = get_star_data_for_plotting(
        star_index_to_plot_arg, full_star_data, TIME_COLUMN, DIFF_MAG_COLUMN, DIFF_MAG_ERROR_COLUMN
    )

    if len(raw_mags) < SG_WINDOW_LENGTH_PLOT or len(diff_mags) < SG_WINDOW_LENGTH_PLOT:
        print(f"Not enough data to plot star {star_index_to_plot_arg}.")
        return

    fig = plt.figure(figsize=(22, 14))

    ra_str = f"{star_info[RA_COLUMN]:.6f}"
    dec_str = f"{star_info[DEC_COLUMN]:.6f}"
    plt.suptitle(
        f"Star Profile: Index {star_index_to_plot_arg} - {simbad_name} "
        f"({simbad_type}, Var: {is_variable})\nRA: {ra_str}  Dec: {dec_str}",
        fontsize=18
    )

    # --- Top-Left: Raw Magnitudes ---
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.errorbar(raw_times, raw_mags, yerr=raw_errs, fmt='o', color='blue',
                 label=f'Star {star_index_to_plot_arg} (Target)', markersize=3)
    colors = plt.cm.tab20.colors
    for i, comp_idx in enumerate(target_comparison_stars):
        comp_times, comp_mags, comp_errs, _ = get_star_data_for_plotting(
            comp_idx, full_star_data, TIME_COLUMN, MAG_COLUMN, MAG_ERROR_COLUMN
        )
        if len(comp_mags) > 0:
            ax1.errorbar(comp_times, comp_mags, yerr=comp_errs,
                         fmt='o', color=colors[i % len(colors)],
                         label=f'Star {comp_idx}', markersize=3)
    ax1.set_ylabel('Raw Magnitude')
    ax1.set_xlabel('Julian Date')
    ax1.set_title('Raw Magnitude Light Curves')
    ax1.legend(fontsize=8)
    ax1.invert_yaxis()
    ax1.grid(True)

    # --- Top-Middle: Differential LC ---
    ax2 = fig.add_subplot(2, 3, 2)
    smoothed_diff = savgol_filter(diff_mags, SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
    ax2.errorbar(diff_times, diff_mags, yerr=diff_errs, fmt='o', color='black', label='Diff Mag', markersize=3)
    ax2.plot(diff_times, smoothed_diff, '-', color='green', label='Smoothed')
    ax2.set_title('Differential Light Curve')
    ax2.set_xlabel('Julian Date')
    ax2.set_ylabel('Differential Magnitude')
    ax2.legend()
    ax2.invert_yaxis()
    ax2.grid(True)

    # --- Top-Right: 2D VIGNET (pixels only, no grid) ---
    ax3 = fig.add_subplot(2, 3, 3)
    vignet = None
    if os.path.exists(fits_file):
        with fits.open(fits_file) as hdul:
            data = hdul[0].data
            wcs = WCS(hdul[0].header)
            coord = SkyCoord(star_info[RA_COLUMN]*u.deg, star_info[DEC_COLUMN]*u.deg)
            x, y = wcs.world_to_pixel(coord)
            cut = Cutout2D(data, position=(x, y), size=(20, 20), wcs=wcs)
            vignet = cut.data
            im = ax3.imshow(vignet, cmap='viridis', origin='lower')
            ax3.set_title('2D VIGNET (First Image)')
            ax3.set_xlabel("X (pixels)")
            ax3.set_ylabel("Y (pixels)")
            ax3.grid(False)  # no grid, just ticks

    # --- Bottom-Left: FITS WCS Cutout cropped (WCS axes, no grid) ---
    ax4 = fig.add_subplot(2, 3, 4, projection=wcs)
    if os.path.exists(fits_file):
        with fits.open(fits_file) as hdul:
            data = hdul[0].data
            ax4.imshow(data, origin='lower', cmap='gray',
                       vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
            coords_list = []
            for idx in [star_index_to_plot_arg] + target_comparison_stars:
                row = full_star_data[full_star_data['star_index'] == idx].iloc[0]
                coord = SkyCoord(row[RA_COLUMN]*u.deg, row[DEC_COLUMN]*u.deg)
                coords_list.append(coord)
                px, py = wcs.world_to_pixel(coord)
                color = 'red' if idx == star_index_to_plot_arg else 'cyan'
                ax4.add_patch(Circle((px, py), 10, edgecolor=color, facecolor='none', linewidth=1.5))

            # Crop dynamically around all stars
            ra_vals = [c.ra.deg for c in coords_list]
            dec_vals = [c.dec.deg for c in coords_list]
            pixel_coords = np.array(wcs.world_to_pixel(SkyCoord(ra_vals*u.deg, dec_vals*u.deg)))
            x_min, x_max = pixel_coords[0].min() - 40, pixel_coords[0].max() + 40
            y_min, y_max = pixel_coords[1].min() - 40, pixel_coords[1].max() + 40
            ax4.set_xlim(x_min, x_max)
            ax4.set_ylim(y_min, y_max)
            ax4.coords[0].set_axislabel("RA (J2000)")
            ax4.coords[1].set_axislabel("Dec (J2000)")
            ax4.invert_yaxis()
            ax4.grid(False)  # no grid
            ax4.set_title('FITS Cutout (Target + Comparisons)')

    # --- Bottom-Middle: RA/Dec Scatter with star labels ---
    ax5 = fig.add_subplot(2, 3, 5)
    tr = star_info
    ax5.plot(tr[RA_COLUMN], tr[DEC_COLUMN], 'ro', label='Target')
    ax5.text(tr[RA_COLUMN]+0.0005, tr[DEC_COLUMN]+0.0005, f'{star_index_to_plot_arg}', color='red', fontsize=10)
    for i in target_comparison_stars:
        row = full_star_data[full_star_data['star_index'] == i].iloc[0]
        ax5.plot(row[RA_COLUMN], row[DEC_COLUMN], 'bx')
        ax5.text(row[RA_COLUMN]+0.0005, row[DEC_COLUMN]+0.0005, f'{i}', color='blue', fontsize=10)
    ax5.set_title('Target and Comparison Stars')
    ax5.set_xlabel('RA (deg)')
    ax5.set_ylabel('Dec (deg)')
    ax5.invert_xaxis()
    ax5.grid(True)
    ax5.legend()

    # --- Bottom-Right: Radial Profile ---
    ax6 = fig.add_subplot(2, 3, 6)
    if vignet is not None:
        distances, intensities = compute_radial_profile(vignet)
        ax6.scatter(distances, intensities, s=10, color='blue')
        if len(distances) > 4:
            coeffs = np.polyfit(distances, intensities, deg=4)
            poly = np.poly1d(coeffs)
            r_fit = np.linspace(0, np.max(distances), 100)
            ax6.plot(r_fit, poly(r_fit), color='black')
        ax6.set_title('Radial Profile from VIGNET')
        ax6.set_xlabel('Radius (pixels)')
        ax6.set_ylabel('Counts (ADU)')
        ax6.grid(True)

    output_file = os.path.join(
        output_plots_dir_arg,
        f"star_profile_star_{star_index_to_plot_arg}_{target_id_arg}_{MAG_COLUMN}_{w}.png"
    )
    os.makedirs(output_plots_dir_arg, exist_ok=True)
    # plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file, dpi=300)
    print(f"Saved star profile to {output_file}")
    plt.close(fig)

if __name__ == "__main__":
    plot_star_profile()






# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.ticker import MultipleLocator
# from scipy.signal import savgol_filter
# from astropy.io import fits
# from astropy.wcs import WCS
# from astropy.coordinates import SkyCoord
# import astropy.units as u
# from matplotlib.patches import Circle
# import seaborn as sns
# import json
# import sys
# import os

# # Set matplotlib style
# from matplotlib import rcParams
# rcParams["figure.figsize"] = (8, 8)
# rcParams["font.size"] = 10
# sns.set(style="whitegrid")

# # Parse command line arguments
# input_data_folder_arg = sys.argv[1]
# output_plots_dir_arg = sys.argv[2]
# target_id_arg = sys.argv[3]
# date_string_arg = sys.argv[4]
# star_index_to_plot_arg = int(sys.argv[5])
# fits_folder_arg = sys.argv[6]
# w = sys.argv[9]
# # File paths
# input_data_path = os.path.join(input_data_folder_arg, f"tracked_stars_with_differential_mags_{target_id_arg}_{date_string_arg}.csv")
# comparison_stars_json_path = os.path.join(input_data_folder_arg, f"comparison_star_map_{target_id_arg}_{date_string_arg}.json")

# # Prompt for magnitude type
# while True:
#     print(">>> Use MAG_AUTO or MAG_APER?")
#     print("(1) MAG_AUTO")
#     print("(2) MAG_APER\n")
#     which_mag = int(input("Answer: "))
#     print()
#     if which_mag == 1:
#         MAG_COLUMN = 'MAG_AUTO'
#         MAG_ERROR_COLUMN = 'MAGERR_AUTO'
#         break
#     if which_mag == 2:
#         MAG_COLUMN = 'MAG_APER'
#         MAG_ERROR_COLUMN = 'MAGERR_APER'
#         break
#     else:
#         print("Error: Please select either 1 or 2\n")

# print(">>> Magnitude option selection")
# print(f"    MAG_COLUMN: {MAG_COLUMN}")
# print(f"    MAG_ERROR_COLUMN: {MAG_ERROR_COLUMN}\n")

# # Column definitions
# DIFF_MAG_COLUMN = sys.argv[7]
# DIFF_MAG_ERROR_COLUMN = sys.argv[8]
# TIME_COLUMN = 'julian_date'
# STAR_NAME_COLUMN = 'Simbad_Name'
# RA_COLUMN = 'ALPHA_J2000'
# DEC_COLUMN = 'DELTA_J2000'
# FITS_IDENTIFIER_COLUMN = 'QUADRANT_FILE'

# # Plotting parameters
# SG_WINDOW_LENGTH_PLOT = 15
# SG_POLYNOMIAL_ORDER_PLOT = 2
# Y_TICK_INCREMENT = 0.5  # Increased from 0.2 for wider tick spacing

# # Load CSV data
# full_star_data = pd.read_csv(input_data_path)
# full_star_data['star_index'] = full_star_data['star_index'].astype(int)

# # Load and process comparison star mapping
# target_comparison_stars = []
# if os.path.exists(comparison_stars_json_path):
#     with open(comparison_stars_json_path, 'r') as f:
#         raw_comparison_map_indices = json.load(f)
#     converted_comparison_map = {int(k): list(map(int, v)) for k, v in raw_comparison_map_indices.items()}
#     target_comparison_stars = converted_comparison_map.get(star_index_to_plot_arg, [])

# # Helper function to extract clean light curve data
# def get_star_data_for_plotting(star_idx, df, time_col, mag_col, err_col):
#     star_data = df[df['star_index'] == star_idx].sort_values(by=time_col)
#     valid = star_data[mag_col].notna() & star_data[time_col].notna() & star_data[err_col].notna()
#     return (star_data[time_col][valid], star_data[mag_col][valid], star_data[err_col][valid], star_data.iloc[0])

# # Begin plotting
# def plot_star_profile():
#     star_info = full_star_data[full_star_data['star_index'] == star_index_to_plot_arg].iloc[0]
#     simbad_name = star_info[STAR_NAME_COLUMN]
#     simbad_type = star_info['Simbad_Type']
#     is_variable = star_info['is_variable']
#     quadrant_file = star_info[FITS_IDENTIFIER_COLUMN]
#     fits_file = os.path.join(fits_folder_arg, quadrant_file.replace('_se.ldac', '.fit'))

#     # Get data for target star
#     raw_times, raw_mags, raw_errs, _ = get_star_data_for_plotting(star_index_to_plot_arg, full_star_data, TIME_COLUMN, MAG_COLUMN, MAG_ERROR_COLUMN)
#     diff_times, diff_mags, diff_errs, _ = get_star_data_for_plotting(star_index_to_plot_arg, full_star_data, TIME_COLUMN, DIFF_MAG_COLUMN, DIFF_MAG_ERROR_COLUMN)

#     if len(raw_mags) < SG_WINDOW_LENGTH_PLOT or len(diff_mags) < SG_WINDOW_LENGTH_PLOT:
#         print(f"Not enough data to plot star {star_index_to_plot_arg}.")
#         return

#     fig, axs = plt.subplots(2, 2, figsize=(18, 15))
#     # plt.suptitle(f"Star Profile: Index {star_index_to_plot_arg} - {simbad_name} ({simbad_type}, Var: {is_variable})", fontsize=18)
#     ra_str = f"{star_info[RA_COLUMN]:.6f}"
#     dec_str = f"{star_info[DEC_COLUMN]:.6f}"
#     plt.suptitle(f"Star Profile: Index {star_index_to_plot_arg} - {simbad_name} ({simbad_type}, Var: {is_variable})\nRA: {ra_str}  Dec: {dec_str}", fontsize=18)


#     # --- Top-Left: Raw Magnitude Light Curve for Target and Comparison Stars ---
#     ax1 = axs[0, 0]
    
#     # Collect all magnitudes for y-axis scaling
#     all_mags = [raw_mags]
    
#     # Plot target star (data points and error bars only)
#     ax1.errorbar(raw_times, raw_mags, yerr=raw_errs, fmt='o', color='blue', alpha=0.6, label=f'Star {star_index_to_plot_arg} (Target)', markersize=3)

#     # Plot comparison stars
#     colors = ['green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan', 'indigo', 'maroon', 'black', 'teal', 'limegreen', 'blueviolet', 'crimson']  # Cycle through colors
#     for i, comp_star_idx in enumerate(target_comparison_stars):
#         comp_times, comp_mags, comp_errs, _ = get_star_data_for_plotting(comp_star_idx, full_star_data, TIME_COLUMN, MAG_COLUMN, MAG_ERROR_COLUMN)
#         if len(comp_mags) > 0:  # Plot even if < SG_WINDOW_LENGTH_PLOT since no smoothing
#             color = colors[i % len(colors)]
#             ax1.errorbar(comp_times, comp_mags, yerr=comp_errs, fmt='o', color=color, alpha=0.6, label=f'Star {comp_star_idx}', markersize=3)
#             all_mags.append(comp_mags)

#     # Set y-axis for combined raw magnitudes
#     if all_mags:
#         combined_mags = pd.concat([pd.Series(m) for m in all_mags if not m.empty])
#         if not combined_mags.empty:
#             raw_median = combined_mags.median()
#             raw_padding = max(3 * combined_mags.std(), Y_TICK_INCREMENT * 2) if combined_mags.std() > 0 else Y_TICK_INCREMENT * 2
#             raw_y_min = np.floor((raw_median - raw_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
#             raw_y_max = np.ceil((raw_median + raw_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
#             ax1.set_ylim(raw_y_max, raw_y_min)  # Inverted for magnitudes
#             ax1.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))
#         else:
#             ax1.set_ylim(0.5, -0.5)
#             ax1.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))
    
#     ax1.set_ylabel('Raw Magnitude', fontsize=12)
#     ax1.set_xlabel('Julian Date', fontsize=12)
#     ax1.set_title('Raw Magnitude Light Curves', fontsize=14)
#     ax1.legend(loc='upper right', fontsize=10)

#     # --- Top-Right: Differential Light Curve with Error Bars ---
#     ax3 = axs[0, 1]
#     smoothed_diff = savgol_filter(diff_mags, SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
#     ax3.errorbar(diff_times, diff_mags, yerr=diff_errs, fmt='o', color='black', alpha=0.6, label='Diff Mag', markersize=3)
#     ax3.plot(diff_times, smoothed_diff, '-', color='green', label='Smoothed Diff')
#     ax3.set_title('Differential Light Curve', fontsize=14)
#     ax3.set_xlabel('Julian Date', fontsize=12)
#     ax3.set_ylabel('Differential Magnitude', fontsize=12)
#     ax3.legend(loc='upper right')
#     ax3.invert_yaxis()

#     # --- Bottom-Left: FITS Image with Target and Comparisons ---
#     if os.path.exists(fits_file):
#         with fits.open(fits_file) as hdul:
#             data = hdul[0].data
#             wcs = WCS(hdul[0].header)
#             ax4 = plt.subplot(2, 2, 3, projection=wcs)
            

#             ax4.imshow(data, origin='lower', cmap='gray', 
#                        vmin=np.percentile(data, 5), vmax=np.percentile(data, 99))
#             ax4.set_title('FITS Image with Target/Comparisons', fontsize=14)

#             # WCS ticks and labels only — no default 0.0 to 1.0 ticks
#             ax4.coords.grid(True, color='white', ls='--', alpha=0.5)
#             ax4.coords[0].set_axislabel('RA', fontsize=12)
#             ax4.coords[1].set_axislabel('Dec', fontsize=12)
#             ax4.coords[0].set_ticklabel(size=11)
#             ax4.coords[1].set_ticklabel(size=11)
#             ax4.set_facecolor('black')

#             # Turn off the native (0.0–1.0) ticks completely
#             ax4.get_xaxis().set_visible(False)
#             ax4.get_yaxis().set_visible(False)

#             # ax4.set_xticks([])
#             # ax4.set_yticks([])

#             # Plot target and comparison stars
#             coords_list = []
#             for idx in [star_index_to_plot_arg] + target_comparison_stars:
#                 row = full_star_data[full_star_data['star_index'] == idx].iloc[0]
#                 coord = SkyCoord(row[RA_COLUMN]*u.deg, row[DEC_COLUMN]*u.deg)
#                 coords_list.append(coord)
#                 px, py = wcs.world_to_pixel(coord)
#                 color = 'red' if idx == star_index_to_plot_arg else 'cyan'
#                 ax4.add_patch(Circle((px, py), 10, edgecolor=color, facecolor='none', linewidth=1.5))

#             # Dynamic zoom to fit all stars
#             ra_vals = [c.ra.deg for c in coords_list]
#             dec_vals = [c.dec.deg for c in coords_list]
#             pixel_coords = np.array(wcs.world_to_pixel(SkyCoord(ra_vals*u.deg, dec_vals*u.deg)))
#             x_min, x_max = pixel_coords[0].min() - 40, pixel_coords[0].max() + 40
#             y_min, y_max = pixel_coords[1].min() - 40, pixel_coords[1].max() + 40
#             ax4.set_xlim(x_min, x_max)
#             ax4.set_ylim(y_min, y_max)
#             ax4.invert_yaxis()

#     # --- Bottom-Right: RA/Dec Location Map ---
#     ax5 = axs[1, 1]
#     tr = full_star_data[full_star_data['star_index'] == star_index_to_plot_arg].iloc[0]
#     ax5.plot(tr[RA_COLUMN], tr[DEC_COLUMN], 'ro', label='Target')
#     ax5.text(tr[RA_COLUMN], tr[DEC_COLUMN]+0.0005, f'{star_index_to_plot_arg}', color='red', fontsize=10)
#     for i in target_comparison_stars:
#         row = full_star_data[full_star_data['star_index'] == i].iloc[0]
#         ax5.plot(row[RA_COLUMN], row[DEC_COLUMN], 'bx')
#         ax5.text(row[RA_COLUMN], row[DEC_COLUMN]+0.0005, f'{i}', color='blue', fontsize=11)
#     ax5.set_title(f'Target and Comparison Stars', fontsize=14)
#     ax5.set_xlabel('RA (deg)', fontsize=12)
#     ax5.set_ylabel('Dec (deg)', fontsize=12)
#     ax5.invert_xaxis()  # Flip RA to increase leftward
#     ax5.grid(True)
#     ax5.legend()

#     output_file = os.path.join(output_plots_dir_arg, f"star_profile_star_{star_index_to_plot_arg}_{target_id_arg}_{MAG_COLUMN}_{w}.png")
#     os.makedirs(output_plots_dir_arg, exist_ok=True)
#     plt.savefig(output_file, dpi=300)
#     print(f"Saved star profile to {output_file}")
#     plt.close(fig)

# if __name__ == "__main__":
#     plot_star_profile()


