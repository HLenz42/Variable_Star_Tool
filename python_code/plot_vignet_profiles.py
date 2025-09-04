import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def parse_vignet(vignet_str, vignet_shape=(20, 20)):
    """Parse VIGNET string into a 2D array."""
    try:
        rows = vignet_str.strip('()').split(');(')
        vignet = np.array([row.split(';') for row in rows], dtype=float)
        if vignet.shape != vignet_shape:
            raise ValueError(f"VIGNET shape {vignet.shape} does not match expected {vignet_shape}")
        return vignet
    except Exception as e:
        print(f"Error parsing VIGNET: {e}")
        return None

def compute_radial_profile(vignet, center=None):
    """Compute radial profile with distances from center to each pixel."""
    if vignet is None:
        return None, None
    height, width = vignet.shape
    if center is None:
        center = ((width - 1) / 2, (height - 1) / 2)  # Center for 20x20: (9.5, 9.5)
    y, x = np.indices((height, width))
    distances = np.sqrt((x - center[0])**2 + (y - center[1])**2).flatten()
    intensities = vignet.flatten()
    return distances, intensities

def plot_vignet_profile(csv_file, start_row, num_rows=10, vignet_shape=(20, 20)):
    """Generate multi-panel plots for sources in CSV with pixel-based radial profiles and polynomial fit."""
    # Read CSV
    try:
        df = pd.read_csv(csv_file, sep=',')
    except Exception as e:
        print(f"Error reading CSV {csv_file}: {e}")
        return
    
    # Validate rows
    total_rows = len(df)
    print(f"CSV {csv_file} contains {total_rows} sources")
    if start_row < 0 or start_row >= total_rows:
        print(f"Error: Start row {start_row} is out of range (0 to {total_rows-1})")
        return
    end_row = min(start_row + num_rows, total_rows)
    print(f"Processing sources from row {start_row} to {end_row-1}")

    # Log file
    log_file = os.path.join(os.path.dirname(csv_file), f"vignet_plot_log_{os.path.basename(csv_file).replace('.csv', '.txt')}")
    log_lines = [f"===== VIGNET Plot Log for {os.path.basename(csv_file)} ====="]

    for idx in range(start_row, end_row):
        row = df.iloc[idx]
        source_id = row['INDEX']
        ra = row['ALPHA_J2000']
        dec = row['DELTA_J2000']
        a_image = row['A_IMAGE'] * 6       # Profile RMS along major axis [pixel]
        b_image = row['B_IMAGE'] * 6      # Profile RMS along minor axis [pixel]
        background = row['BACKGROUND'] # Background at centroid position [count]
        threshold = row['THRESHOLD']   # Detection threshold above background [count]
        mag_auto = row['MAG_AUTO']
        mag_aper = row.get('MAG_APER', np.nan)  # Handle multiple MAG_APER if needed
        vignet_str = row['VIGNET']
        quadrant_file = row['QUADRANT_FILE']

        # Parse VIGNET
        vignet = parse_vignet(vignet_str, vignet_shape)
        if vignet is None:
            log_lines.append(f"Source {source_id} (row {idx}): Failed to parse VIGNET")
            print(f"Skipping source {source_id} (row {idx}): Invalid VIGNET")
            continue

        # Compute radial profile
        distances, intensities = compute_radial_profile(vignet)
        if distances is None:
            log_lines.append(f"Source {source_id} (row {idx}): Failed to compute radial profile")
            print(f"Skipping source {source_id} (row {idx}): No radial profile")
            continue

        # Polynomial fit (2nd degree)
        valid_mask = ~np.isnan(distances) & ~np.isnan(intensities)
        if valid_mask.sum() > 3:  # Need at least 4 points for a 2nd-degree fit
            poly_coeffs = np.polyfit(distances[valid_mask], intensities[valid_mask], deg=4)
            poly_func = np.poly1d(poly_coeffs)
            r_fine = np.linspace(0, np.max(distances[valid_mask]), 100)
            poly_intensities = poly_func(r_fine)
        else:
            r_fine, poly_intensities = None, None
            log_lines.append(f"Source {source_id} (row {idx}): Insufficient valid data for polynomial fit")

        # Create three-panel plot
        fig = plt.figure(figsize=(18, 5))
        fig.suptitle(f"Source {source_id}: RA={ra:.6f}, Dec={dec:.6f}, MAG_AUTO={mag_auto:.2f}, MAG_APER={mag_aper:.2f}", fontsize=14)

        # Left: 2D image
        ax1 = fig.add_subplot(131)
        vmin, vmax = np.nanmin(vignet), np.nanmax(vignet)
        im = ax1.imshow(vignet, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
        ax1.set_xlabel('X Pixel')
        ax1.set_ylabel('Y Pixel')
        ax1.set_title('VIGNET Cutout (2D)')
        plt.colorbar(im, ax=ax1, label='Counts (ADU)')

        # Middle: 3D surface
        ax2 = fig.add_subplot(132, projection='3d')
        x, y = np.arange(vignet_shape[1]), np.arange(vignet_shape[0])
        X, Y = np.meshgrid(x, y)
        surf = ax2.plot_surface(X, Y, vignet, cmap='viridis', vmin=vmin, vmax=vmax)
        ax2.set_xlabel('X Pixel')
        ax2.set_ylabel('Y Pixel')
        ax2.set_zlabel('Counts (ADU)')
        ax2.set_title('VIGNET Surface (3D)')
        plt.colorbar(surf, ax=ax2, label='Counts (ADU)')

        # Right: Radial profile
        ax3 = fig.add_subplot(133)
        avg_kron_radii = (a_image + b_image) / 2
        thresh = background + threshold
        ax3.scatter(distances, intensities, color='blue', s=10, label='Pixel Intensities')
        if r_fine is not None:
            ax3.plot(r_fine, poly_intensities, color='black', linestyle='-', label='Polynomial Fit')
        ax3.axvline(x=3, color='red', linestyle='--', label='Fixed Aperture')
        ax3.axvline(x=avg_kron_radii, color='green', linestyle='--', label='Auto Aperture AVG')
        ax3.axvline(x=a_image, color='orange', linestyle='--', label='Major Axis')
        ax3.axvline(x=b_image, color='purple', linestyle='--', label='Minor Axis')
        # ax3.axhline(y=background, color='dimgrey', linestyle='--', label=f'Background: {background:.0f}')
        # ax3.axhline(y=thresh, color='cyan', linestyle='--', label=f'Threshold: {thresh:.0f}')
        ax3.set_xlabel('Radius (pixels)')
        ax3.set_ylabel('Counts (ADU)')
        ax3.set_title('Radial Profile')
        ax3.grid(True, linestyle='--', alpha=0.5)
        ax3.legend()

        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust for suptitle
        os.makedirs(os.path.join(os.path.dirname(csv_file), f"vignet_profiles_{star_filtered}"), exist_ok=True)
        output_plot = os.path.join(os.path.dirname(csv_file), f"vignet_profiles_{star_filtered}", f"vignet_profile_source_{source_id}_row_{idx}.png")
        plt.savefig(output_plot, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot for source {source_id}: {output_plot}")
        log_lines.append(f"Source {source_id} (row {idx}): Saved plot to {output_plot}")

    # Save log
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Saved log: {log_file}")

if __name__ == "__main__":
    output_folder = sys.argv[1]
    star_filtered = sys.argv[2]
    csv_files = sorted(glob.glob(os.path.join(output_folder, "*.csv")))
    if not csv_files:
        print(f"Error: No CSV files found in {output_folder}")
        sys.exit(1)

    csv_file = csv_files[0]
    print(f"Using first CSV: {csv_file}")

    # Get number of rows
    df = pd.read_csv(csv_file, sep=',')
    total_rows = len(df)
    print(f"CSV contains {total_rows} sources")

    # Prompt user for start row
    try:
        start_row = int(input(f"Enter start row (0 to {total_rows-1}): "))
    except ValueError:
        print("Error: Invalid input, must be an integer")
        sys.exit(1)

    plot_vignet_profile(csv_file, start_row, num_rows=10, vignet_shape=(20, 20))

# import os
# import glob
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import sys

# def parse_vignet(vignet_str, vignet_shape=(20,20)):
#     """Parse VIGNET string into a 2D array."""
#     try:
#         rows = vignet_str.strip('()').split(');(')
#         vignet = np.array([row.split(';') for row in rows], dtype=float)
#         if vignet.shape != vignet_shape:
#             raise ValueError(f"VIGNET shape {vignet.shape} does not match expected {vignet_shape}")
#         return vignet
#     except Exception as e:
#         print(f"Error parsing VIGNET: {e}")
#         return None

# def compute_radial_profile(vignet, center=None):
#     """Compute radial profile from VIGNET array."""
#     if vignet is None:
#         return None, None
#     height, width = vignet.shape
#     if center is None:
#         center = ((width - 1) / 2, (height - 1) / 2)  # Center for 4x4: (1.5, 1.5)
#     y, x = np.indices((height, width))
#     distances = np.sqrt((x - center[0])**2 + (y - center[1])**2).flatten()
#     intensities = vignet.flatten()
    
#     bins = np.arange(0, np.ceil(np.max(distances)) + 0.25, 0.25)  # 0 to max radius, 0.25-pixel steps
#     radial_profile = []
#     for r_min, r_max in zip(bins[:-1], bins[1:]):
#         mask = (distances >= r_min) & (distances < r_max)
#         if mask.any():
#             radial_profile.append(np.mean(intensities[mask]))
#         else:
#             radial_profile.append(np.nan)
#     radii = (bins[:-1] + bins[1:]) / 2
#     return radii, radial_profile

# def plot_vignet_profile(csv_file, start_row, num_rows=10, vignet_shape=(4, 4)):
#     """Generate multi-panel plots for sources in CSV."""
#     # Read CSV
#     try:
#         df = pd.read_csv(csv_file, sep=',')
#     except Exception as e:
#         print(f"Error reading CSV {csv_file}: {e}")
#         return
    
#     # Validate rows
#     total_rows = len(df)
#     print(f"CSV {csv_file} contains {total_rows} sources")
#     if start_row < 0 or start_row >= total_rows:
#         print(f"Error: Start row {start_row} is out of range (0 to {total_rows-1})")
#         return
#     end_row = min(start_row + num_rows, total_rows)
#     print(f"Processing sources from row {start_row} to {end_row-1}")

#     # Log file
#     log_file = os.path.join(os.path.dirname(csv_file), f"vignet_plot_log_{os.path.basename(csv_file).replace('.csv', '.txt')}")
#     log_lines = [f"===== VIGNET Plot Log for {os.path.basename(csv_file)} ====="]

#     for idx in range(start_row, end_row):
#         row = df.iloc[idx]
#         source_id = row['INDEX']
#         ra = row['ALPHA_J2000']
#         dec = row['DELTA_J2000']
#         a_image = row['A_IMAGE']       # Profile RMS along major axis                              [pixel]
#         b_image = row['B_IMAGE']       # Profile RMS along minor axis                              [pixel]
#         background = row['BACKGROUND'] # Background at centroid position                           [count]
#         threshold = row['THRESHOLD']   # Detection threshold above background                      [count]
#         mag_auto = row['MAG_AUTO']
#         mag_aper = row.get('MAG_APER', np.nan)  # Handle multiple MAG_APER if needed
#         vignet_str = row['VIGNET']
#         quadrant_file = row['QUADRANT_FILE']

#         # Parse VIGNET
#         vignet = parse_vignet(vignet_str, vignet_shape)
#         if vignet is None:
#             log_lines.append(f"Source {source_id} (row {idx}): Failed to parse VIGNET")
#             print(f"Skipping source {source_id} (row {idx}): Invalid VIGNET")
#             continue

#         # Compute radial profile
#         radii, radial_profile = compute_radial_profile(vignet)
#         if radii is None:
#             log_lines.append(f"Source {source_id} (row {idx}): Failed to compute radial profile")
#             print(f"Skipping source {source_id} (row {idx}): No radial profile")
#             continue

#         # Create three-panel plot
#         fig = plt.figure(figsize=(18, 5))
#         fig.suptitle(f"Source {source_id}: RA={ra:.6f}, Dec={dec:.6f}, MAG_AUTO={mag_auto:.2f}, MAG_APER={mag_aper:.2f}", fontsize=14)

#         # Left: 2D image
#         ax1 = fig.add_subplot(131)
#         vmin, vmax = np.nanmin(vignet), np.nanmax(vignet)
#         im = ax1.imshow(vignet, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
#         ax1.set_xlabel('X Pixel')
#         ax1.set_ylabel('Y Pixel')
#         ax1.set_title('VIGNET Cutout (2D)')
#         plt.colorbar(im, ax=ax1, label='Counts (ADU)')

#         # Middle: 3D surface
#         ax2 = fig.add_subplot(132, projection='3d')
#         x, y = np.arange(vignet_shape[1]), np.arange(vignet_shape[0])
#         X, Y = np.meshgrid(x, y)
#         surf = ax2.plot_surface(X, Y, vignet, cmap='viridis', vmin=vmin, vmax=vmax)
#         ax2.set_xlabel('X Pixel')
#         ax2.set_ylabel('Y Pixel')
#         ax2.set_zlabel('Counts (ADU)')
#         ax2.set_title('VIGNET Surface (3D)')
#         plt.colorbar(surf, ax=ax2, label='Counts (ADU)')

#         # Right: Radial profile
#         ax3 = fig.add_subplot(133)
#         avg_kron_radii = (a_image + b_image) / 2
#         thresh = background + threshold
#         ax3.scatter(radii, radial_profile, color='blue', s=20, label='Radial Profile')
#         ax3.axvline(x=3, color='red', linestyle='--', label='Fixed Aperture')
#         ax3.axvline(x=avg_kron_radii, color='green', linestyle='--', label='Auto Aperture AVG')
#         ax3.axvline(x=a_image, color='orange', linestyle='--', label='Major Axis')
#         ax3.axvline(x=b_image, color='purple', linestyle='--', label='Minor Axis')
#         ax3.axhline(y=background, color='dimgrey', linestyle='--', label=f'Background: {background:.0f}')
#         ax3.axhline(y=thresh, color='cyan', linestyle='--', label=f'Threshold: {thresh:.0f}')
#         ax3.set_xlabel('Radius (pixels)')
#         ax3.set_ylabel('Mean Counts')
#         ax3.set_title('Radial Profile')
#         ax3.grid(True, linestyle='--', alpha=0.5)
#         ax3.legend()

#         plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust for suptitle
#         os.makedirs(os.path.join(os.path.dirname(csv_file), f"vignet_profiles_{star_filtered}"), exist_ok=True)
#         output_plot = os.path.join(os.path.dirname(csv_file), f"vignet_profiles_{star_filtered}", f"vignet_profile_source_{source_id}_row_{idx}.png")
#         plt.savefig(output_plot, dpi=300, bbox_inches='tight')
#         plt.close()
#         print(f"Saved plot for source {source_id}: {output_plot}")
#         log_lines.append(f"Source {source_id} (row {idx}): Saved plot to {output_plot}")

#     # Save log
#     with open(log_file, 'w') as f:
#         f.write('\n'.join(log_lines))
#     print(f"Saved log: {log_file}")

# if __name__ == "__main__":
#     # if len(sys.argv) < 2:
#     #     print("Usage: python plot_vignet_profiles.py [output_folder]")
#     #     sys.exit(1)

#     output_folder = sys.argv[1]
#     star_filtered = sys.argv[2]
#     csv_files = sorted(glob.glob(os.path.join(output_folder, "*.csv")))
#     if not csv_files:
#         print(f"Error: No CSV files found in {output_folder}")
#         sys.exit(1)

#     csv_file = csv_files[0]
#     print(f"Using first CSV: {csv_file}")

#     # Get number of rows
#     df = pd.read_csv(csv_file, sep=',')
#     total_rows = len(df)
#     print(f"CSV contains {total_rows} sources")

#     # Prompt user for start row
#     try:
#         start_row = int(input(f"Enter start row (0 to {total_rows-1}): "))
#     except ValueError:
#         print("Error: Invalid input, must be an integer")
#         sys.exit(1)

#     plot_vignet_profile(csv_file, start_row, num_rows=10, vignet_shape=(20, 20))