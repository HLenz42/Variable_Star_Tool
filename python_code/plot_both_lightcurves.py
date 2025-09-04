import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.signal import savgol_filter
import os
import sys

# --- Configuration ---
input_data_folder = sys.argv[1]
target_id_for_filename = sys.argv[3]
date_string_for_filename = sys.argv[4]
input_data_path = os.path.join(input_data_folder, f"tracked_stars_with_differential_mags_{target_id_for_filename}_{date_string_for_filename}.csv")
w = sys.argv[7]

# Output directory for saving all plots
output_plots_dir = sys.argv[2]

# General Plotting Grid Configuration
GRID_PLOTS_PER_PAGE = 50
GRID_NUM_COLS = 5
GRID_NUM_ROWS = 10

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
        print("Error: Please select either 1, 2" + "\n")

print(">>> Magnitude option selection")
print(f"    MAG_COLUMN: {MAG_COLUMN}")
print(f"    MAG_ERROR_COLUMN: {MAG_ERROR_COLUMN}\n")


# Define column names from your CSV
DIFF_MAG_COLUMN = sys.argv[5]
DIFF_MAG_ERROR_COLUMN = sys.argv[6]
TIME_COLUMN = 'julian_date'
STAR_NAME_COLUMN = 'Simbad_Name'

# Y-axis tick mark increment (same for raw and differential magnitudes)
Y_TICK_INCREMENT = 0.2  # Magnitude units between major tick marks

# Savitzky-Golay filter parameters for smoothing
SG_WINDOW_LENGTH_PLOT = 15
SG_POLYNOMIAL_ORDER_PLOT = 2

# Define Simbad variable star types for grouping context
VARIABLE_SIMBAD_TYPES_FOR_GROUPING = [
    'EB*', 'RR*', 'BY*', 'V*', 'LP?', 'LP*', 'Ro*', 'Er*', 'RS*', 'cC*'
]

print(f"\n--- Plotting Configuration ---")
print(f"Grid layout: {GRID_NUM_ROWS} rows x {GRID_NUM_COLS} cols ({GRID_PLOTS_PER_PAGE} plots/page)")
print(f"Y-axis tick increment: {Y_TICK_INCREMENT} magnitudes")
print(f"Smoothing filter window length: {SG_WINDOW_LENGTH_PLOT}")
print(f"Smoothing filter polynomial order: {SG_POLYNOMIAL_ORDER_PLOT}")
print(f"----------------------------\n")

# --- Utility to clear output ---
try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        print("\n" * 50)

def plot_single_axis_grid(star_indices, data_df, plot_title_prefix, output_filename, plot_raw_mag=False):
    """
    Plots light curves (either raw or differential) for a given list of star indices in a grid,
    including a shaded error region and a smoothed line.

    Args:
        star_indices (list): A list of star_index values to plot.
        data_df (pd.DataFrame): The DataFrame containing all star data.
        plot_title_prefix (str): Prefix for the plot title.
        output_filename (str): The full path for saving the plot.
        plot_raw_mag (bool): If True, plots MAG_AUTO. If False, plots differential_magnitude.
    """
    if not star_indices:
        print(f"No stars to plot for '{plot_title_prefix}'. Skipping plot generation.")
        return

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    fig, axs = plt.subplots(GRID_NUM_ROWS, GRID_NUM_COLS, figsize=(GRID_NUM_COLS * 3, GRID_NUM_ROWS * 2.5),
                            sharex=False, sharey=False)

    axs = axs.flatten()

    plotted_count = 0
    for i, star_idx in enumerate(star_indices[:GRID_PLOTS_PER_PAGE]):
        if plotted_count >= GRID_PLOTS_PER_PAGE:
            break

        ax = axs[plotted_count]

        current_mag_col = MAG_COLUMN if plot_raw_mag else DIFF_MAG_COLUMN
        current_err_col = MAG_ERROR_COLUMN if plot_raw_mag else DIFF_MAG_ERROR_COLUMN
        y_label = 'Magnitude' if plot_raw_mag else 'Diff. Magnitude'
        plot_label_prefix = 'Raw' if plot_raw_mag else 'Diff.'

        star_data = data_df[data_df['star_index'] == star_idx].sort_values(by=TIME_COLUMN)

        times = star_data[TIME_COLUMN]
        magnitudes = star_data[current_mag_col]
        magnitude_errors = star_data[current_err_col]

        valid_indices = magnitudes.notna() & times.notna() & magnitude_errors.notna()
        valid_times = times[valid_indices]
        valid_magnitudes = magnitudes[valid_indices]
        valid_errors = magnitude_errors[valid_indices]

        if valid_magnitudes.empty or len(valid_magnitudes) < SG_WINDOW_LENGTH_PLOT:
            ax.set_visible(False)
            plotted_count += 1
            continue

        ax.plot(valid_times, valid_magnitudes, 'o', markersize=1, label=f'{plot_label_prefix} Data', alpha=0.7)
        ax.fill_between(valid_times, valid_magnitudes - valid_errors, valid_magnitudes + valid_errors,
                        color='lightgray', alpha=0.5, label='Error')

        smoothed = savgol_filter(valid_magnitudes.to_numpy(), SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
        ax.plot(valid_times, smoothed, color='black', alpha=0.9, linewidth=1, label='Smoothed')

        simbad_type = star_data['Simbad_Type'].iloc[0]
        simbad_name = star_data[STAR_NAME_COLUMN].iloc[0]
        is_variable_flag = star_data['is_variable'].iloc[0]

        ax.set_title(f"Star {star_idx}\n({simbad_type}, {simbad_name}, Var: {is_variable_flag})", fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.set_xlabel('Julian Date', fontsize=7)
        ax.set_ylabel(y_label, fontsize=7)

        # Set y-axis limits and tick marks
        if not valid_magnitudes.empty:
            median_val = valid_magnitudes.median()
            padding = max(3 * valid_magnitudes.std(), Y_TICK_INCREMENT * 2) if valid_magnitudes.std() > 0 else Y_TICK_INCREMENT * 2
            y_min = np.floor((median_val - padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            y_max = np.ceil((median_val + padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            ax.set_ylim(y_max, y_min)  # Inverted for magnitudes
            ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))
        else:
            ax.set_ylim(0.5, -0.5)
            ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))

        plotted_count += 1

    for j in range(plotted_count, len(axs)):
        axs[j].set_visible(False)

    fig.suptitle(plot_title_prefix, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to: {output_filename}")

def plot_dual_axis_grid(star_indices, data_df, plot_title_prefix, output_filename):
    """
    Plots raw and differential light curves for a list of stars in a grid,
    each subplot having dual Y-axes with consistent tick mark spacing.

    Args:
        star_indices (list): A list of star_index values to plot.
        data_df (pd.DataFrame): The DataFrame containing all star data.
        plot_title_prefix (str): Prefix for the plot title.
        output_filename (str): The full path for saving the plot.
    """
    if not star_indices:
        print(f"No stars to plot for '{plot_title_prefix}'. Skipping plot generation.")
        return

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    fig, axs = plt.subplots(GRID_NUM_ROWS, GRID_NUM_COLS, figsize=(GRID_NUM_COLS * 4, GRID_NUM_ROWS * 3),
                            sharex=False)

    axs = axs.flatten()

    plotted_count = 0
    for i, star_idx in enumerate(star_indices[:GRID_PLOTS_PER_PAGE]):
        if plotted_count >= GRID_PLOTS_PER_PAGE:
            break

        ax1 = axs[plotted_count]  # Primary axis for raw magnitudes
        ax2 = ax1.twinx()         # Secondary axis for differential magnitudes

        star_data = data_df[data_df['star_index'] == star_idx].sort_values(by=TIME_COLUMN)

        # Raw Magnitudes
        raw_magnitudes = star_data[MAG_COLUMN]
        raw_magnitudes_errors = star_data[MAG_ERROR_COLUMN]
        valid_raw_indices = raw_magnitudes.notna() & raw_magnitudes_errors.notna() & star_data[TIME_COLUMN].notna()
        valid_raw_times = star_data[TIME_COLUMN][valid_raw_indices]
        valid_raw_magnitudes = raw_magnitudes[valid_raw_indices]
        valid_raw_errors = raw_magnitudes_errors[valid_raw_indices]

        # Differential Magnitudes
        differential_magnitudes = star_data[DIFF_MAG_COLUMN]
        differential_magnitude_errors = star_data[DIFF_MAG_ERROR_COLUMN]
        valid_diff_indices = differential_magnitudes.notna() & differential_magnitude_errors.notna() & star_data[TIME_COLUMN].notna()
        valid_diff_times = star_data[TIME_COLUMN][valid_diff_indices]
        valid_diff_magnitudes = differential_magnitudes[valid_diff_indices]
        valid_diff_errors = differential_magnitude_errors[valid_diff_indices]

        if (valid_raw_magnitudes.empty or len(valid_raw_magnitudes) < SG_WINDOW_LENGTH_PLOT) or \
           (valid_diff_magnitudes.empty or len(valid_diff_magnitudes) < SG_WINDOW_LENGTH_PLOT):
            ax1.set_visible(False)
            plotted_count += 1
            continue

        # Plot Raw Magnitudes on ax1
        ax1.plot(valid_raw_times, valid_raw_magnitudes, 'o', markersize=1, color='blue', label='Raw Mag Data', alpha=0.7)
        ax1.fill_between(valid_raw_times, valid_raw_magnitudes - valid_raw_errors, valid_raw_magnitudes + valid_raw_errors,
                         color='lightblue', alpha=0.3, label='Raw Mag Error')
        smoothed_raw = savgol_filter(valid_raw_magnitudes.to_numpy(), SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
        ax1.plot(valid_raw_times, smoothed_raw, color='darkblue', linewidth=1, label='Smoothed Raw Mag')

        ax1.set_xlabel('Julian Date', fontsize=7)
        ax1.set_ylabel('Raw Magnitude', color='blue', fontsize=7)
        ax1.tick_params(axis='y', labelcolor='blue', labelsize=7)
        ax1.tick_params(axis='x', labelsize=7)

        # Set y-axis for raw magnitudes
        if not valid_raw_magnitudes.empty:
            raw_median = valid_raw_magnitudes.median()
            raw_padding = max(3 * valid_raw_magnitudes.std(), Y_TICK_INCREMENT * 2) if valid_raw_magnitudes.std() > 0 else Y_TICK_INCREMENT * 2
            raw_y_min = np.floor((raw_median - raw_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            raw_y_max = np.ceil((raw_median + raw_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            ax1.set_ylim(raw_y_max, raw_y_min)  # Inverted for magnitudes
            ax1.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))
        else:
            ax1.set_ylim(0.5, -0.5)
            ax1.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))

        # Plot Differential Magnitudes on ax2
        ax2.plot(valid_diff_times, valid_diff_magnitudes, 's', markersize=1, color='red', label='Diff. Mag Data', alpha=0.7)
        ax2.fill_between(valid_diff_times, valid_diff_magnitudes - valid_diff_errors, valid_diff_magnitudes + valid_diff_errors,
                         color='lightcoral', alpha=0.3, label='Diff. Mag Error')
        smoothed_diff = savgol_filter(valid_diff_magnitudes.to_numpy(), SG_WINDOW_LENGTH_PLOT, SG_POLYNOMIAL_ORDER_PLOT)
        ax2.plot(valid_diff_times, smoothed_diff, color='darkred', linewidth=1, label='Smoothed Diff. Mag')

        ax2.set_ylabel('Diff. Magnitude', color='red', fontsize=7)
        ax2.tick_params(axis='y', labelcolor='red', labelsize=7)

        # Set y-axis for differential magnitudes
        if not valid_diff_magnitudes.empty:
            diff_median = valid_diff_magnitudes.median()
            diff_padding = max(3 * valid_diff_magnitudes.std(), Y_TICK_INCREMENT * 2) if valid_diff_magnitudes.std() > 0 else Y_TICK_INCREMENT * 2
            diff_y_min = np.floor((diff_median - diff_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            diff_y_max = np.ceil((diff_median + diff_padding) / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            ax2.set_ylim(diff_y_max, diff_y_min)  # Inverted for magnitudes
            ax2.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))
        else:
            ax2.set_ylim(0.5, -0.5)
            ax2.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))

        # Get metadata for title
        star_info = star_data.drop_duplicates(subset=['star_index']).iloc[0]
        simbad_type = star_info['Simbad_Type']
        simbad_name = star_info[STAR_NAME_COLUMN]
        is_variable_flag = star_info['is_variable']

        ax1.set_title(f"Star {star_idx}\n({simbad_type}, {simbad_name}, Var: {is_variable_flag})", fontsize=8)

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=6)

        plotted_count += 1

    for j in range(plotted_count, len(axs)):
        axs[j].set_visible(False)

    fig.suptitle(plot_title_prefix, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Dual-axis grid plot saved to: {output_filename}")

if __name__ == "__main__":
    # Load the data
    print(f"Loading data from: {input_data_path}")
    try:
        full_star_data = pd.read_csv(input_data_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_data_path}")
        sys.exit(1)

    full_star_data['star_index'] = full_star_data['star_index'].astype(int)
    full_star_data = full_star_data.sort_values(by=[TIME_COLUMN, 'star_index']).reset_index(drop=True)

    # Group stars for plotting
    variable_stars_flagged = full_star_data[full_star_data['is_variable'] == 'yes']['star_index'].unique().tolist()
    unclassified_and_not_variable = full_star_data[
        (full_star_data['Simbad_Type'] == 'UNCLASSIFIED') &
        (full_star_data['is_variable'] == 'no')
    ]['star_index'].unique().tolist()
    other_classified_and_not_variable = full_star_data[
        (full_star_data['Simbad_Type'] != 'UNCLASSIFIED') &
        (full_star_data['is_variable'] == 'no')
    ]['star_index'].unique().tolist()

    print("----------------------------------------------------")
    print(f"Summary of Star Groups for Plotting:")
    print(f"  - {len(variable_stars_flagged)} stars flagged as VARIABLE.")
    print(f"  - {len(unclassified_and_not_variable)} stars UNCLASSIFIED by Simbad AND NOT flagged as variable.")
    print(f"  - {len(other_classified_and_not_variable)} stars OTHER CLASSIFIED by Simbad AND NOT flagged as variable.")
    print("----------------------------------------------------")

    # Generate plots
    output_subdir = os.path.join(output_plots_dir, "dual_axis_photometry_plots")
    os.makedirs(output_subdir, exist_ok=True)
    print(f"\nSaving plots to: {output_subdir}")

    # Single-Axis Differential Magnitude Plots
    print("\n--- Generating Single-Axis Differential Magnitude Plots ---")
    plot_single_axis_grid(
        variable_stars_flagged,
        full_star_data,
        "VARIABLE Stars Differential Light Curves (with Smoothing)",
        os.path.join(output_subdir, f"{target_id_for_filename}_variable_diff_light_curves_smoothed_{MAG_COLUMN}_{w}.png"),
        plot_raw_mag=False
    )

    plot_single_axis_grid(
        unclassified_and_not_variable,
        full_star_data,
        "UNCLASSIFIED (Non-Variable Flagged) Stars Differential Light Curves (with Smoothing)",
        os.path.join(output_subdir, f"{target_id_for_filename}_unclassified_non_variable_diff_light_curves_smoothed_{MAG_COLUMN}_{w}.png"),
        plot_raw_mag=False
    )

    plot_single_axis_grid(
        other_classified_and_not_variable,
        full_star_data,
        "Other Classified (Non-Variable Flagged) Stars Differential Light Curves (with Smoothing)",
        os.path.join(output_subdir, f"{target_id_for_filename}_classified_non_variable_diff_light_curves_smoothed_{MAG_COLUMN}_{w}.png"),
        plot_raw_mag=False
    )

    # Dual-Axis Raw vs. Differential Magnitude Plots
    print("\n--- Generating Dual-Axis Raw vs. Differential Magnitude Plots ---")
    plot_dual_axis_grid(
        variable_stars_flagged,
        full_star_data,
        "VARIABLE Stars: Raw vs. Differential Light Curves",
        os.path.join(output_subdir, f"{target_id_for_filename}_variable_raw_vs_diff_light_curves_{MAG_COLUMN}.png")
    )

    plot_dual_axis_grid(
        unclassified_and_not_variable,
        full_star_data,
        "UNCLASSIFIED (Non-Variable Flagged) Stars: Raw vs. Differential Light Curves",
        os.path.join(output_subdir, f"{target_id_for_filename}_unclassified_raw_vs_diff_light_curves_{MAG_COLUMN}.png")
    )

    plot_dual_axis_grid(
        other_classified_and_not_variable,
        full_star_data,
        "Other Classified (Non-Variable Flagged) Stars: Raw vs. Differential Light Curves",
        os.path.join(output_subdir, f"{target_id_for_filename}_classified_raw_vs_diff_light_curves_{MAG_COLUMN}.png")
    )

    print("\nScript finished generating all plots.")