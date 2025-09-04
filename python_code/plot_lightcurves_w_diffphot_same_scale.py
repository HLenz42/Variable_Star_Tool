import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import os
import sys

target = sys.argv[3]
date = sys.argv[4]
w = sys.argv[7]

# --- Configuration ---
# Input path for your data with differential magnitudes
input_data_path = os.path.join(sys.argv[1], f"tracked_stars_with_differential_mags_{target}_{date}.csv")

# Output directory for saving plots
output_plots_dir = sys.argv[2]

# Number of plots per grid (5 columns * 10 rows)
plots_per_grid = 50
num_cols = 5
num_rows = 10

# Define the differential magnitude, time, and error column names
DIFF_MAG_COLUMN = sys.argv[5]
DIFF_MAG_ERROR_COLUMN = sys.argv[6]
TIME_COLUMN = 'julian_date'
STAR_NAME_COLUMN = 'Simbad_Name'

# Define Simbad variable star types for grouping context
VARIABLE_SIMBAD_TYPES_FOR_GROUPING = [
    'EB*', 'RR*', 'BY*', 'V*', 'LP?', 'LP*', 'Ro*', 'Er*', 'RS*', 'cC*'
]

# Y-axis tick mark increment (same for all subplots)
Y_TICK_INCREMENT = 0.2  # Magnitude units between major tick marks

# --- End Configuration ---

try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        print("\n" * 50)

def plot_differential_light_curves(star_indices, data_df, plot_title_prefix, output_filename):
    """
    Plots differential light curves for a given list of star indices in a grid,
    including a shaded error region, with consistent y-axis tick mark spacing.

    Args:
        star_indices (list): A list of star_index values to plot.
        data_df (pd.DataFrame): The DataFrame containing all star data.
        plot_title_prefix (str): Prefix for the plot title (e.g., "Variable Stars").
        output_filename (str): The full path including filename for saving the plot.
    """
    if not star_indices:
        print(f"No stars to plot for '{plot_title_prefix}'. Skipping plot generation.")
        return

    # Ensure output directory exists
    os.makedirs(output_plots_dir, exist_ok=True)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols * 3, num_rows * 2.5),
                            sharex=False, sharey=False)  # Individual scaling for each subplot

    axs = axs.flatten()  # Flatten for easy iteration

    plotted_count = 0
    for i, star_idx in enumerate(star_indices[:plots_per_grid]):
        if plotted_count >= plots_per_grid:
            break

        ax = axs[plotted_count]

        # Get data for the current star, sorted by time
        star_data = data_df[data_df['star_index'] == star_idx].sort_values(by=TIME_COLUMN)

        # Skip if insufficient data
        if star_data[DIFF_MAG_COLUMN].isnull().all() or len(star_data) < 2:
            ax.set_visible(False)
            plotted_count += 1
            continue

        # Extract data for plotting
        times = star_data[TIME_COLUMN]
        differential_magnitudes = star_data[DIFF_MAG_COLUMN]
        differential_magnitude_errors = star_data[DIFF_MAG_ERROR_COLUMN]

        # Plot differential magnitude vs time
        ax.plot(times, differential_magnitudes, 'o-', markersize=2, linewidth=0.5, label='Differential Magnitude')

        # Add shaded error region
        ax.fill_between(times,
                        differential_magnitudes - differential_magnitude_errors,
                        differential_magnitudes + differential_magnitude_errors,
                        color='lightgray', alpha=0.5, label='Error')

        # Get metadata for title
        simbad_type = star_data['Simbad_Type'].iloc[0]
        simbad_name = star_data[STAR_NAME_COLUMN].iloc[0]
        is_variable_flag = star_data['is_variable'].iloc[0]

        # Set title
        ax.set_title(f"Star {star_idx}\n({simbad_type}, {simbad_name}, Var: {is_variable_flag})", fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.set_xlabel('Julian Date', fontsize=7)
        ax.set_ylabel('Diff. Magnitude', fontsize=7)

        # Set y-axis limits and tick marks
        mag_range = differential_magnitudes.max() - differential_magnitudes.min()
        if not np.isnan(mag_range) and mag_range > 0:
            median_val = differential_magnitudes.median()
            # Use 3 * std_dev or a minimum range to ensure visibility
            padding = max(3 * differential_magnitudes.std(), Y_TICK_INCREMENT * 2)  # At least 2 tick increments
            y_min = median_val - padding
            y_max = median_val + padding
            # Round limits to nearest tick increment for clean tick marks
            y_min = np.floor(y_min / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            y_max = np.ceil(y_max / Y_TICK_INCREMENT) * Y_TICK_INCREMENT
            ax.set_ylim(y_max, y_min)  # Inverted for magnitudes (larger values lower)
        else:
            # Default range if no variation
            ax.set_ylim(0.5, -0.5)

        # Set consistent y-axis tick mark spacing
        ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_INCREMENT))

        plotted_count += 1

    # Hide unused subplots
    for j in range(plotted_count, len(axs)):
        axs[j].set_visible(False)

    fig.suptitle(plot_title_prefix, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to: {output_filename}")

if __name__ == "__main__":
    print(f"Loading data from: {input_data_path}")
    try:
        full_star_data = pd.read_csv(input_data_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_data_path}")
        sys.exit(1)

    # Prepare data
    full_star_data['star_index'] = full_star_data['star_index'].astype(int)

    # Group stars based on is_variable
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
    print(f"Found {len(variable_stars_flagged)} stars flagged as VARIABLE.")
    print(f"Found {len(unclassified_and_not_variable)} stars UNCLASSIFIED by Simbad AND NOT flagged as variable.")
    print(f"Found {len(other_classified_and_not_variable)} stars OTHER CLASSIFIED by Simbad AND NOT flagged as variable.")
    print("----------------------------------------------------")

    # Generate plots
    output_dir = os.path.join(output_plots_dir, "light_curves_same_scale")
    os.makedirs(output_dir, exist_ok=True)

    plot_differential_light_curves(
        variable_stars_flagged,
        full_star_data,
        "VARIABLE Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_variable_diff_light_curves_{w}.png")
    )

    plot_differential_light_curves(
        unclassified_and_not_variable,
        full_star_data,
        "UNCLASSIFIED (Non-Variable Flagged) Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_unclassified_diff_light_curves_{w}.png")
    )

    plot_differential_light_curves(
        other_classified_and_not_variable,
        full_star_data,
        "Other Classified (Non-Variable Flagged) Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_classified_diff_light_curves_{w}.png")
    )

    print("\nScript finished generating differential light curve plots.")