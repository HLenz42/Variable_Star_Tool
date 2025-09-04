import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import os
import sys

target = sys.argv[3]
date = sys.argv[4]
w = sys.argv[7]

# --- Configuration ---
# Input path for your data with differential magnitudes
# This should be the output from your differential photometry script
input_data_path = os.path.join(sys.argv[1], f"tracked_stars_with_differential_mags_{target}_{date}.csv")

# Output directory for saving plots
output_plots_dir = sys.argv[2]



# Number of plots per grid (5 columns * 10 rows)
plots_per_grid = 50
num_cols = 5
num_rows = 10

# Define the differential magnitude, time, and error column names
# These are the NEW columns from your differential photometry script
DIFF_MAG_COLUMN = sys.argv[5]
DIFF_MAG_ERROR_COLUMN = sys.argv[6]
TIME_COLUMN = 'julian_date'
STAR_NAME_COLUMN = 'Simbad_Name' # To be used in plot titles

# Define Simbad variable star types for the 'is_variable' flag logic
# This list should match what you used to set the 'is_variable' column in the previous script.
# We'll use the 'is_variable' column primarily, but this is kept for context.
VARIABLE_SIMBAD_TYPES_FOR_GROUPING = [
    'EB*', 'RR*', 'BY*', 'V*', 'LP?', 'LP*', 'Ro*', 'Er*', 'RS*', 'cC*'
]
# --- End Configuration ---

try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        # Define a dummy clear_output if not in IPython
        # This will just print newlines to push content down
        print("\n" * 50)

def plot_differential_light_curves(star_indices, data_df, plot_title_prefix, output_filename):
    """
    Plots differential light curves for a given list of star indices in a grid,
    including a shaded error region.

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
                            sharex=False, sharey=False) # Keep sharex/sharey as False for individual scaling

    axs = axs.flatten() # Flatten the 2D array of axes for easy iteration

    plotted_count = 0
    # Iterate only through the first 'plots_per_grid' stars
    for i, star_idx in enumerate(star_indices[:plots_per_grid]):
        if plotted_count >= plots_per_grid:
            break

        ax = axs[plotted_count]

        # Get all observations for the current star_index, sorted by time
        star_data = data_df[data_df['star_index'] == star_idx].sort_values(by=TIME_COLUMN)

        # Skip if there's not enough data to plot the light curve
        # (e.g., if differential magnitude was NaN for all points)
        if star_data[DIFF_MAG_COLUMN].isnull().all() or len(star_data) < 2:
            ax.set_visible(False) # Hide empty subplot
            plotted_count += 1
            continue

        # Extract data for plotting differential magnitudes and errors
        times = star_data[TIME_COLUMN]
        differential_magnitudes = star_data[DIFF_MAG_COLUMN]
        differential_magnitude_errors = star_data[DIFF_MAG_ERROR_COLUMN]

        # Plot differential magnitude vs time
        ax.plot(times, differential_magnitudes, 'o-', markersize=2, linewidth=0.5, label='Differential Magnitude')

        # Add shaded error region for differential magnitude
        ax.fill_between(times,
                        differential_magnitudes - differential_magnitude_errors,
                        differential_magnitudes + differential_magnitude_errors,
                        color='lightgray', alpha=0.5, label='Error')

        # Get Simbad_Type and Simbad_Name for the title
        simbad_type = star_data['Simbad_Type'].iloc[0]
        simbad_name = star_data[STAR_NAME_COLUMN].iloc[0] # Using STAR_NAME_COLUMN for consistency
        is_variable_flag = star_data['is_variable'].iloc[0] # Get the 'is_variable' flag

        # Set title with star index, Simbad type, Simbad name, and variability flag
        ax.set_title(f"Star {star_idx}\n({simbad_type}, {simbad_name}, Var: {is_variable_flag})", fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.set_xlabel('Julian Date', fontsize=7) # Add X-axis label
        ax.set_ylabel('Diff. Magnitude', fontsize=7) # Add Y-axis label

        # Adjust y-axis limits for differential magnitudes
        # Differential magnitudes should ideally hover around 0.
        # You might want to set symmetric limits around 0 for better comparison.
        # Calculate a reasonable range based on the data's spread
        mag_range = differential_magnitudes.max() - differential_magnitudes.min()
        if not np.isnan(mag_range) and mag_range > 0:
            # Set limits to be symmetric around 0 or the median, with some padding
            median_val = differential_magnitudes.median()
            # Use a fixed padding or a multiple of the standard deviation
            padding = 3 * differential_magnitudes.std() if differential_magnitudes.std() > 0 else 0.1
            ax.set_ylim(median_val + padding, median_val - padding) # Inverted for magnitudes
        else:
            # Fallback if no variation or all NaNs
            ax.set_ylim(0.5, -0.5) # A default small range around zero

        plotted_count += 1

    # Hide any unused subplots
    for j in range(plotted_count, len(axs)):
        axs[j].set_visible(False)

    fig.suptitle(plot_title_prefix, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust layout to make space for suptitle
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to: {output_filename}")


if __name__ == "__main__":
    # if len(sys.argv) != 4:
    #     print("Usage: python plot_differential_lightcurves.py <path_to_data_folder> <output_directory_for_plots>")
    #     print("Example: python plot_differential_lightcurves.py /data/processed /plots/differential_lcs")
    #     sys.exit(1)

    # Load the data
    print(f"Loading data from: {input_data_path}")
    try:
        full_star_data = pd.read_csv(input_data_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_data_path}")
        sys.exit(1)

    # --- Prepare data for plotting ---
    full_star_data['star_index'] = full_star_data['star_index'].astype(int)

    # Group stars based on the 'is_variable' column
    variable_stars_flagged = full_star_data[full_star_data['is_variable'] == 'yes']['star_index'].unique().tolist()
    
    # Group for "Unclassified" based on Simbad_Type AND 'is_variable' == 'No'
    unclassified_and_not_variable = full_star_data[
        (full_star_data['Simbad_Type'] == 'UNCLASSIFIED') &
        (full_star_data['is_variable'] == 'no')
    ]['star_index'].unique().tolist()

    # Group for "Other Classified (non-variable)"
    # These are stars that SIMBAD classified as something other than UNCLASSIFIED,
    # AND were not flagged as 'is_variable' == 'Yes' by your script.
    other_classified_and_not_variable = full_star_data[
        (full_star_data['Simbad_Type'] != 'UNCLASSIFIED') &
        (full_star_data['is_variable'] == 'no')
    ]['star_index'].unique().tolist()

    print("----------------------------------------------------")
    print(f"Found {len(variable_stars_flagged)} stars flagged as VARIABLE.")
    print(f"Found {len(unclassified_and_not_variable)} stars UNCLASSIFIED by Simbad AND NOT flagged as variable.")
    print(f"Found {len(other_classified_and_not_variable)} stars OTHER CLASSIFIED by Simbad AND NOT flagged as variable.")
    print("----------------------------------------------------")
    
    # --- Generate Plots ---

    output_dir = os.path.join(output_plots_dir, "light_curves_not_scaled")  #Output for this target

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Plot stars flagged as VARIABLE
    plot_differential_light_curves(
        variable_stars_flagged,
        full_star_data,
        "VARIABLE Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_variable_diff_light_curves_{w}.png")
    )

    # Plot stars that are UNCLASSIFIED by Simbad AND not flagged as variable
    plot_differential_light_curves(
        unclassified_and_not_variable,
        full_star_data,
        "UNCLASSIFIED (Non-Variable Flagged) Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_unclassified_diff_light_curves_{w}.png")
    )

    # Plot stars that are OTHER CLASSIFIED by Simbad AND not flagged as variable
    plot_differential_light_curves(
        other_classified_and_not_variable,
        full_star_data,
        "Other Classified (Non-Variable Flagged) Stars Differential Light Curves",
        os.path.join(output_dir, f"{target}_classified_diff_light_curves_{w}.png")
    )

    print("\nScript finished generating differential light curve plots.")