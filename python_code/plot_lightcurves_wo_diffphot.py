import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# --- Configuration ---
# Input path for your classified data
input_classified_data_path = sys.argv[1] # Expected to be the full path to tracked_stars_with_simbad.csv

# Output directory for saving plots
output_plots_dir = sys.argv[2] # Expected to be the directory to save plots

target = sys.argv[3]
date = sys.argv[4]


# Number of plots per grid (5 columns * 10 rows)
plots_per_grid = 50
num_cols = 5
num_rows = 10

# Define Simbad variable star types based on your previous output
VARIABLE_STAR_TYPES = ['EB*', 'RR*', 'BY*', 'V*', 'LP?', 'LP*', 'Ro*', 'Er*', 'RS*', 'cC*']

# Define the magnitude, time, and error column names

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

TIME_COLUMN = 'julian_date'
# --- End Configuration ---

try:
    from IPython.display import clear_output
except ImportError:
    def clear_output(wait=False):
        pass

def plot_light_curves(star_indices, data_df, plot_title_prefix, output_filename):
    """
    Plots light curves for a given list of star indices in a grid,
    including a shaded error region.

    Args:
        star_indices (list): A list of star_index values to plot.
        data_df (pd.DataFrame): The DataFrame containing all star data.
        plot_title_prefix (str): Prefix for the plot title (e.g., "UNCLASSIFIED Stars").
        output_filename (str): The full path including filename for saving the plot.
    """
    if not star_indices:
        print(f"No stars to plot for '{plot_title_prefix}'. Skipping plot generation.")
        return

    # Ensure output directory exists
    os.makedirs(output_plots_dir, exist_ok=True)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols * 3, num_rows * 2.5),
                            sharex=False, sharey=False)

    axs = axs.flatten()

    plotted_count = 0
    for i, star_idx in enumerate(star_indices[:plots_per_grid]):
        if plotted_count >= plots_per_grid:
            break

        ax = axs[plotted_count]

        star_data = data_df[data_df['star_index'] == star_idx].sort_values(by=TIME_COLUMN)

        if len(star_data) < 2:
            ax.set_visible(False)
            plotted_count += 1
            continue

        # Extract data for plotting
        times = star_data[TIME_COLUMN]
        magnitudes = star_data[MAG_COLUMN]
        magnitude_errors = star_data[MAG_ERROR_COLUMN] # NEW: Get magnitude errors

        # Plot magnitude vs time
        ax.plot(times, magnitudes, 'o-', markersize=2, linewidth=0.5, label='Magnitude')

        # NEW: Add shaded error region
        ax.fill_between(times,
                        magnitudes - magnitude_errors,
                        magnitudes + magnitude_errors,
                        color='lightgray', alpha=0.5, label='Error') # You can choose your desired color and transparency

        # Get Simbad_Type and Simbad_Name for the title
        simbad_type = star_data['Simbad_Type'].iloc[0]
        simbad_name = star_data['Simbad_Name'].iloc[0]

        ax.set_title(f"Star {star_idx}\n({simbad_type}, {simbad_name})", fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=7)

        # Optional: Adjust y-axis limits (consider if you want this for all plots or per-plot)
        # You might need to experiment with these values or make them dynamic
        # For magnitudes, lower values mean brighter, so typically ylim is inverted:
        # ax.set_ylim(magnitudes.max() + 0.1, magnitudes.min() - 0.1)
        # Or, if you want a fixed range around the mean (like your original snippet):
        # plotScale = 0.5 # Example value for magnitude range
        # if not magnitudes.empty:
        #     mean_mag = magnitudes.mean()
        #     ax.set_ylim(mean_mag + plotScale, mean_mag - plotScale)


        plotted_count += 1

    # Hide any unused subplots
    for j in range(plotted_count, len(axs)):
        axs[j].set_visible(False)

    fig.suptitle(plot_title_prefix, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to: {output_filename}")


if __name__ == "__main__":
    # if len(sys.argv) != 4:
    #     print("Usage: python your_script_name.py <path_to_classified_csv> <output_directory_for_plots>")
    #     sys.exit(1)

    # Load the data
    print(f"Loading data from: {input_classified_data_path}")
    try:
        full_star_data = pd.read_csv(input_classified_data_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_classified_data_path}")
        sys.exit(1)

    # --- Prepare data for plotting ---
    full_star_data['star_index'] = full_star_data['star_index'].astype(int)

    unclassified_stars = full_star_data[full_star_data['Simbad_Type'] == 'UNCLASSIFIED']['star_index'].unique().tolist()
    
    variable_stars = full_star_data[full_star_data['Simbad_Type'].isin(VARIABLE_STAR_TYPES)]['star_index'].unique().tolist()

    other_classified_stars = full_star_data[
        (full_star_data['Simbad_Type'] != 'UNCLASSIFIED') &
        (~full_star_data['Simbad_Type'].isin(VARIABLE_STAR_TYPES))
    ]['star_index'].unique().tolist()

    print("----------------------------------------------------")
    print(f"Found {len(unclassified_stars)} UNCLASSIFIED stars.")
    print(f"Found {len(variable_stars)} VARIABLE stars.")
    print(f"Found {len(other_classified_stars)} OTHER CLASSIFIED stars.")
    print("----------------------------------------------------")
    

    # --- Generate Plots ---

    output_dir = os.path.join(output_plots_dir)  #Output for this target

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plot_light_curves(
        unclassified_stars,
        full_star_data,
        "UNCLASSIFIED Stars Light Curves (Magnitude vs. Julian Date with Error)",
        os.path.join(output_dir, f"{target}_unclassified_light_curves_{MAG_COLUMN}.png")
    )

    plot_light_curves(
        variable_stars,
        full_star_data,
        "Variable Stars Light Curves (Magnitude vs. Julian Date with Error)",
        os.path.join(output_dir, f"{target}_variable_light_curves_{MAG_COLUMN}.png")
    )

    plot_light_curves(
        other_classified_stars,
        full_star_data,
        "Other Classified Stars Light Curves (Magnitude vs. Julian Date with Error)",
        os.path.join(output_dir, f"{target}_classified_light_curves_{MAG_COLUMN}.png")
    )

    print("\nScript finished generating plots.")