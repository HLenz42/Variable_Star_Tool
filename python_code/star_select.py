import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


def filter_stars_from_catalogue(input_csv_folder, output_csv_folder, mag_faint, mag_bright, star_cat_plots):
    os.makedirs(output_csv_folder, exist_ok=True)
    csv_files = [f for f in os.listdir(input_csv_folder) if f.endswith(".csv")]

    for csv_file in csv_files:
        input_path = os.path.join(input_csv_folder, csv_file)
        df = pd.read_csv(input_path)

        # Check for required columns
        required_columns = {"FLAGS", "MAG_AUTO", "FLUX_RADIUS", "ELONGATION"}
        if not required_columns.issubset(df.columns):
            print(f"[✗] Skipping {csv_file}: missing one or more required columns: {required_columns}")
            continue

        # Basic filters
        filtered = df[df["FLAGS"] <= flags]  # Higher the number the more sceptical sources allowed though
        filtered = filtered[(filtered["MAG_AUTO"] >= mag_bright) & (filtered["MAG_AUTO"] <= mag_faint)] #Only sources in a certain range of magnitudeds allowed to be stars

        if filtered.empty:
            print(f"[!] No stars passed mag filter in {csv_file}")
            continue

        median_flux_radius = np.median(filtered["FLUX_RADIUS"])
        filtered = filtered[
            (filtered["FLUX_RADIUS"] >= min_flux_radiuse * median_flux_radius) &
            (filtered["FLUX_RADIUS"] <= max_flux_radiuse * median_flux_radius)
        ]
        filtered = filtered[filtered["ELONGATION"] >= min_elongation]

        # Plotting
        x_limits = (-1,15) # FLUX_RADIUSE
        y_limits = (2.2,14) # MAG_AUTO
        plt.figure()
        plt.plot(df["FLUX_RADIUS"], df["MAG_AUTO"], "ko", ms=2, label="All")
        plt.plot(filtered["FLUX_RADIUS"], filtered["MAG_AUTO"], "ro", ms=2, label="Selected Stars")
        plt.axhline(mag_faint, color="blue", linestyle=":", label="Mag limits")
        plt.axhline(mag_bright, color="blue", linestyle=":")
        plt.xlabel("FLUX_RADIUS")
        plt.ylabel("MAG_AUTO")
        plt.xlim(x_limits)
        plt.ylim(y_limits)
        plt.gca().invert_yaxis()
        plt.legend()
        plt.title(f"Star Selection: {csv_file}")
        plot_path = os.path.join(star_cat_plots, csv_file.replace(".csv", "_stars_plot.png"))
        plt.savefig(plot_path, bbox_inches="tight")
        plt.close()

        # Save filtered catalog
        output_path = os.path.join(output_csv_folder, csv_file.replace(".csv", "_stars.csv"))
        filtered.to_csv(output_path, index=False)
        print(f">>> Saved star catalog: {output_path}")

all_sources_files = sys.argv[1]
output_stars_files = sys.argv[2]
mag_faint = float(sys.argv[3])
mag_bright = float(sys.argv[4])
star_cat_plots = sys.argv[2] + "/star_cat_plots"
flags = float(sys.argv[5])
min_flux_radiuse = float(sys.argv[6])
max_flux_radiuse = float(sys.argv[7])
min_elongation = float(sys.argv[8])


print(f"\n======== filtering values ========")
print(f"mag_faint        = {mag_faint}")
print(f"mag_bright       = {mag_bright}")
print(f"max_flags        = {flags}")
print(f"min_flux_radiuse = {min_flux_radiuse}")
print(f"max_flux_radiuse = {max_flux_radiuse}")
print(f"min_elongation   = {min_elongation}")
print(f"==================================")
print()


os.makedirs(star_cat_plots, exist_ok=True)

filter_stars_from_catalogue(all_sources_files, output_stars_files, mag_faint, mag_bright, star_cat_plots)

