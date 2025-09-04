import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

# Command-line arguments
input_dir = sys.argv[1]
target = sys.argv[2]
date = sys.argv[3]
star_index = sys.argv[4]
version = sys.argv[5]
w = sys.argv[8]
if version == "1":
    data_csv = os.path.join(input_dir, f"star_data_{star_index}_{target}_{date}.csv")
if version == "2":
    data_csv = os.path.join(input_dir, f"star_data_filtered_{star_index}_{target}_{date}.csv")
output_file = os.path.join(input_dir, "light_curve_plot", f"diff_curve_{star_index}_{target}_{date}_{w}_V{version}.png")

# Load data
try:
    data = pd.read_csv(data_csv, sep=',', engine='python')
except FileNotFoundError:
    print(f"Error: CSV file not found at {data_csv}")
    sys.exit(1)

# Column names (matching extract_star_data.py and pipeline)
image_index = data["image_index"]
julian_date = data["julian_date"]
simbad_name = data["Simbad_Name"].iloc[0]  # Single value for title
is_variable = data["is_variable"].iloc[0]  # Single value for title
diff_mag = data[sys.argv[6]]
diff_mag_err = data[sys.argv[7]]

# Create plot
fig, ax = plt.subplots(figsize=(11, 8))

# Plot differential magnitude with error bars
ax.errorbar(julian_date, diff_mag, yerr=diff_mag_err, fmt='bo', alpha=0.8, markersize=5, label='Diff Mag', capsize=3)

# # Add image_index labels next to each point
# for idx, row in data.iterrows():
#     ax.text(row["julian_date"], row["differential_magnitude"] + 0.06, f'{int(row["image_index"])}', 
#             color='blue', fontsize=6, ha='center', va='bottom')

# Set axis limits with padding
if not diff_mag.empty:
    jd_min, jd_max = julian_date.min(), julian_date.max()
    jd_padding = (jd_max - jd_min) * 0.05  # 5% padding
    mag_median = diff_mag.median()
    mag_min = diff_mag.min()
    mag_max = diff_mag.max()
    mag_padding = 0.04
    # mag_padding = max(0.5 * diff_mag.std(), 0.3) if diff_mag.std() > 0 else 0.1
    ax.set_xlim(jd_min - jd_padding, jd_max + jd_padding)
    ax.set_ylim(mag_max + mag_padding, mag_min - mag_padding)  # Inverted for magnitudes
else:
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, -0.5)

# Customize axes
ax.set_xlabel("Julian Date", fontsize=14)
ax.set_ylabel("Differential Magnitude", fontsize=14)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(loc='upper right')

# Title with star info
plt.title(f"Differential Magnitude Light Curve for {target}\nStar Index = {star_index}, Simbad Name = {simbad_name}, Variable = {is_variable}", 
          fontsize=14, pad=10)

# Save plot
os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Saved differential light curve to {output_file}")
plt.close(fig)