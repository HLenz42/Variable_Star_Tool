import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import seaborn as sns
from matplotlib.gridspec import GridSpec

# Input CSV filename
tracked_csv = "tracked_stars.csv"
# tracked_csv = "tracked_stars_all_lower_sigma.csv"

# Matplotlib settings
rcParams["font.size"] = 10
sns.set(style="whitegrid")

def plot_tracking_map(stars_df, output_path, quadrant_name, target, date, total_images):
    """
    Faster scatter plot: one call for all points in the quadrant.
    Points plotted for non-NaN ALPHA_J2000/DELTA_J2000, colored by tracking_percentage.
    Y-axis ticks evenly spaced for each star_index.
    """
    # Dynamically set figure height based on number of stars
    star_indices = stars_df['star_index'].unique()
    num_stars = len(star_indices)
    fig_height = max(300, num_stars / 50)  # ~2,000 stars -> ~20 inches
    fig = plt.figure(figsize=(35, fig_height))
    gs = GridSpec(1, 2, width_ratios=[80, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    cmap = plt.cm.viridis_r

    # Map star_index to evenly spaced positions (0, 1, 2, ...)
    star_index_map = {idx: i for i, idx in enumerate(sorted(star_indices))}
    stars_df = stars_df.copy()
    stars_df['plot_y'] = stars_df['star_index'].map(star_index_map)

    # Filter non-NaN entries
    valid_points = stars_df[stars_df['ALPHA_J2000'].notna() & stars_df['DELTA_J2000'].notna()]

    # Scatter all points at once
    scatter = ax.scatter(
        valid_points['image_number'],
        valid_points['plot_y'],
        c=valid_points['tracking_percentage'],
        cmap=cmap,
        norm=plt.Normalize(80, 100),
        s=8,
        # alpha=0.5
    )

    # Set axis labels and title
    ax.set_xlabel('Image Number')
    ax.set_ylabel('Star Index')
    ax.set_title(f"Tracking Map: {quadrant_name}")

    # Set axis limits
    ax.set_xlim(0.5, total_images + 0.5)
    ax.set_ylim(-0.5, len(star_indices) - 0.5)

    # X-axis ticks for every whole number
    ax.set_xticks(range(1, total_images + 1))
    ax.set_xticklabels([str(i) for i in range(1, total_images + 1)], rotation=45)

    # Y-axis ticks for every star_index, evenly spaced
    ax.set_yticks(range(len(star_indices)))
    ax.set_yticklabels([str(int(i)) for i in sorted(star_indices)])

    # Add gridlines
    ax.grid(True, which='both', linestyle='--', alpha=0.7)

    # Color bar
    cax = fig.add_subplot(gs[1])
    cbar = fig.colorbar(scatter, cax=cax, label='Tracking Percentage (%)')
    cbar.set_ticks(np.linspace(80, 100, 5))

    # Save plot
    save_path = os.path.join(output_path, f"{quadrant_name}_tracking_map_{target}_{date}.png")
    plt.savefig(save_path, bbox_inches='tight', dpi=200)
    plt.close()

def summarize_tracking_percentages(df):
    """
    Print a summary of tracking percentages and NaN diagnostics.
    """
    tracking_percentages = df.groupby('star_index')['tracking_percentage'].first()
    print("\n===== Tracking Percentage Summary =====")
    print(f"Total stars: {len(tracking_percentages)}")
    print(f"Mean: {tracking_percentages.mean():.1f}% | Median: {tracking_percentages.median():.1f}%")
    print(f"Min: {tracking_percentages.min():.1f}% | Max: {tracking_percentages.max():.1f}%")
    bins = np.histogram(tracking_percentages, bins=[0, 20, 40, 60, 80, 100])[0]
    for i, count in enumerate(bins):
        print(f"{i*20}-{(i+1)*20}%: {count} stars")
    print("\nNaN Diagnostics:")
    nan_quadrants = df['QUADRANT_FILE'].isna().sum()
    print(f"Rows with NaN QUADRANT_FILE: {nan_quadrants} ({100 * nan_quadrants / len(df):.1f}%)")
    print(f"Total unique julian_date values: {len(df['julian_date'].unique())}")
    print("Top 5 julian_date counts:\n", df['julian_date'].value_counts().head())
    print("=======================================\n")

def main(csv_path, output_folder, target, date):
    """
    Generate tracking map scatter plots for 16 quadrants.
    """
    os.makedirs(output_folder, exist_ok=True)
    df = pd.read_csv(os.path.join(csv_path, tracked_csv))

    # Debug: Initial inspection
    print("\n===== DEBUG: INITIAL DATA INSPECTION =====")
    print(f"Total rows: {len(df)}")
    print("Columns:", df.columns.tolist())
    print(f"Unique julian_date count: {df['julian_date'].nunique()}")
    print("===========================================\n")

    # Assign image_number based on julian_date
    unique_dates = sorted(df['julian_date'].unique())
    image_numbers_map = {jd: i + 1 for i, jd in enumerate(unique_dates)}
    df['image_number'] = df['julian_date'].map(image_numbers_map)
    total_images = len(unique_dates)
    print(f"[Info] Total images: {total_images}")

    # Extract quadrant number
    df['quadrant_num'] = df['QUADRANT_FILE'].str.extract(r'_q(\d{1,2})').astype(float)

    # Compute tracking percentage per star
    non_nan_counts = df.groupby('star_index').apply(
        lambda x: x[['ALPHA_J2000', 'DELTA_J2000']].notna().all(axis=1).sum(),
        include_groups=False
    )
    df = df.merge(non_nan_counts.rename('tracking_count').to_frame(),
                  left_on='star_index', right_index=True)
    df['tracking_percentage'] = (df['tracking_count'] / total_images) * 100

    summarize_tracking_percentages(df)

    # Process each quadrant (0–15)
    for q_num in range(16):
        stars_in_quadrant = df[df['quadrant_num'] == q_num]

        if stars_in_quadrant.empty:
            print(f"[!] No stars in quadrant {q_num}. Skipping.")
            continue

        print(f">>> Quadrant {q_num}: {stars_in_quadrant['star_index'].nunique()} stars, "
              f"{stars_in_quadrant['image_number'].nunique()} images")

        plot_tracking_map(stars_in_quadrant, output_folder, f"quadrant_{q_num}", target, date, total_images)

    print(f"[✓] All tracking map plots generated for 16 quadrants.")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python plot_tracking_map.py <csv_dir> <output_dir> <target> <date>")
        sys.exit(1)

    csv_dir = sys.argv[1]
    output_dir = sys.argv[2]
    target = sys.argv[3]
    date = sys.argv[4]

    main(csv_path=csv_dir, output_folder=output_dir, target=target, date=date)



# import os
# import sys
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib import rcParams
# import seaborn as sns
# from matplotlib.gridspec import GridSpec

# # Input CSV filename
# tracked_csv = "tracked_stars_all_lower_sigma.csv"

# # Matplotlib settings
# rcParams["figure.figsize"] = (20, 12)
# rcParams["font.size"] = 10
# sns.set(style="whitegrid")


# def plot_tracking_map(stars_df, output_path, quadrant_name, target, date, total_images):
#     """
#     Faster scatter plot: one call for all points in the quadrant.
#     """
#     fig = plt.figure(figsize=(20, 12))
#     gs = GridSpec(1, 2, width_ratios=[10, 1], wspace=0.05)
#     ax = fig.add_subplot(gs[0])
#     cmap = plt.cm.hsv_r

#     # Scatter all points at once
#     scatter = ax.scatter(
#         stars_df['image_number'],
#         stars_df['star_index'],
#         c=stars_df['tracking_percentage'],
#         cmap=cmap,
#         norm=plt.Normalize(80, 100),  # Normalize between 80 and 100%
#         s=2,
#         alpha=0.5
#     )

#     ax.set_xlabel('Image Number')
#     ax.set_ylabel('Star Index')
#     ax.set_title(f"Tracking Map: {quadrant_name}")
#     ax.set_xlim(0.5, total_images + 0.5)

#     if total_images <= 50:
#         ax.set_xticks(range(1, total_images + 1))
#         ax.set_xticklabels([str(i) for i in range(1, total_images + 1)], rotation=45)

#     # Subsample y-axis ticks for clarity
#     star_indices = stars_df['star_index'].unique()
#     if len(star_indices) > 100:
#         tick_indices = star_indices[::50]
#     else:
#         tick_indices = star_indices
#     ax.set_yticks(tick_indices)
#     ax.set_yticklabels([str(int(i)) for i in tick_indices])

#     ax.grid(True, which='both', linestyle='--', alpha=0.7)

#     # Color bar
#     cax = fig.add_subplot(gs[1])
#     cbar = fig.colorbar(scatter, cax=cax, label='Tracking Percentage (%)')
#     cbar.set_ticks(np.linspace(80, 100, 5))

#     save_path = os.path.join(output_path, f"{quadrant_name}_tracking_map_{target}_{date}.png")
#     plt.savefig(save_path, bbox_inches='tight', dpi=200)
#     plt.close()


# def summarize_tracking_percentages(df):
#     tracking_percentages = df.groupby('star_index')['tracking_percentage'].first()
#     print("\n===== Tracking Percentage Summary =====")
#     print(f"Total stars: {len(tracking_percentages)}")
#     print(f"Mean: {tracking_percentages.mean():.1f}% | Median: {tracking_percentages.median():.1f}%")
#     print(f"Min: {tracking_percentages.min():.1f}% | Max: {tracking_percentages.max():.1f}%")
#     bins = np.histogram(tracking_percentages, bins=[0, 20, 40, 60, 80, 100])[0]
#     for i, count in enumerate(bins):
#         print(f"{i*20}-{(i+1)*20}%: {count} stars")
#     print("=======================================\n")


# def main(csv_path, output_folder, target, date):
#     os.makedirs(output_folder, exist_ok=True)
#     df = pd.read_csv(os.path.join(csv_path, tracked_csv))

#     # Debug: Initial inspection
#     print("\n===== DEBUG: INITIAL DATA INSPECTION =====")
#     print(f"Total rows: {len(df)}")
#     print("Columns:", df.columns.tolist())
#     print(f"Unique julian_date count: {df['julian_date'].nunique()}")
#     print("===========================================\n")

#     # Assign image_number based on julian_date
#     unique_dates = sorted(df['julian_date'].unique())
#     image_numbers_map = {jd: i + 1 for i, jd in enumerate(unique_dates)}
#     df['image_number'] = df['julian_date'].map(image_numbers_map)
#     total_images = len(unique_dates)
#     print(f"[Info] Total images: {total_images}")

#     # Extract quadrant number
#     df['quadrant_num'] = df['QUADRANT_FILE'].str.extract(r'_q(\d{1,2})').astype(float)

#     # Compute tracking percentage per star
#     non_nan_counts = df.groupby('star_index').apply(
#         lambda x: x[['ALPHA_J2000', 'DELTA_J2000']].notna().all(axis=1).sum(),
#         include_groups=False
#     )
#     df = df.merge(non_nan_counts.rename('tracking_count').to_frame(),
#                   left_on='star_index', right_index=True)
#     df['tracking_percentage'] = (df['tracking_count'] / total_images) * 100

#     summarize_tracking_percentages(df)

#     # Process each quadrant (0–15)
#     for q_num in range(16):
#         stars_in_quadrant = df[df['quadrant_num'] == q_num]

#         if stars_in_quadrant.empty:
#             print(f"[!] No stars in quadrant {q_num}. Skipping.")
#             continue

#         print(f">>> Quadrant {q_num}: {stars_in_quadrant['star_index'].nunique()} stars, "
#               f"{stars_in_quadrant['image_number'].nunique()} images")

#         plot_tracking_map(stars_in_quadrant, output_folder, f"quadrant_{q_num}", target, date, total_images)

#     print(f"[✓] All tracking map plots generated for 16 quadrants.")


# if __name__ == "__main__":
#     if len(sys.argv) != 5:
#         print("Usage: python plot_tracking_map.py <csv_dir> <output_dir> <target> <date>")
#         sys.exit(1)

#     csv_dir = sys.argv[1]
#     output_dir = sys.argv[2]
#     target = sys.argv[3]
#     date = sys.argv[4]

#     main(csv_path=csv_dir, output_folder=output_dir, target=target, date=date)
