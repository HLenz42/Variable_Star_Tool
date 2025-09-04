import pandas as pd
import numpy as np
import os
import sys

# === Configuration ===
input_path = sys.argv[1] + f"/tracked_stars_with_variability_flag_{sys.argv[3]}_{sys.argv[4]}.csv"
output_dir = sys.argv[2]
output_log = f"analysis_log_{sys.argv[3]}_{sys.argv[4]}.txt"
TIME_COLUMN = 'julian_date'
QUADRANT_COLUMN = 'QUADRANT_FILE'
STAR_INDEX_COLUMN = 'star_index'
RA_COLUMN = 'ALPHA_J2000'
DEC_COLUMN = 'DELTA_J2000'
MAG_COLUMN = 'MAG_AUTO'

# === Main Script ===
if __name__ == "__main__":
    print("Loading data...")
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        sys.exit(1)

    df[STAR_INDEX_COLUMN] = df[STAR_INDEX_COLUMN].astype(int)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, output_log)

    # Initialize log content
    log_lines = []

    # === General Dataset Summary ===
    log_lines.append(f"Log for {sys.argv[3]} observed on {sys.argv[4]}\n\n")
    log_lines.append("===== Dataset Summary =====")
    log_lines.append(f"Total rows: {len(df)}")
    log_lines.append(f"Unique stars (star_index): {df[STAR_INDEX_COLUMN].nunique()}")
    log_lines.append(f"Unique images (julian_date): {df[TIME_COLUMN].nunique()}")
    log_lines.append("==========================\n")

    # === NaN Diagnostics ===
    log_lines.append("===== NaN Diagnostics =====")
    for col in [RA_COLUMN, DEC_COLUMN, MAG_COLUMN, QUADRANT_COLUMN]:
        nan_count = df[col].isna().sum()
        nan_pct = 100 * nan_count / len(df)
        log_lines.append(f"{col} NaN count: {nan_count} ({nan_pct:.1f}%)")
    log_lines.append("==========================\n")

    # === Tracking Percentage ===
    total_images = df[TIME_COLUMN].nunique()
    non_nan_counts = df.groupby(STAR_INDEX_COLUMN).apply(
        lambda x: x[[RA_COLUMN, DEC_COLUMN, MAG_COLUMN]].notna().all(axis=1).sum(),
        include_groups=False
    )
    tracking_percentage = (non_nan_counts / total_images) * 100

    log_lines.append("===== Tracking Percentage Summary =====")
    log_lines.append(f"Mean tracking percentage: {tracking_percentage.mean():.1f}%")
    log_lines.append(f"Median tracking percentage: {tracking_percentage.median():.1f}%")
    log_lines.append(f"Min tracking percentage: {tracking_percentage.min():.1f}%")
    log_lines.append(f"Max tracking percentage: {tracking_percentage.max():.1f}%")
    bins = np.histogram(tracking_percentage, bins=[0, 20, 40, 60, 80, 100])[0]
    for i, count in enumerate(bins):
        log_lines.append(f"{i*20}-{(i+1)*20}%: {count} stars")
    log_lines.append(f"Stars with 100% tracking: {(tracking_percentage == 100).sum()}")
    log_lines.append("======================================\n")

    # === Quadrant Analysis ===
    df['quadrant_num'] = df[QUADRANT_COLUMN].str.extract(r'_q(\d{1,2})').astype(float)
    log_lines.append("===== Quadrant Analysis =====")
    for q_num in range(16):
        quad_df = df[df['quadrant_num'] == q_num]
        if quad_df.empty:
            log_lines.append(f"Quadrant {q_num}: No stars")
            continue
        star_count = quad_df[STAR_INDEX_COLUMN].nunique()
        log_lines.append(f"Quadrant {q_num}: {star_count} unique stars")
    log_lines.append("============================\n")

    # === First Image Analysis ===
    first_image_date = df[TIME_COLUMN].min()
    first_image_df = df[df[TIME_COLUMN] == first_image_date]
    log_lines.append("===== Star Indices in First Image by Quadrant =====")
    log_lines.append(f"First image julian_date: {first_image_date}")
    for q_num in range(16):
        quad_stars = first_image_df[first_image_df['quadrant_num'] == q_num][STAR_INDEX_COLUMN].unique()
        quad_stars_str = ', '.join(map(str, sorted(quad_stars)))
        log_lines.append(f"Quadrant {q_num}: {len(quad_stars)} stars\n{quad_stars_str}")
    log_lines.append("==============================================\n")

    # === Variable Stars ===
    variable_stars = df[df['is_variable'].str.lower() == 'yes'][STAR_INDEX_COLUMN].unique()
    log_lines.append("===== Variable Stars (is_variable == 'yes') =====")
    log_lines.append(f"Total variable stars: {len(variable_stars)}")
    log_lines.append(f"Variable star indices: {', '.join(map(str, sorted(variable_stars)))}")
    log_lines.append("\nVariable star types (Simbad_Type):")
    variable_types = df[df['is_variable'].str.lower() == 'yes']['Simbad_Type'].value_counts()
    for v_type, count in variable_types.items():
        log_lines.append(f"  {v_type}: {count} stars")
    log_lines.append("============================================\n")

    # === Exoplanet, Uncertain, and Unclassified Stars ===
    exoplanet_stars = df[df['Simbad_Type'] == 'EXOPLANET'][STAR_INDEX_COLUMN].unique()
    uncertain_stars = df[df['Simbad_Type'] == 'UNCERTAIN'][STAR_INDEX_COLUMN].unique()
    unclassified_stars = df[df['Simbad_Type'] == 'UNCLASSIFIED'][STAR_INDEX_COLUMN].unique()

    log_lines.append("===== Exoplanet Systems =====")
    log_lines.append(f"Total exoplanet stars: {len(exoplanet_stars)}")
    log_lines.append(f"Exoplanet star indices: {', '.join(map(str, sorted(exoplanet_stars)))}")
    log_lines.append("============================\n")

    log_lines.append("===== Uncertain Stars =====")
    log_lines.append(f"Total uncertain stars: {len(uncertain_stars)}")
    log_lines.append(f"Uncertain star indices: {', '.join(map(str, sorted(uncertain_stars)))}")
    log_lines.append("==========================\n")

    log_lines.append("===== Unclassified Stars =====")
    log_lines.append(f"Total unclassified stars: {len(unclassified_stars)}")
    log_lines.append(f"Unclassified star indices: {', '.join(map(str, sorted(unclassified_stars)))}")
    log_lines.append("=============================\n")

    # === Simbad_Type Distribution ===
    log_lines.append("===== Simbad_Type Distribution =====")
    simbad_types = df['Simbad_Type'].value_counts()
    for s_type, count in simbad_types.items():
        log_lines.append(f"{s_type}: {count} stars")
    log_lines.append("==================================\n")

    # === Save Log ===
    with open(log_path, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Analysis log saved to: {log_path}")