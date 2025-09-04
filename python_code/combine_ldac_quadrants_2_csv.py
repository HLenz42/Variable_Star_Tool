import os
import glob
import pandas as pd
import numpy as np
from astropy.io import fits
from collections import defaultdict
import sys

def convert_ldac_to_csv(ldac_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    log_file = os.path.join(output_folder, f"combine_ldac_log_{os.path.basename(ldac_folder)}.txt")
    log_lines = ["===== Combine LDAC to CSV Log ====="]

    # Find LDAC files
    ldac_files = sorted(glob.glob(os.path.join(ldac_folder, "*.ldac")))
    if not ldac_files:
        log_lines.append(f"ERROR: No .ldac files found in: {ldac_folder}")
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_lines))
        print(f"ERROR: No .ldac files found in: {ldac_folder}")
        return

    # Group files by base image name (before '_q#')
    groups = defaultdict(list)
    for f in ldac_files:
        fname = os.path.basename(f)
        if "_q" in fname:
            base = fname.split("_q")[0]
            groups[base].append(f)
        else:
            log_lines.append(f"[!] Skipping non-quadrant file: {fname}")

    for base_name, files in groups.items():
        print(f"\nProcessing {base_name} with {len(files)} quadrants")
        log_lines.append(f"\nProcessing {base_name} with {len(files)} quadrants")
        combined_df = pd.DataFrame()

        for f in sorted(files):
            try:
                with fits.open(f) as hdul:
                    data = hdul[2].data
                    if data is None:
                        raise ValueError("No data in HDU[2]")
                    
                    # Convert FITS data to DataFrame, handling VIGNET separately
                    df_data = {}
                    for col in data.names:
                        if col == 'VIGNET':
                            # Format VIGNET as ((row1;row2;...);(row3;...);...)
                            vignet_data = data['VIGNET']  # Shape: (n_sources, height, width)
                            vignet_strings = []
                            for v in vignet_data:
                                # Convert each row to semicolon-separated string
                                row_strings = [';'.join(map(str, row)) for row in v]
                                # Join rows with parentheses
                                vignet_str = '(' + ';'.join(f'({row})' for row in row_strings) + ')'
                                vignet_strings.append(vignet_str)
                            df_data['VIGNET'] = vignet_strings
                        else:
                            df_data[col] = data[col]
                    
                    df = pd.DataFrame(df_data)
                    df["QUADRANT_FILE"] = os.path.basename(f)
                    combined_df = pd.concat([combined_df, df], ignore_index=True)
                    print(f"   Added: {os.path.basename(f)} ({len(df)} sources)")
                    log_lines.append(f"   Added: {os.path.basename(f)} ({len(df)} sources)")
                    if 'VIGNET' in data.names:
                        log_lines.append(f"   VIGNET shape: {data['VIGNET'].shape}")
            except Exception as e:
                print(f"   ERROR: Failed to read {f}: {e}")
                log_lines.append(f"   ERROR: Failed to read {f}: {e}")

        if not combined_df.empty:
            # Add INDEX column
            combined_df.insert(0, "INDEX", range(len(combined_df)))
            # Ensure consistent column types
            for col in combined_df.columns:
                if col not in ['QUADRANT_FILE', 'VIGNET']:
                    combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')

            # Save CSV with comma separator
            output_path = os.path.join(output_folder, f"{base_name}.csv")
            combined_df.to_csv(output_path, sep=',', index=False)
            print(f">>> Saved combined CSV: {output_path} ({len(combined_df)} sources)")
            log_lines.append(f">>> Saved combined CSV: {output_path} ({len(combined_df)} sources)")
        else:
            print(f"ERROR: No data extracted for {base_name}")
            log_lines.append(f"ERROR: No data extracted for {base_name}")

    # Save log
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f">>> Saved log: {log_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python combine_ldac_quadrants_2_csv.py [ldac_folder] [output_folder]")
        sys.exit(1)

    ldac_folder = sys.argv[1]
    output_folder = sys.argv[2]
    convert_ldac_to_csv(ldac_folder, output_folder)

# import os
# import glob
# import pandas as pd
# from astropy.io import fits
# from collections import defaultdict
# import sys 

# def convert_ldac_to_csv(ldac_folder, output_folder):
#     os.makedirs(output_folder, exist_ok=True)

#     ldac_files = sorted(glob.glob(os.path.join(ldac_folder, "*.ldac")))
#     if not ldac_files:
#         print(f"[!] No .ldac files found in: {ldac_folder}")
#         return

#     # Group files by base image name (everything before '_q#')
#     groups = defaultdict(list)
#     for f in ldac_files:
#         fname = os.path.basename(f)
#         if "_q" in fname:
#             base = fname.split("_q")[0]
#             groups[base].append(f)

#     for base_name, files in groups.items():
#         print(f"\n[⧗] Processing {base_name} with {len(files)} quadrants")

#         combined_df = pd.DataFrame()

#         for f in sorted(files):
#             try:
#                 with fits.open(f) as hdul:
#                     data = hdul[2].data
#                     if data is None:
#                         raise ValueError("No data in HDU[2]")
#                     df = pd.DataFrame(data)  # No newbyteorder() needed
#                     df["QUADRANT_FILE"] = os.path.basename(f)
#                     combined_df = pd.concat([combined_df, df], ignore_index=True)
#                     print(f"   [✓] Added: {os.path.basename(f)}")
#             except Exception as e:
#                 print(f"   [✗] Failed to read {f}: {e}")

#         if not combined_df.empty:
#             combined_df.insert(0, "INDEX", range(len(combined_df)))  # Add new column
#             output_path = os.path.join(output_folder, base_name + ".csv")
#             combined_df.to_csv(output_path, index=False)
#             print(f">>> Saved combined CSV: {output_path}")
#         else:
#             print(f"[!] No data extracted for {base_name}")

# if __name__ == "__main__":
#     # Update these paths as needed:
#     ldac_folder = sys.argv[1]
#     output_folder = sys.argv[2]

#     convert_ldac_to_csv(ldac_folder, output_folder)


