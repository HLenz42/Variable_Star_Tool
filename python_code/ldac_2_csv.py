import os
import glob
import pandas as pd
from astropy.io import fits
from tqdm import tqdm
import sys

def convert_ldac_to_csv(ldac_folder, output_folder):
    """
    Convert .ldac files in ldac_folder to CSV and save to output_folder.
    Each .ldac file is converted to a separate CSV file with an INDEX column.
    Output filenames remove '_se' from the input filename.
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Find all .ldac files
    ldac_files = sorted(glob.glob(os.path.join(ldac_folder, "*.ldac")))
    if not ldac_files:
        print(f"[!] No .ldac files found in: {ldac_folder}")
        return

    print(f"[⧗] Found {len(ldac_files)} .ldac files to process.")

    # Process each .ldac file
    for ldac_file in tqdm(ldac_files, desc="Converting LDAC to CSV"):
        # Get base name and remove '_se' suffix if present
        base_name = os.path.splitext(os.path.basename(ldac_file))[0]
        if base_name.endswith('_wcs_se'):
            base_name = base_name[:-7]  # Remove '_wcs_se'
        output_path = os.path.join(output_folder, f"{base_name}.csv")

        try:
            with fits.open(ldac_file) as hdul:
                # Check if HDU[2] contains data (typically the catalog table)
                if len(hdul) < 3 or hdul[2].data is None:
                    raise ValueError("No data in HDU[2]")
                
                # Convert FITS table to pandas DataFrame
                df = pd.DataFrame(hdul[2].data)
                
                # Add INDEX column as the first column
                df.insert(0, "INDEX", range(len(df)))
                
                # Add source file name as a column
                df["LDAC_FILE"] = os.path.basename(ldac_file)
                
                # Save to CSV
                df.to_csv(output_path, index=False)
                print(f"[✓] Converted {os.path.basename(ldac_file)} to {output_path}")

        except Exception as e:
            print(f"[✗] Failed to convert {os.path.basename(ldac_file)}: {e}")

    print(f"\n[✓] Conversion complete. Processed {len(ldac_files)} files.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ldac_to_csv.py <ldac_folder> <output_folder>")
        sys.exit(1)

    ldac_folder = sys.argv[1]
    output_folder = sys.argv[2]

    # Verify input folder exists
    if not os.path.isdir(ldac_folder):
        print(f"[!] Input folder does not exist: {ldac_folder}")
        sys.exit(1)

    convert_ldac_to_csv(ldac_folder, output_folder)