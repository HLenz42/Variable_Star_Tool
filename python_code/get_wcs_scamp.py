# wcs_pipeline.py

import os
import glob
import subprocess
import numpy as np
import sys
from astropy.io import fits
from tqdm import tqdm

# ---------- Config -----------
input_fits_dir = sys.argv[1]      # Path to processed images
output_fits_dir = sys.argv[2]     # Path to output WCS-ed images
scamp_config_dir = sys.argv[3]    # Path to SCAMP config dir
scamp_output_dir = sys.argv[4]    # Path to SCAMP output dir
se_config_dir = sys.argv[5]       # Path to SExtractor config
ldac_dir = sys.argv[6]            # Path to SExtractor output (.ldac files)
path_2_swarp_config = sys.argv[7]    # Path to SWarp config dir


# Create directories if they don't exist
os.makedirs(input_fits_dir, exist_ok=True)
os.makedirs(output_fits_dir, exist_ok=True)
os.makedirs(ldac_dir, exist_ok=True)
os.makedirs(se_config_dir, exist_ok=True)
os.makedirs(scamp_config_dir, exist_ok=True)
os.makedirs(scamp_output_dir, exist_ok=True)

# SExtractor configuration files
se_config_file = os.path.join(se_config_dir, "kgmt.config")
se_param_file = os.path.join(se_config_dir, "kgmt.param")
se_filter_file = os.path.join(se_config_dir, "default.conv")

# Verify configuration files exist
for config_file in [se_config_file, se_param_file, se_filter_file]:
    if not os.path.exists(config_file):
        print(f"Error: SExtractor config file not found: {config_file}")
        sys.exit(1)

# SExtractor parameters
se_command = "source-extractor"
im_gain = "1.42"        # e/ADU
im_pixscale = "3.92"    # arcsec/pixel
satur_level = "65535"

# =============== Step 1: Run SExtractor to generate .ldac catalogs ===============
print("\n>>> Run SExtractor or skip to next step?")
print("(1) Run SExtractor")
print("(2) Use previous runs output")
choice = int(input("Answer: "))

run_se = False
if choice == 1:
    run_se = True
    print(">>> Running SExtractor\n")
elif choice == 2:
    print(">>> Skipping to next step\n")
else:
    print("Error: Please select 1 or 2")
    sys.exit(1)

if run_se:
    # Gather FITS files
    im_files = sorted(glob.glob(os.path.join(input_fits_dir, '*.[fF][iI][tT][sS]')) + 
    				  glob.glob(os.path.join(input_fits_dir, '*.[fF][iI][tT]')))
    if not im_files:
        print(f"Error: No FITS files found in {input_fits_dir}")
        sys.exit(1)

    all_se_cat_paths = []
    for image in tqdm(im_files, desc="Running SExtractor"):
        base = os.path.splitext(os.path.basename(image))[0]
        output_cat = os.path.join(ldac_dir, f"{base}_se.ldac")

        # SExtractor command
        command = [
            se_command, image,
            "-c", se_config_file,
            "-CATALOG_NAME", output_cat,
            "-CATALOG_TYPE", "FITS_LDAC",
            "-PARAMETERS_NAME", se_param_file,
            "-FILTER_NAME", se_filter_file,
            "-SATUR_LEVEL", satur_level,
            "-GAIN", im_gain,
            "-PIXEL_SCALE", im_pixscale,
            "-CHECKIMAGE_TYPE", "NONE"
        ]

        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            if os.path.exists(output_cat):
                print(f">>> SExtractor catalog created: {output_cat}")
                all_se_cat_paths.append(output_cat)
            else:
                print(f">>> Warning: SExtractor failed to create {output_cat}")
        except subprocess.CalledProcessError as e:
            print(f">>> Error: SExtractor failed for {image}: {e.stderr}")

    # Optionally create ldac_list.txt for SCAMP (if needed)
    ldac_list_path = os.path.join(scamp_config_dir, "ldac_list.txt")
    with open(ldac_list_path, "w") as f:
        for cat_path in all_se_cat_paths:
            f.write(os.path.abspath(cat_path) + "\n")
    print(f">>> Created SExtractor catalog list for SCAMP: {ldac_list_path}")

	# print(f"\n>>> SExtractor complete: {len(all_se_cat_paths)} catalogs created.\n")

# =============== Step 2: Run SCAMP on all .ldac catalogs ===============
print(">>> Run SCAMP or reuse previous outputs?")
print("(1) Run SCAMP")
print("(2) Reuse outputs")
run_scamp_choice = int(input("Answer: "))

run_scamp = run_scamp_choice == 1
if run_scamp:
    print(">>> Running SCAMP...")
    while True:
        subprocess.run([
            "python3", "scamp.py",
            input_fits_dir,    # Fits we just applyed the rough WCS to
            ldac_dir,	 	   # Output of Source Extractor
            scamp_config_dir,  # Scamp Config
            scamp_output_dir   # Scamp Output
        ])

        print("\n>>> SCAMP complete.")
        print(">>> Look at the checkplots & star cat plots")
        print(">>> Run SCAMP again?")
        print("(1) Yes")
        print("(2) No")
        print("(3) Exit")
        choice = int(input("Answer: "))

        if choice == 1:
            print(">>> Going back...")
        elif choice == 2:
            print(">>> Going to the next step...")
            break
        elif choice == 3:
            print(">>> Goodbye!")
            sys.exit(0)
        else:
            print("Error: Please select 1, 2, or 3")

else:
    print(">>> Skipping to next step...")

# ============ Step 3: Apply .ahead WCS headers to FITS images ============
print(">>> Apply generated .ahead files to bdf processed fits headers?")
print("(1) yes")
print("(2) Reuse outputs")
run_ahead_choice = int(input("Answer: "))

if run_ahead_choice == 1:
    print(">>> Applying WCS headers in-place...")

    im_files = sorted(glob.glob(os.path.join(input_fits_dir, '*.[fF][iI][tT][sS]')) + 
                          glob.glob(os.path.join(input_fits_dir, '*.[fF][iI][tT]')))
    if not im_files:
        print(f"Error: No FITS files found in {input_fits_dir}")
        sys.exit(1)

    fits_files = im_files
    input_folder = input_fits_dir
    ahead_folder = scamp_output_dir

    for fits_file in tqdm(fits_files, desc="Applying WCS headers"):
        # Handle _wcs suffix if add_rough_wcs.py used --no-overwrite
        base = os.path.splitext(os.path.basename(fits_file))[0]
        if base.endswith('_wcs'):
            base = base[:-4]  # Remove _wcs suffix to match scamp.py's .ldac and .ahead files

        ahead_filename = base + ".ahead"
        ahead_path = os.path.join(ahead_folder, ahead_filename)

        if not os.path.exists(ahead_path):
            print(f">>> Missing .ahead file: {ahead_filename}, skipping...")
            continue

        input_path = os.path.join(input_folder, fits_file)
        try:
            hdul = fits.open(input_path, mode='update')  # Open in update mode for in-place modification

            # Delete existing WCS keywords before applying SCAMP's
            wcs_keywords_to_delete = [
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                'CTYPE1', 'CTYPE2', 'WCSAXES', 'WCSINFO',
                'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
                'CROTA1', 'CROTA2',
                # SCAMP distortion keywords
                'PV1_0', 'PV1_1', 'PV1_2', 'PV1_3', 'PV1_4', 'PV1_5',
                'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_4', 'PV2_5',
                'A_ORDER', 'B_ORDER', 'AP_ORDER', 'BP_ORDER',
                'A_0_0', 'A_0_1', 'A_1_0', 'A_0_2', 'A_1_1', 'A_2_0',
                'B_0_0', 'B_0_1', 'B_1_0', 'B_0_2', 'B_1_1', 'B_2_0'
            ]

            for key in wcs_keywords_to_delete:
                if key in hdul[0].header:
                    del hdul[0].header[key]

            # Read and apply .ahead file
            try:
                with open(ahead_path, "r") as f:
                    header_lines = f.readlines()

                for line in header_lines:
                    if line.strip() and "=" in line:
                        key_val = line.split("=")
                        key = key_val[0].strip()
                        val = key_val[1].split("/")[0].strip()
                        try:
                            hdul[0].header[key] = eval(val)
                        except:
                            hdul[0].header[key] = val.strip("'")
            except Exception as e:
                print(f">>> Error reading .ahead file {ahead_filename}: {e}")
                hdul.close()
                continue

            try:
                hdul.flush()  # Save changes to the original file
                print(f">>> Updated WCS in-place: {input_path}")
            except Exception as e:
                print(f">>> Error updating FITS file {input_path}: {e}")
            finally:
                hdul.close()

        except Exception as e:
            print(f">>> Error processing FITS file {fits_file}: {e}")

# ============ Step 4: Run SWarp on each of the images to 'flatten' them ============
print(">>> Run SWarp to create final WCS'ed images?")
print("(1) yes")
print("(2) Reuse outputs")
run_swarp_choice = int(input("Answer: "))

if run_swarp_choice == 1:

    print(">>> Running swarp.py ...")

    subprocess.run([
            "python3", "swarp.py",
            input_fits_dir,      # Fits we just applyed the .ahead wcs to
            path_2_swarp_config, # SWarp Config
            scamp_output_dir,    # Scamp Output
            output_fits_dir      # SWarp Output
        ])

    print()
    print(">>> SWarp Complete")

else:
    print(">>> Going to next step...")
