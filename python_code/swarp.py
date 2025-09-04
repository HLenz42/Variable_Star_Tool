# swarp.py

import glob
import subprocess
import os
import sys
from tqdm import tqdm
import numpy as np

# ----------------- Read from command line -----------------
if len(sys.argv) != 5:
    print("Usage: python swarp.py <sci_im_dir> <swarp_config_dir> <scamp_output_dir> <swarp_output_dir>")
    sys.exit(1)

sci_im_dir = sys.argv[1]          # Where are the science images?
swarp_config_dir = sys.argv[2]    # Path to SWarp config (for swarp.conf)
scamp_output_dir = sys.argv[3]    # Path to SCAMP output (for _stars.head files)
swarp_output_dir = sys.argv[4]    # Path to SWarp output (final WCS-aligned images)
fits_pattern = "*.fit"  # Default to *.fit

# Verify directories
for d in [sci_im_dir, swarp_config_dir, scamp_output_dir, swarp_output_dir]:
    if not os.path.isdir(d):
        print(f"Error: Directory not found: {d}")
        sys.exit(1)

# Check write permissions for swarp_output_dir
if not os.access(swarp_output_dir, os.W_OK):
    print(f"Error: No write permissions for {swarp_output_dir}")
    sys.exit(1)

# Verify SWarp config file
swarp_config_file = os.path.join(swarp_config_dir, "swarp.conf")
if not os.path.exists(swarp_config_file):
    print(f"Error: SWarp config file not found: {swarp_config_file}")
    sys.exit(1)

# ----------------- Constants -----------------
swarp_command = "SWarp"  # Lowercase for case-sensitive systems
pixel_scale_arcsec = 3.92  # Match input pixel scale from scamp.py
os.makedirs(swarp_output_dir, exist_ok=True)  # Ensure output directory exists

# ----------------- Image list -----------------
im_files = np.sort(glob.glob(os.path.join(sci_im_dir, fits_pattern)))
if len(im_files) == 0:
    print(f"Error: No FITS files found in {sci_im_dir} with pattern {fits_pattern}")
    sys.exit(1)

# ----------------- Run SWarp on each image individually -----------------
print("\n>>> Running SWarp on individual images:")
for image in tqdm(im_files, desc="Running SWarp"):
    base = os.path.splitext(os.path.basename(image))[0]
    head_file = os.path.join(scamp_output_dir, base + "_stars.head")  # Use _stars.head
    swarp_output_file = os.path.join(swarp_output_dir, base + "_wcs.fit")

    if not os.path.exists(head_file):
        print(f">>> No _stars.head file found for {image} at {head_file}, skipping SWarp...")
        continue

    # Run SWarp command
    swarp_command_list = [
        swarp_command, image,
        "-c", swarp_config_file,
        "-IMAGEOUT_NAME", swarp_output_file,
        "-WEIGHTOUT_NAME", "none",  # Explicitly disable weight file output
        "-WEIGHT_TYPE", "NONE",  # Disable weight map generation
        "-PIXEL_SCALE", str(pixel_scale_arcsec),
        "-RESAMPLE", "Y",  # Enable resampling of the input image to the WCS defined in the _stars.head file
        "-RESAMPLING_TYPE", "LANCZOS3",  # Specifies resampling kernel Lanczos-3 for high-quality interpolation
        "-RESAMPLE_DIR", swarp_output_dir,
        "-PROJECTION_TYPE", "TAN",  # Specifies the gnomonic tangent plane projection used in the _stars.head files
        "-COMBINE", "N",  # Disables combining multiple images into one output
        "-SUBTRACT_BACK", "N",  # Disables background subtraction for differential photometry
        "-HEADER_SUFFIX", "_stars.head",  # Use _stars.head files
        "-RESAMPLE_SUFFIX", "_wcs.fit"  # Prevent .resamp suffix
    ]

    print(f">>> Running SWarp for {image}:")
    print(" ".join(swarp_command_list))
    try:
        result = subprocess.run(swarp_command_list, capture_output=True, text=True)
        if result.returncode != 0:
            print(f">>> SWarp failed for {image}: {result.stderr}")
            continue
        if os.path.exists(swarp_output_file):
            print(f">>> Wrote SWarp output: {swarp_output_file}")
        else:
            print(f">>> Failed to write SWarp output: {swarp_output_file}")
            print(f">>> Check if {swarp_output_dir} is writable and swarp.conf settings")
    except Exception as e:
        print(f">>> Error running SWarp for {image}: {e}")