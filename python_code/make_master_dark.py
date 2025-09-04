import sys
import os
import glob
from astropy import units as u
from ccdproc import CCDData, combine
import numpy as np

def make_master_dark(dark_folder, master_bias_path, output_folder, output_name):
    dark_files = sorted(glob.glob(f"{dark_folder}/*.fit") + 
                        glob.glob(f"{dark_folder}/*.fits") +
                        glob.glob(f"{dark_folder}/*.FITS"))

    print("dark_files:")
    print(f"{dark_files}")
    print()
    if not dark_files:
        print(f"Error: Not FITS files found in: {dark_folder}")
        print()
        sys.exit(1)


    bias = CCDData.read(master_bias_path, unit='adu', format='fits')
    dark_ccds = [CCDData.read(f, unit='adu', format='fits').subtract(bias) for f in dark_files]

    master_dark = combine(dark_ccds, method='median', sigma_clip=True, unit=u.adu)

        # Get dark exposure time
    dark_exptime = CCDData.read(dark_files[0], unit='adu', format='fits').header['EXPTIME']
    
    # Get image exposure time from environment variable (set by main script)
    image_exptime = float(os.environ.get("IMAGE_EXPTIME", dark_exptime))

    if dark_exptime != image_exptime:
        scale_factor = image_exptime / dark_exptime
        print()
        print(f">>> Scaling master dark by factor {scale_factor:.2f} to match image exposure.")
        master_dark = master_dark.multiply(scale_factor)
    
    master_dark.header['EXPTIME'] = image_exptime  # Set to match image

    master_dark.data=master_dard.data.astype(np.uint32)

    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_name)
    master_dark.write(output_path, overwrite=True)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python make_master_dark.py [dark_folder] [master_bias_path] [output_folder] [output_name]")
        sys.exit(1)

    dark_folder = sys.argv[1]
    master_bias_path = sys.argv[2]
    output_folder = sys.argv[3]
    output_name = sys.argv[4]

    make_master_dark(dark_folder, master_bias_path, output_folder, output_name)
