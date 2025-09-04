import sys
import os
import glob
import numpy as np
from astropy import units as u
from ccdproc import CCDData, combine

def make_master_bias(bias_folder, output_folder, output_name):
    bias_files = sorted(glob.glob(f"{bias_folder}/*.fit") +     # Get all fits in the bias folder 
    					glob.glob(f"{bias_folder}/*.fits") +
    					glob.glob(f"{bias_folder}/*.FITS"))

    if not bias_files:
    	print(f"Error: No FITS files found in: {bias_folder}")
    	print()
    	sys.exit(1)

    bias_ccds = [CCDData.read(f, unit='adu', format='fits') for f in bias_files] # read all found fits 

    master_bias = combine(bias_ccds, method='median', sigma_clip=True, unit=u.adu) # median combine all loaded fits 

    master_bias.data = master_bias.data.astype(np.uint32)

    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_name)
    master_bias.write(output_path, overwrite=True) # Write combined fit

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python make_master_bias.py [bias_folder] [output_folder] [output_name]")
        sys.exit(1)

    bias_folder = sys.argv[1]
    output_folder = sys.argv[2]
    output_name = sys.argv[3]

    make_master_bias(bias_folder, output_folder, output_name)

