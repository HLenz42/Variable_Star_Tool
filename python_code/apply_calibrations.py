# apply_calibrations.py
import sys
import glob
import os
from astropy import units as u
from ccdproc import CCDData, subtract_bias, subtract_dark, flat_correct
from astropy.io import fits

def apply_calibrations(raw_folder, master_bias_path, master_dark_path, master_flat_path, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    # Read master calibration frames
    try:
        bias = CCDData.read(master_bias_path, unit='adu', format='fits')
        dark = CCDData.read(master_dark_path, unit='adu', format='fits')
    except Exception as e:
        print(f"Error: Could not read master calibration files. Bias: {master_bias_path}, Dark: {master_dark_path}. Error: {e}")
        sys.exit(1)

    # Get master dark exposure time
    dark_exptime = float(dark.header.get('EXPTIME', 1.0))  # Default to 1.0 if missing
    print(f">>> Master dark exposure time: {dark_exptime:.2f} s")

    # Read master flat
    flat_files = sorted(glob.glob(f"{master_flat_path}/*.fit") +
                        glob.glob(f"{master_flat_path}/*.fits") +
                        glob.glob(f"{master_flat_path}/*.FITS"))
    if not flat_files:
        print(f"Error: No flat files found in {master_flat_path}")
        sys.exit(1)

    try:
        flat = CCDData.read(flat_files[0], unit='adu', format='fits')
    except Exception as e:
        print(f"Error: Could not read flat file {flat_files[0]}. Error: {e}")
        sys.exit(1)

    # Get raw image files
    raw_files = sorted(glob.glob(f"{raw_folder}/*.fit") +
                       glob.glob(f"{raw_folder}/*.fits") +
                       glob.glob(f"{raw_folder}/*.FITS"))
    if not raw_files:
        print(f"Error: No raw files found in {raw_folder}")
        sys.exit(1)

    for raw_file in raw_files:
        print(f"\n>>> Processing {raw_file}")
        try:
            img = CCDData.read(raw_file, unit='adu', format='fits')
        except Exception as e:
            print(f"Warning: Could not read {raw_file}. Skipping. Error: {e}")
            continue

        # Get image exposure time
        image_exptime = float(img.header.get('EXPTIME', os.environ.get('IMAGE_EXPTIME', dark_exptime)))
        print(f">>> Image exposure time: {image_exptime:.2f} s")

        # Scale master dark if needed
        scaled_dark = dark.copy()  # Create a copy to avoid modifying the original
        if dark_exptime != image_exptime:
            scale_factor = image_exptime / dark_exptime
            print(f">>> Scaling master dark by factor {scale_factor:.2f} to match image exposure")
            scaled_dark = scaled_dark.multiply(scale_factor)
            scaled_dark.header['EXPTIME'] = image_exptime  # Update header
        else:
            print(">>> No scaling needed (dark and image exposure times match)")

        # Apply calibrations
        try:
            img = subtract_bias(img, bias)
            img = subtract_dark(img, scaled_dark, exposure_time='EXPTIME', exposure_unit=u.s, scale=False)  # Scaling already done
            img = flat_correct(img, flat)
        except Exception as e:
            print(f"Warning: Calibration failed for {raw_file}. Skipping. Error: {e}")
            continue

        # Modify filename to add "_bdf"
        base = os.path.basename(raw_file)
        name, ext = os.path.splitext(base)
        output_filename = f"{name}_bdf{ext}"
        output_path = os.path.join(output_folder, output_filename)

        # Save as single-HDU FITS
        try:
            primary_hdu = fits.PrimaryHDU(data=img.data, header=img.header)
            hdul = fits.HDUList([primary_hdu])
            hdul.writeto(output_path, overwrite=True)
            print(f">>> Saved calibrated image: {output_path}")
        except Exception as e:
            print(f"Warning: Could not save {output_path}. Skipping. Error: {e}")
            continue

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python apply_calibrations.py [raw_folder] [master_bias_path] [master_dark_path] [master_flat_path] [output_folder]")
        sys.exit(1)

    raw_folder = sys.argv[1]
    master_bias_path = sys.argv[2]
    master_dark_path = sys.argv[3]
    master_flat_path = sys.argv[4]
    output_folder = sys.argv[5]

    apply_calibrations(raw_folder, master_bias_path, master_dark_path, master_flat_path, output_folder)