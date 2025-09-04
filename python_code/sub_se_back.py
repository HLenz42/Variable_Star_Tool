import sys
import glob
import os
from astropy import units as u
from ccdproc import CCDData
from astropy.io import fits

def subtract_background(wcs_folder, backgrounds_folder, output_folder):
    """
    Subtract Source Extractor background images from WCS-calibrated FITS images,
    preserving the WCS header from the input image.

    Args:
        wcs_folder (str): Directory containing WCS-calibrated FITS images (e.g., Wasp12b_20-001b_bdf_q2_wcs.fit).
        backgrounds_folder (str): Directory containing background FITS images (e.g., Wasp12b_20-001b_bdf_q2_wcs_back.fit).
        output_folder (str): Directory to save background-subtracted FITS images.
    """
    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Get WCS image files
    wcs_files = sorted(
        glob.glob(f"{wcs_folder}/*.fit") +
        glob.glob(f"{wcs_folder}/*.fits") +
        glob.glob(f"{wcs_folder}/*.FITS")
    )
    if not wcs_files:
        print(f"Error: No WCS FITS files found in {wcs_folder}")
        sys.exit(1)
    print(f">>> Found {len(wcs_files)} WCS FITS files in {wcs_folder}")

    # Get background files
    back_files = sorted(
        glob.glob(f"{backgrounds_folder}/*.fit") +
        glob.glob(f"{backgrounds_folder}/*.fits") +
        glob.glob(f"{backgrounds_folder}/*.FITS")
    )
    if not back_files:
        print(f"Error: No background FITS files found in {backgrounds_folder}")
        sys.exit(1)
    print(f">>> Found {len(back_files)} background FITS files in {backgrounds_folder}")

    # Create a dictionary mapping base names to background file paths
    back_dict = {}
    for back_file in back_files:
        base = os.path.basename(back_file)
        # Remove '_back.fit' or similar to get the base name
        base_name = base.replace('_back.fit', '').replace('_back.fits', '').replace('_back.FITS', '')
        back_dict[base_name] = back_file

    # Process each WCS image
    for wcs_file in wcs_files:
        base = os.path.basename(wcs_file)
        base_name, ext = os.path.splitext(base)
        print(f"\n>>> Processing {base}")

        # Find corresponding background file
        back_file = back_dict.get(base_name)
        if not back_file:
            print(f"Warning: No background file found for {base}. Skipping.")
            continue

        # Read WCS image and its original header
        try:
            with fits.open(wcs_file) as hdul:
                wcs_header = hdul[0].header  # Store original header
                wcs_img = CCDData.read(wcs_file, unit='adu')
        except Exception as e:
            print(f"Warning: Could not read WCS image {wcs_file}. Skipping. Error: {e}")
            continue

        # Read background image
        try:
            back_img = CCDData.read(back_file, unit='adu')
        except Exception as e:
            print(f"Warning: Could not read background image {back_file}. Skipping. Error: {e}")
            continue

        # Verify dimensions match
        if wcs_img.shape != back_img.shape:
            print(f"Warning: Dimension mismatch between {wcs_file} ({wcs_img.shape}) and {back_file} ({back_img.shape}). Skipping.")
            continue

        # Subtract background
        try:
            result_img = wcs_img.subtract(back_img)
        except Exception as e:
            print(f"Warning: Background subtraction failed for {wcs_file}. Skipping. Error: {e}")
            continue

        # Create output filename
        output_filename = f"{base_name}{ext}"
        output_path = os.path.join(output_folder, output_filename)

        # Save result as single-HDU FITS with original WCS header
        try:
            primary_hdu = fits.PrimaryHDU(data=result_img.data, header=wcs_header)
            hdul = fits.HDUList([primary_hdu])
            hdul.writeto(output_path, overwrite=True)
            print(f">>> Saved background-subtracted image with WCS: {output_path}")
        except Exception as e:
            print(f"Warning: Could not save {output_path}. Skipping. Error: {e}")
            continue

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python subtract_se_background_wcs_preserved.py <wcs_folder> <backgrounds_folder> <output_folder>")
        sys.exit(1)

    wcs_folder = sys.argv[1]
    backgrounds_folder = sys.argv[2]
    output_folder = sys.argv[3]

    subtract_background(wcs_folder, backgrounds_folder, output_folder)