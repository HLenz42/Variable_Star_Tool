import os
import glob
import argparse
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import name_resolve
from tqdm import tqdm

def get_target_coords_from_simbad(target_name):
    """Query SIMBAD for the target's RA/Dec coordinates."""
    try:
        coord = SkyCoord.from_name(target_name)
        print(f">>> Found target '{target_name}': RA = {coord.ra.deg:.6f}, Dec = {coord.dec.deg:.6f}")
        return coord.ra.deg, coord.dec.deg
    except name_resolve.NameResolveError:
        raise ValueError(f">>> Could not resolve target name '{target_name}' using SIMBAD.")

def add_rough_wcs_to_image(header, ra_deg, dec_deg, x_center, y_center, pixel_scale_deg):
    """Add rough WCS to the FITS header."""
    hdr = header.copy()
    hdr['CRVAL1'] = ra_deg
    hdr['CRVAL2'] = dec_deg
    hdr['CRPIX1'] = x_center
    hdr['CRPIX2'] = y_center
    hdr['CD1_1'] = -pixel_scale_deg  # Negative for RA increasing leftward
    hdr['CD1_2'] = 0.0
    hdr['CD2_1'] = 0.0
    hdr['CD2_2'] = pixel_scale_deg
    hdr['CTYPE1'] = 'RA---TAN'
    hdr['CTYPE2'] = 'DEC--TAN'
    hdr['PIXELSCALE'] = pixel_scale_deg * 3600  # Store pixel scale in arcsec/pixel
    hdr['WCSINFO'] = 'Rough WCS from SIMBAD target for SCAMP preprocessing'
    return hdr

def process_images(input_folder, ra_center, dec_center, pixel_scale_arcsec=3.92, overwrite=True):
    """Add rough WCS to all FITS files in the input folder."""
    image_paths = sorted(glob.glob(os.path.join(input_folder, '*.[fF][iI][tT][sS]')) +
                         glob.glob(os.path.join(input_folder, '*.[fF][iI][tT]')))
    if not image_paths:
        print(">>> No FITS files found in input folder.")
        return

    print(f">>> Found {len(image_paths)} FITS files.")
    for path in tqdm(image_paths, desc="Processing FITS files"):
        base_name = os.path.splitext(os.path.basename(path))[0]
        print(f">>> Processing image: {base_name}")

        try:
            with fits.open(path, mode='update' if overwrite else 'readonly') as hdul:
                data = hdul[0].data
                header = hdul[0].header
                ny, nx = data.shape

                # Check for binning in the FITS header
                binning_x = header.get('XBINNING', 1)
                binning_y = header.get('YBINNING', 1)
                if binning_x != binning_y:
                    print(f">>> Warning: Asymmetric binning (XBINNING={binning_x}, YBINNING={binning_y}) in {base_name}")
                binning = max(binning_x, binning_y)  # Use the larger binning factor
                effective_pixel_scale = pixel_scale_arcsec * binning
                pixel_scale_deg = effective_pixel_scale / 3600.0
                print(f">>> Using effective pixel scale: {effective_pixel_scale:.2f} arcsec/pixel (binning={binning})")

                # Set reference pixel to image center
                x_center = nx / 2
                y_center = ny / 2

                # Add rough WCS to header
                new_header = add_rough_wcs_to_image(header, ra_center, dec_center, x_center, y_center, pixel_scale_deg)

                # Update the FITS file
                if overwrite:
                    hdul[0].header = new_header
                    hdul.flush()  # Save changes to the original file
                    print(f">>> Updated {base_name} with rough WCS (pixel scale={effective_pixel_scale:.2f} arcsec/pixel)")
                else:
                    # Save to a new file in the same directory
                    output_path = os.path.join(input_folder, f"{base_name}_wcs.fits")
                    fits.writeto(output_path, data, new_header, overwrite=True)
                    print(f">>> Wrote {os.path.basename(output_path)} with rough WCS (pixel scale={effective_pixel_scale:.2f} arcsec/pixel)")

        except Exception as e:
            print(f">>> Failed to process {base_name}: {e}")

    print("\n>>> All images processed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add rough WCS to FITS images using SIMBAD coordinates for SCAMP preprocessing.")
    parser.add_argument("input_folder", help="Folder containing full-field FITS images")
    parser.add_argument("--pixel-scale", type=float, default=3.92, help="Unbinned pixel scale in arcseconds per pixel (default: 3.92)")
    parser.add_argument("--no-overwrite", action="store_false", dest="overwrite",
                        help="Save new FITS files with WCS instead of overwriting originals")
    args = parser.parse_args()

    target_name = input(">>> Enter SIMBAD target name (e.g., WASP-12): ").strip()
    try:
        ra, dec = get_target_coords_from_simbad(target_name)
        process_images(args.input_folder, ra, dec, pixel_scale_arcsec=args.pixel_scale, overwrite=args.overwrite)
    except ValueError as e:
        print(e)
        exit(1)