import os
import glob
import sys
from astropy.io import fits
from astroquery.astrometry_net import AstrometryNet
from tqdm import tqdm

# ------------------ Configuration ------------------
API_KEY = 'rjkkkkxfpcmklonl'  # Dr. L's key

input_dir = sys.argv[1]
output_dir = sys.argv[2]
solve_timeout = 600  # Default to 10 minutes

# Validate arguments
if not os.path.isdir(input_dir):
    print(f"Error: Input directory '{input_dir}' does not exist.")
    sys.exit(1)
if solve_timeout <= 0:
    print("Error: solve_timeout must be a positive integer.")
    sys.exit(1)

# Create output directories
os.makedirs(output_dir, exist_ok=True)
unsolved_dir = os.path.join(output_dir, "unsolved_quadrants")
os.makedirs(unsolved_dir, exist_ok=True)

print(f"Input directory: {input_dir}")
print(f"Output directory: {output_dir}")
print(f"Unsolved quadrants directory: {unsolved_dir}")
print(f"Solve timeout: {solve_timeout} seconds")

# Set up Astrometry.net
ast = AstrometryNet()
ast.api_key = API_KEY

# Gather FITS images
fits_files = sorted(glob.glob(os.path.join(input_dir, '*.[fF][iI][tT][sS]')) +
                    glob.glob(os.path.join(input_dir, '*.[fF][iI][tT]')))
print(f"Found {len(fits_files)} FITS files in {input_dir}")
if not fits_files:
    print(f"Error: No FITS files found in {input_dir}. Check previous pipeline steps.")
    sys.exit(1)

# Process each FITS file
for fits_path in tqdm(fits_files, desc="Processing FITS files"):
    base_name = os.path.splitext(os.path.basename(fits_path))[0]
    print(f"\nProcessing image: {base_name}")

    # Read image data and header
    try:
        data, header = fits.getdata(fits_path, header=True)
        ny, nx = data.shape
    except Exception as e:
        print(f"     Failed to read FITS file {fits_path}: {e}")
        continue

    # Treat entire image as quadrant 0
    q_num = 0
    print(f"  -> Solving WCS for full image (q{q_num}) of {base_name}")

    # Temporary image for solving (use full image)
    temp_path = os.path.join(output_dir, f"{base_name}_q{q_num}.fit")
    try:
        fits.writeto(temp_path, data, header, overwrite=True)

        # Solve WCS
        try:
            wcs_header = ast.solve_from_image(
                temp_path,
                force_image_upload=True,
                solve_timeout=solve_timeout,
                detect_threshold=3
            )

            # Save WCS-enhanced image
            wcs_path = os.path.join(output_dir, f"{base_name}_q{q_num}_wcs.fit")
            fits.writeto(wcs_path, data, wcs_header, overwrite=True)
            print(f"     WCS solution saved as: {os.path.basename(wcs_path)}")

        except TimeoutError:
            print(f"     WCS solve timed out for full image (q{q_num}) after {solve_timeout} seconds")
            # Save unsolved image to subfolder
            unsolved_path = os.path.join(unsolved_dir, f"{base_name}_q{q_num}_unsolved.fit")
            fits.writeto(unsolved_path, data, header, overwrite=True)
            print(f"     Saved unsolved image to: {os.path.basename(unsolved_path)}")

        except Exception as e:
            print(f"     WCS solve failed for full image (q{q_num}): {e}")
            # Save unsolved image to subfolder
            unsolved_path = os.path.join(unsolved_dir, f"{base_name}_q{q_num}_unsolved.fit")
            fits.writeto(unsolved_path, data, header, overwrite=True)
            print(f"     Saved unsolved image to: {os.path.basename(unsolved_path)}")

        finally:
            # Clean up temporary file
            if os.path.exists(temp_path):
                os.remove(temp_path)
                print(f"     Removed temporary file: {os.path.basename(temp_path)}")

    except Exception as e:
        print(f"     Failed to process full image (q{q_num}): {e}")