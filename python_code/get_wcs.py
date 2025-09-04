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

    # Quadrant boundaries for 4x4 split
    dx = nx // 4
    dy = ny // 4
    quadrants = {}
    q_num = 0
    for i in range(4):
        for j in range(4):
            y1 = i * dy
            y2 = (i + 1) * dy if i < 3 else ny
            x1 = j * dx
            x2 = (j + 1) * dx if j < 3 else nx
            quadrants[q_num] = (y1, y2, x1, x2)
            q_num += 1
            print(f"     Defined quadrant {q_num-1}: y1={y1}, y2={y2}, x1={x1}, x2={x2}")

    for q_num, (y1, y2, x1, x2) in quadrants.items():
        print(f"  -> Solving WCS for quadrant {q_num} of {base_name}")
        quad_data = data[y1:y2, x1:x2]
        quad_header = header.copy()

        # Temporary cropped image for solving
        quad_temp_path = os.path.join(output_dir, f"{base_name}_q{q_num}.fit")
        try:
            fits.writeto(quad_temp_path, quad_data, quad_header, overwrite=True)

            # Solve WCS
            try:
                wcs_header = ast.solve_from_image(
                    quad_temp_path,
                    force_image_upload=True,
                    solve_timeout=solve_timeout,
                    detect_threshold=3
                )

                # Save WCS-enhanced image
                wcs_path = os.path.join(output_dir, f"{base_name}_q{q_num}_wcs.fit")
                fits.writeto(wcs_path, quad_data, wcs_header, overwrite=True)
                print(f"     WCS solution saved as: {os.path.basename(wcs_path)}")

            except TimeoutError:
                print(f"     WCS solve timed out for quadrant {q_num} after {solve_timeout} seconds")
                # Save unsolved quadrant to subfolder
                unsolved_path = os.path.join(unsolved_dir, f"{base_name}_q{q_num}_unsolved.fit")
                fits.writeto(unsolved_path, quad_data, quad_header, overwrite=True)
                print(f"     Saved unsolved quadrant to: {os.path.basename(unsolved_path)}")

            except Exception as e:
                print(f"     WCS solve failed for quadrant {q_num}: {e}")
                # Save unsolved quadrant to subfolder
                unsolved_path = os.path.join(unsolved_dir, f"{base_name}_q{q_num}_unsolved.fit")
                fits.writeto(unsolved_path, quad_data, quad_header, overwrite=True)
                print(f"     Saved unsolved quadrant to: {os.path.basename(unsolved_path)}")

            finally:
                # Clean up temporary file
                if os.path.exists(quad_temp_path):
                    os.remove(quad_temp_path)
                    print(f"     Removed temporary file: {os.path.basename(quad_temp_path)}")

        except Exception as e:
            print(f"     Failed to process quadrant {q_num}: {e}")

# import os
# import glob
# import sys
# from astropy.io import fits
# from astroquery.astrometry_net import AstrometryNet



# # ------------------ Configuration ------------------
# API_KEY = 'rjkkkkxfpcmklonl'  # Dr. L's key

# # Check for command-line arguments
# if len(sys.argv) != 3:
#     print("Usage: python3 get_wcs.py <input_dir> <output_dir>")
#     sys.exit(1)

# input_dir = sys.argv[1]
# output_dir = sys.argv[2]

# # fits.open(input_dir + "/HATP-P-3b_Apr17_25s_8C-001_o_p.fit").info()

# print(f"input_dir: {input_dir}")
# print(f"output_dir: {output_dir}")

# # Validate directories
# if not os.path.isdir(input_dir):
#     print(f"Error: Input directory '{input_dir}' does not exist.")
#     sys.exit(1)

# os.makedirs(output_dir, exist_ok=True)
# # ---------------------------------------------------

# # Set up Astrometry.net
# ast = AstrometryNet()
# ast.api_key = API_KEY

# # Gather FITS images
# # fits_files = sorted(glob.glob(os.path.join(input_dir, '*.fit') + 
# #                               os.path.join(input_dir, '*.fits') +
# #                               os.path.join(input_dir, '*.FITS')))

# fits_files = sorted(
#     glob.glob(os.path.join(input_dir, '*.[fF][iI][tT][sS]')) +
#     glob.glob(os.path.join(input_dir, '*.[fF][iI][tT]'))
# )
# print(f"Found {len(fits_files)} FITS files in {input_dir}: {fits_files}")
# if not fits_files:
#     print(f"Error: No FITS files found in {input_dir}. Check previous pipeline steps.")
#     sys.exit(1)

# print(fits_files)

# for fits_path in fits_files:
#     base_name = os.path.splitext(os.path.basename(fits_path))[0]
#     print(f"\nProcessing image: {base_name}")

#     # Read image data and header
#     try:
#         data, header = fits.getdata(fits_path, header=True)
#         ny, nx = data.shape
#     except Exception as e:
#         print(f"     Failed to read FITS file {fits_path}: {e}")
#         continue

#     # Quadrant boundaries for 4x4 split
#     dx = nx // 4
#     dy = ny // 4
#     quadrants = {}
#     q_num = 0
#     for i in range(4):
#         for j in range(4):
#             y1 = i * dy
#             y2 = (i + 1) * dy if i < 3 else ny
#             x1 = j * dx
#             x2 = (j + 1) * dx if j < 3 else nx
#             quadrants[q_num] = (y1, y2, x1, x2)
#             q_num += 1
#             print(">>> Split quadrants")

#     for q_num, (y1, y2, x1, x2) in quadrants.items():
#         print(f"  -> Solving WCS for quadrant {q_num} of {base_name}")
#         quad_data = data[y1:y2, x1:x2]
#         quad_header = header.copy()

#         # Temporary cropped image for solving
#         quad_temp_path = os.path.join(output_dir, f"{base_name}_q{q_num}.fits")
#         try:
#             fits.writeto(quad_temp_path, quad_data, quad_header, overwrite=True)

#             # Solve WCS
#             try:
#                 wcs_header = ast.solve_from_image(
#                     quad_temp_path,
#                     force_image_upload=True,
#                     solve_timeout=300,
#                     detect_threshold=3
#                 )

#                 # Save WCS-enhanced image
#                 wcs_path = os.path.join(output_dir, f"{base_name}_q{q_num}_wcs.fits")
#                 fits.writeto(wcs_path, quad_data, wcs_header, overwrite=True)
#                 print(f"     WCS solution saved as: {os.path.basename(wcs_path)}")

#             except Exception as e:
#                 print(f"     WCS solve failed for quadrant {q_num}: {e}")

#             finally:
#                 # Clean up temporary file
#                 if os.path.exists(quad_temp_path):
#                     os.remove(quad_temp_path)
#                     print(f"     Removed temporary file: {os.path.basename(quad_temp_path)}")

#         except Exception as e:
#             print(f"     Failed to process quadrant {q_num}: {e}")