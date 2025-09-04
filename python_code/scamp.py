# scamp.py

import glob
import subprocess
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm

# ----------------- Read from command line -----------------
sci_im_dir = sys.argv[1]          # Where are the science images?
se_output_dir = sys.argv[2]       # Where are the SE catalogs?
scamp_config_dir = sys.argv[3]    # Path to SCAMP config
scamp_output_dir = sys.argv[4]    # Path to SCAMP output
fits_pattern = "*.fit"  # Default to *.fit

# Verify directories
for d in [sci_im_dir, se_output_dir, scamp_config_dir, scamp_output_dir]:
    if not os.path.isdir(d):
        print(f"Error: Directory not found: {d}")
        sys.exit(1)

# Verify SCAMP config file
scamp_config_file = os.path.join(scamp_config_dir, "scamp.conf")
if not os.path.exists(scamp_config_file):
    print(f"Error: SCAMP config file not found: {scamp_config_file}")
    sys.exit(1)

# ----------------- Constants -----------------
scamp_command = "scamp"
scamp_refcat_dir = os.path.join(scamp_output_dir, "refcat")
os.makedirs(scamp_refcat_dir, exist_ok=True)

scamp_cat_type_out = "FITS_LDAC"
scamp_check_type = "DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,PHOT_ERROR"
scamp_check_file = ",".join([
    os.path.join(scamp_output_dir, name) for name in [
        "distort.pdf", "astr_interror2d.pdf", "astr_interror1d.pdf",
        "astr_referror2d.pdf", "astr_referror1d.pdf", "phot_error.pdf"
    ]
])

scamp_distort_deg = "4"
scamp_astref_cat = "GAIA-DR2"
scamp_astref_band = "DEFAULT"

pixel_scale_arcsec = 3.92  # Match SExtractor and telescope specs
pixel_scale_deg = pixel_scale_arcsec / 3600.0

# ----------------- Image list -----------------
im_files = np.sort(glob.glob(os.path.join(sci_im_dir, fits_pattern)))
if len(im_files) == 0:
    print(f"Error: No FITS files found in {sci_im_dir} with pattern {fits_pattern}")
    sys.exit(1)

scamp_list_file = os.path.join(scamp_config_dir, "scamp_kgmt_images.list")
scamp_cat_file_out = os.path.join(scamp_output_dir, "scamp_kgmt.ldac")

# Reset SCAMP input list
with open(scamp_list_file, "w") as f:
    pass

# ----------------- Process each image -----------------
for image in tqdm(im_files, desc="Processing images"):
    base = os.path.splitext(os.path.basename(image))[0]
    se_cat_file = os.path.join(se_output_dir, base + "_se.ldac")
    scamp_cat_file = os.path.join(scamp_output_dir, base + "_stars.ldac")
    scamp_ahead_file = os.path.join(scamp_output_dir, base + ".ahead")
    plot_dir = os.path.join(scamp_output_dir, "star_cat_plots", base)
    os.makedirs(plot_dir, exist_ok=True)

    if not os.path.exists(se_cat_file):
        print(f">>> Missing SE catalog: {se_cat_file}, skipping...")
        continue

    print(f">>> Using SE catalog: {se_cat_file}")

    # Get image size and rough WCS from FITS header
    try:
        with fits.open(image) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            ny, nx = data.shape
            x_center = nx / 2
            y_center = ny / 2
            if 'CRVAL1' in header and 'CRVAL2' in header:
                ra_deg = header['CRVAL1']
                dec_deg = header['CRVAL2']
                print(f">>> Using FITS header WCS: RA = {ra_deg:.6f}, Dec = {dec_deg:.6f}")
            else:
                print(f">>> Error: No WCS (CRVAL1, CRVAL2) found in {image}, skipping...")
                continue
    except Exception as e:
        print(f">>> Error reading FITS file {image}: {e}")
        continue

    mag_faint = 10
    mag_bright = 6

    try:
        temp_se_cat = fits.open(se_cat_file, ignore_missing_end=True)
        temp_star_cat = fits.HDUList()

        for jj in range(len(temp_se_cat)):
            if jj == 0 or jj % 2 != 0:
                temp_star_cat.append(temp_se_cat[jj])
            else:
                temp_data = temp_se_cat[jj].data[(temp_se_cat[jj].data["FLAGS"] <= 1)]

                plt.figure()
                plt.plot(temp_data["FLUX_RADIUS"], temp_data["MAG_AUTO"], "ko", ms=2)

                temp_data = temp_data[
                    np.logical_and(temp_data["MAG_AUTO"] >= mag_bright,
                                   temp_data["MAG_AUTO"] <= mag_faint)
                ]
                scamp_rad = np.median(temp_data["FLUX_RADIUS"])
                temp_data = temp_data[
                    np.logical_and(temp_data["FLUX_RADIUS"] >= 0.5 * scamp_rad,
                                   temp_data["FLUX_RADIUS"] <= 1.0 * scamp_rad)
                ]
                temp_data = temp_data[(temp_data["ELONGATION"] >= 0.7)]

                plt.plot(temp_data["FLUX_RADIUS"], temp_data["MAG_AUTO"], "ro", ms=2)
                plt.plot([0.0, 11], [mag_faint, mag_faint], "b:")
                plt.plot([0.0, 11], [mag_bright, mag_bright], "b:")
                plt.xlabel("FLUX_RADIUS")
                plt.ylabel("MAG_AUTO")
                plt.xlim(0.0, 9)
                plt.ylim(14., 2.)
                plt.savefig(os.path.join(plot_dir, f"stars_chip{int(jj/2)}.png"), bbox_inches="tight")
                plt.close()

                temp_hdu = fits.BinTableHDU.from_columns(temp_se_cat[jj].columns, nrows=len(temp_data))
                temp_hdu.data = temp_data
                temp_hdu.header = temp_se_cat[jj].header
                temp_star_cat.append(temp_hdu)

        temp_star_cat.writeto(scamp_cat_file, overwrite=True)
        abs_scamp_cat_file = os.path.abspath(scamp_cat_file)

        if os.path.exists(scamp_cat_file):
            print(f">>> Wrote LDAC: {scamp_cat_file}")
            with open(scamp_list_file, "a", encoding="utf-8") as f:
                f.write(abs_scamp_cat_file + "\n")
        else:
            print(f">>> Failed to write LDAC: {scamp_cat_file}")
            continue

        # ----------------- Write initial .ahead file -----------------
        with open(scamp_ahead_file, "w") as f:
            f.write(f"CRVAL1  = {ra_deg:.8f} / RA of reference pixel\n")
            f.write(f"CRVAL2  = {dec_deg:.8f} / Dec of reference pixel\n")
            f.write(f"CRPIX1  = {x_center:.1f} / X coordinate of reference pixel\n")
            f.write(f"CRPIX2  = {y_center:.1f} / Y coordinate of reference pixel\n")
            f.write(f"CD1_1   = {-pixel_scale_deg:.10f} / Pixel scale matrix\n")
            f.write("CD1_2   = 0.0\n")
            f.write("CD2_1   = 0.0\n")
            f.write(f"CD2_2   = {pixel_scale_deg:.10f}\n")
            f.write("CTYPE1  = 'RA---TAN'\n")
            f.write("CTYPE2  = 'DEC--TAN'\n")
            f.write("PHOTFLAG= 'T'\n")
            f.write("END\n")

    except Exception as e:
        print(f">>> Error processing {se_cat_file}: {e}")
        continue

# ----------------- Print SCAMP list file -----------------
print("\n>>> SCAMP input list contains:")
with open(scamp_list_file, "r") as f:
    print(f.read())

# ----------------- Run SCAMP -----------------
command = [
    scamp_command, "@" + scamp_list_file,
    "-c", scamp_config_file,
    "-MERGEDOUTCAT_TYPE", scamp_cat_type_out,
    "-MERGEDOUTCAT_NAME", scamp_cat_file_out,
    "-MATCH", "Y",
    "-SAVE_REFCATALOG", "Y",
    "-REFOUT_CATPATH", scamp_refcat_dir,
    "-CHECKPLOT_TYPE", scamp_check_type,
    "-CHECKPLOT_NAME", scamp_check_file,
    "-ASTREF_CATALOG", scamp_astref_cat,
    "-ASTREF_BAND", scamp_astref_band,
    "-DISTORT_DEGREES", scamp_distort_deg,
    "-ASTREFMAG_LIMITS", "12,19",   # mag range in the refrence cataloge's band for selecting refrence stars
    "-CROSSID_RADIUS", "3",         # sets the initial radius in arcsec for cross-matching detected sources
    "-SN_THRESHOLDS", "6,30",       # specifies the snr ratio thresholds for selecting sources for comparing to gaia
    "-FWHM_THRESHOLDS", "1,6.0",    # specifies the FWHM range used to compair sources to gaia
    "-ELLIPTICITY_MAX", "0.3",      # Max ellipticity of a source to be used when compairing to gaia (circle=0)
    "-POSITION_MAXERR", "0.6"       # max position uncertainty in arcmin for matching to catalogue
]

print("\n>>> Running SCAMP with command:")
print(" ".join(command))
subprocess.call(command)

# ----------------- Generate refined .ahead files from .head files -----------------
for image in tqdm(im_files, desc="Generating refined .ahead files"):
    base = os.path.splitext(os.path.basename(image))[0]
    head_file = os.path.join(scamp_output_dir, base + "_stars.head")
    scamp_ahead_file = os.path.join(scamp_output_dir, base + ".ahead")
    
    if os.path.exists(head_file):
        # Read the .head file as a FITS header
        with open(head_file, "r") as f:
            head_content = f.read()
        head_lines = head_content.split("\n")
        new_header = fits.Header()
        for line in head_lines:
            if line.strip() and not line.startswith("COMMENT") and not line.startswith("HISTORY"):
                try:
                    key, value = line.split("=", 1)
                    key = key.strip()
                    value = value.split("/")[0].strip()
                    new_header[key] = value
                except:
                    continue
        
        # Write a new .ahead file with core WCS keywords
        with open(scamp_ahead_file, "w") as f:
            f.write(f"CRVAL1  = {float(new_header['CRVAL1']):.8f} / RA of reference pixel\n")
            f.write(f"CRVAL2  = {float(new_header['CRVAL2']):.8f} / Dec of reference pixel\n")
            f.write(f"CRPIX1  = {float(new_header['CRPIX1']):.1f} / X coordinate of reference pixel\n")
            f.write(f"CRPIX2  = {float(new_header['CRPIX2']):.1f} / Y coordinate of reference pixel\n")
            f.write(f"CD1_1   = {float(new_header['CD1_1']):.10f} / Pixel scale matrix\n")
            f.write(f"CD1_2   = {float(new_header['CD1_2']):.10f}\n")
            f.write(f"CD2_1   = {float(new_header['CD2_1']):.10f}\n")
            f.write(f"CD2_2   = {float(new_header['CD2_2']):.10f}\n")
            f.write(f"CTYPE1  = '{new_header['CTYPE1']}'\n")
            f.write(f"CTYPE2  = '{new_header['CTYPE2']}'\n")
            f.write(f"PHOTFLAG= 'T'\n")
            f.write("END\n")
        print(f">>> Wrote refined .ahead file: {scamp_ahead_file}")
    else:
        print(f">>> No .head file found for {image}, skipping .ahead generation")



# # scamp.py

# import glob
# import subprocess
# import os
# import sys
# import numpy as np
# from astropy.io import fits
# import matplotlib.pyplot as plt
# from tqdm import tqdm

# # ----------------- Read from command line -----------------

# sci_im_dir = sys.argv[1]          # Where are the science images?
# se_output_dir = sys.argv[2]       # Where are the SE catalogs?
# scamp_config_dir = sys.argv[3]    # Path to SCAMP config
# scamp_output_dir = sys.argv[4]    # Path to SCAMP output
# fits_pattern = "*.fit"  # Default to *.fits

# # Verify directories
# for d in [sci_im_dir, se_output_dir, scamp_config_dir, scamp_output_dir]:
#     if not os.path.isdir(d):
#         print(f"Error: Directory not found: {d}")
#         sys.exit(1)

# # Verify SCAMP config file
# scamp_config_file = os.path.join(scamp_config_dir, "scamp.conf")
# if not os.path.exists(scamp_config_file):
#     print(f"Error: SCAMP config file not found: {scamp_config_file}")
#     sys.exit(1)

# # ----------------- Constants -----------------
# scamp_command = "scamp"
# scamp_refcat_dir = os.path.join(scamp_output_dir, "refcat")
# os.makedirs(scamp_refcat_dir, exist_ok=True)

# scamp_cat_type_out = "FITS_LDAC"
# scamp_check_type = "DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,PHOT_ERROR"
# scamp_check_file = ",".join([
#     os.path.join(scamp_output_dir, name) for name in [
#         "distort.pdf", "astr_interror2d.pdf", "astr_interror1d.pdf",
#         "astr_referror2d.pdf", "astr_referror1d.pdf", "phot_error.pdf"
#     ]
# ])

# scamp_distort_deg = "4"
# scamp_astref_cat = "GAIA-DR2"
# scamp_astref_band = "DEFAULT"

# pixel_scale_arcsec = 3.92  # Match SExtractor and telescope specs
# pixel_scale_deg = pixel_scale_arcsec / 3600.0

# # ----------------- Image list -----------------
# im_files = np.sort(glob.glob(os.path.join(sci_im_dir, fits_pattern)))
# if len(im_files) == 0:  # Fixed: Check array length instead of truth value
#     print(f"Error: No FITS files found in {sci_im_dir} with pattern {fits_pattern}")
#     sys.exit(1)

# scamp_list_file = os.path.join(scamp_config_dir, "scamp_kgmt_images.list")
# scamp_cat_file_out = os.path.join(scamp_output_dir, "scamp_kgmt.ldac")

# # Reset SCAMP input list
# with open(scamp_list_file, "w") as f:
#     pass

# # ----------------- Process each image -----------------
# for image in tqdm(im_files, desc="Processing images"):
#     base = os.path.splitext(os.path.basename(image))[0]
#     se_cat_file = os.path.join(se_output_dir, base + "_se.ldac")
#     scamp_cat_file = os.path.join(scamp_output_dir, base + "_stars.ldac")
#     scamp_ahead_file = scamp_cat_file.replace(".ldac", ".ahead")
#     plot_dir = os.path.join(scamp_output_dir, "star_cat_plots", base)
#     os.makedirs(plot_dir, exist_ok=True)

#     if not os.path.exists(se_cat_file):
#         print(f">>> Missing SE catalog: {se_cat_file}, skipping...")
#         continue

#     print(f">>> Using SE catalog: {se_cat_file}")

#     # Get image size and rough WCS from FITS header
#     try:
#         with fits.open(image) as hdul:
#             data = hdul[0].data
#             header = hdul[0].header
#             ny, nx = data.shape
#             x_center = nx / 2
#             y_center = ny / 2
#             if 'CRVAL1' in header and 'CRVAL2' in header:
#                 ra_deg = header['CRVAL1']
#                 dec_deg = header['CRVAL2']
#                 print(f">>> Using FITS header WCS: RA = {ra_deg:.6f}, Dec = {dec_deg:.6f}")
#             else:
#                 print(f">>> Error: No WCS (CRVAL1, CRVAL2) found in {image}, skipping...")
#                 continue
#     except Exception as e:
#         print(f">>> Error reading FITS file {image}: {e}")
#         continue

#     mag_faint = 10
#     mag_bright = 6

#     try:
#         temp_se_cat = fits.open(se_cat_file, ignore_missing_end=True)
#         temp_star_cat = fits.HDUList()

#         for jj in range(len(temp_se_cat)):
#             if jj == 0 or jj % 2 != 0:
#                 temp_star_cat.append(temp_se_cat[jj])
#             else:
#                 temp_data = temp_se_cat[jj].data[(temp_se_cat[jj].data["FLAGS"] <= 1)]

#                 plt.figure()
#                 plt.plot(temp_data["FLUX_RADIUS"], temp_data["MAG_AUTO"], "ko", ms=2)

#                 temp_data = temp_data[
#                     np.logical_and(temp_data["MAG_AUTO"] >= mag_bright,
#                                    temp_data["MAG_AUTO"] <= mag_faint)
#                 ]
#                 scamp_rad = np.median(temp_data["FLUX_RADIUS"])
#                 temp_data = temp_data[
#                     np.logical_and(temp_data["FLUX_RADIUS"] >= 0.5 * scamp_rad,
#                                    temp_data["FLUX_RADIUS"] <= 1.0 * scamp_rad)
#                 ]
#                 temp_data = temp_data[(temp_data["ELONGATION"] >= 0.7)]

#                 plt.plot(temp_data["FLUX_RADIUS"], temp_data["MAG_AUTO"], "ro", ms=2)
#                 plt.plot([0.0, 11], [mag_faint, mag_faint], "b:")
#                 plt.plot([0.0, 11], [mag_bright, mag_bright], "b:")
#                 plt.xlabel("FLUX_RADIUS")
#                 plt.ylabel("MAG_AUTO")
#                 plt.xlim(0.0, 9)
#                 plt.ylim(14., 2.)
#                 plt.savefig(os.path.join(plot_dir, f"stars_chip{int(jj/2)}.png"), bbox_inches="tight")
#                 plt.close()

#                 temp_hdu = fits.BinTableHDU.from_columns(temp_se_cat[jj].columns, nrows=len(temp_data))
#                 temp_hdu.data = temp_data
#                 temp_hdu.header = temp_se_cat[jj].header
#                 temp_star_cat.append(temp_hdu)

#         temp_star_cat.writeto(scamp_cat_file, overwrite=True)
#         abs_scamp_cat_file = os.path.abspath(scamp_cat_file)

#         if os.path.exists(scamp_cat_file):
#             print(f">>> Wrote LDAC: {scamp_cat_file}")
#             with open(scamp_list_file, "a", encoding="utf-8") as f:
#                 f.write(abs_scamp_cat_file + "\n")
#         else:
#             print(f">>> Failed to write LDAC: {scamp_cat_file}")
#             continue

#         # ----------------- Write .ahead file -----------------
#         with open(scamp_ahead_file, "w") as f:
#             f.write(f"CRVAL1  = {ra_deg:.8f} / RA of reference pixel\n")
#             f.write(f"CRVAL2  = {dec_deg:.8f} / Dec of reference pixel\n")
#             f.write(f"CRPIX1  = {x_center:.1f} / X coordinate of reference pixel\n")
#             f.write(f"CRPIX2  = {y_center:.1f} / Y coordinate of reference pixel\n")
#             f.write(f"CD1_1   = {-pixel_scale_deg:.10f} / Pixel scale matrix\n")
#             f.write("CD1_2   = 0.0\n")
#             f.write("CD2_1   = 0.0\n")
#             f.write(f"CD2_2   = {pixel_scale_deg:.10f}\n")
#             f.write("CTYPE1  = 'RA---TAN'\n")
#             f.write("CTYPE2  = 'DEC--TAN'\n")
#             f.write("PHOTFLAG= 'T'\n")
#             f.write("END\n")

#     except Exception as e:
#         print(f">>> Error processing {se_cat_file}: {e}")
#         continue

# # ----------------- Print SCAMP list file -----------------
# print("\n>>> SCAMP input list contains:")
# with open(scamp_list_file, "r") as f:
#     print(f.read())

# # ----------------- Run SCAMP -----------------
# command = [
#     scamp_command, "@" + scamp_list_file,
#     "-c", scamp_config_file,
#     "-MERGEDOUTCAT_TYPE", scamp_cat_type_out,
#     "-MERGEDOUTCAT_NAME", scamp_cat_file_out,
#     "-MATCH", "Y",
#     "-SAVE_REFCATALOG", "Y",
#     "-REFOUT_CATPATH", scamp_refcat_dir,
#     "-CHECKPLOT_TYPE", scamp_check_type,
#     "-CHECKPLOT_NAME", scamp_check_file,
#     "-ASTREF_CATALOG", scamp_astref_cat,
#     "-ASTREF_BAND", scamp_astref_band,
#     "-DISTORT_DEGREES", scamp_distort_deg,
#     "-ASTREFMAG_LIMITS", "13.0,18.5",
#     "-CROSSID_RADIUS", "1.5",
#     "-SN_THRESHOLDS", "10.0,30.0",
#     "-FWHM_THRESHOLDS", "2.5,7.0",
#     "-ELLIPTICITY_MAX", "0.3",
#     "-POSITION_MAXERR", "0.5"
# ]

# print("\n>>> Running SCAMP with command:")
# print(" ".join(command))
# subprocess.call(command)