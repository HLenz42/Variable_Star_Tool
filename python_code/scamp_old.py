import glob
import subprocess
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys


#-----------------Set up the directories-----------------
sci_im_dir = sys.argv[1]          # where are the science images?
se_output_dir = sys.argv[2]       # where are the SE catalogs?
scamp_config_dir = sys.argv[3]    # path to SCAMP config
scamp_output_dir = sys.argv[4]    # path to SCAMP output

scamp_command = "scamp"
scamp_refcat_dir = os.path.join(scamp_output_dir, "refcat")
os.makedirs(scamp_refcat_dir, exist_ok=True)

#-----------------Set up check plots-----------------
scamp_cat_type_out = "FITS_LDAC"
scamp_check_type = "DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,PHOT_ERROR"
scamp_check_file = ",".join([
    os.path.join(scamp_output_dir, name) for name in [
        "distort.pdf", "astr_interror2d.pdf", "astr_interror1d.pdf",
        "astr_referror2d.pdf", "astr_referror1d.pdf", "phot_error.pdf"
    ]
])

#-----------------SCAMP config overrides-----------------
scamp_distort_deg = "4"
scamp_astref_cat = "GAIA-DR2" # "GAIA-EDR3"
scamp_astref_band = "DEFAULT" # "1"   # for G band

#-----------------Image discovery-----------------
im_files = np.sort(glob.glob(os.path.join(sci_im_dir, "*.fit")))

#-----------------SCAMP output organization-----------------
scamp_list_file = os.path.join(scamp_config_dir, "scamp_kgmt_images.list")
scamp_cat_file_out = os.path.join(scamp_config_dir, "scamp_kgmt.ldac")

# Reset the SCAMP list file
with open(scamp_list_file, "w") as f:
    pass  # create or empty it

#-----------------Process each image-----------------
for image in im_files:
    base = os.path.splitext(os.path.basename(image))[0]
    se_cat_file = os.path.join(se_output_dir, base + "_se.ldac")
    scamp_cat_file = os.path.join(scamp_output_dir, base + "_stars.ldac")
    scamp_head_file = scamp_cat_file.replace(".ldac", ".head")
    scamp_ahead_file = scamp_cat_file.replace(".ldac", ".ahead")
    plot_dir = os.path.join(scamp_output_dir, "star_cat_plots", base)
    os.makedirs(plot_dir, exist_ok=True)

    print(f"Using SE catalog: {se_cat_file}")

    mag_faint = 10
    mag_bright = 6
    mag_zp = 20.4

    temp_se_cat = fits.open(se_cat_file, ignore_missing_end=True)
    temp_star_cat = fits.HDUList()
    
    for jj in range(len(temp_se_cat)):
        if jj == 0 or jj % 2 != 0:
            temp_star_cat.append(temp_se_cat[jj])
            temp_star_cat[jj].header = temp_se_cat[jj].header
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
                               temp_data["FLUX_RADIUS"] <= 1 * scamp_rad)
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
            temp_hdu.header = temp_se_cat[jj].header
            temp_hdu.data = temp_data
            temp_star_cat.append(temp_hdu)

    temp_star_cat.writeto(scamp_cat_file, overwrite=True)

    if os.path.exists(scamp_cat_file):
        print(f">>> Wrote LDAC: {scamp_cat_file}")
    else:
        print(f">>> Failed to write LDAC: {scamp_cat_file}")

    # Write .ahead file
    with open(scamp_ahead_file, "w") as f:
        for jj in [x for x in np.arange(2, len(temp_se_cat)) if x % 2 == 0]:
            f.write(f"CRVAL1  = {ra_deg: .8f} / RA of reference pixel\n") # RA in degrees of central pixel
            f.write(f"CRVAL2  = {dec_deg: .8f} / Dec of reference pixel\n") # Dec in degrees of central pixel
            f.write(f"CRPIX1  = {x_center: .1f} / X coordinate of reference pixel\n") # X pixel coordinate used as the refrence point
            f.write(f"CRPIX2  = {y_center: .1f} / Y coordinate of reference pixel\n") # Y pixel coordinate used as the refrence point
            f.write(f"CD1_1   = {-pixel_scale_deg:.10f} / Pixel scale matrix\n")
            f.write( "CD1_2   = 0.0\n")
            f.write( "CD2_1   = 0.0\n")
            f.write(f"CD2_2   = {pixel_scale_deg:.10f}\n")
            f.write( "CTYPE1  = 'RA---TAN'\n")
            f.write( "CTYPE2  = 'DEC--TAN'\n")
            f.write( "PHOTFLAG= 'T'\n")
            f.write( "END\n")

    # Append to SCAMP input list
    with open(scamp_list_file, "a", encoding="utf-8") as f:
    	f.write(scamp_cat_file.strip() + "\n")

#-----------------Print contents of the list file-----------------
if os.path.exists(scamp_list_file):
    print("\n>>> SCAMP input list contains:")
    with open(scamp_list_file, "r") as f:
        print(f.read())
else:
    print(">>> SCAMP list file not found.")

#-----------------Run SCAMP-----------------
command = [
    scamp_command, "@" + scamp_list_file,
    "-c", os.path.join(scamp_config_dir, "scamp.conf"),
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
    "-ASTREFMAG_LIMITS", "13.0,18.5",     # GAIA DR3 does G ~3-20     11-18.5 
	"-CROSSID_RADIUS", "1.5",             # max distance in arcsec between detected sources might increase
	"-SN_THRESHOLDS", "10.0,30.0",		  # might change to 20 50
	"-FWHM_THRESHOLDS", "2.5,7.0", 	 	  # range to be defined a star widen if necessary
	"-ELLIPTICITY_MAX", "0.3",			  # max elipticity to be a star (0=circle)
	"-POSITION_MAXERR", "0.5"]		  	  # max acceptable positional error between detected source and catalogue match


print("\n Running SCAMP with command:")
print(" ".join(command))
# run scamp
subprocess.call(command)



