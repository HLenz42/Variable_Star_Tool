import os
import glob
import subprocess
import numpy as np
import sys
from astropy.io import fits


#-----------------Set the directories-----------------
se_config_dir = sys.argv[2]
se_output_dir = sys.argv[3]
config_name = sys.argv[4]
checkplot = sys.argv[5]
path_2_background_imgs = sys.argv[6]

if checkplot == "y":
    CHECKIMAGE_TYPE = "BACKGROUND" # "APERTURES" # "-BACKGROUND"

elif checkplot == "n":
    CHECKIMAGE_TYPE = "APERTURES" # "NONE"

se_config_file = se_config_dir + config_name
se_param_file = se_config_dir + "/kgmt.param"

sci_im_dir =  sys.argv[1]  # where are the science images?
# wt_im_dir = prefix + "/data/weightmap_k_split"  # optional: weight images

se_command = "source-extractor"

#-----------------Set up KGMT specific parameters-----------------           
im_gain = "1.42"  # e/ADU
im_rdnoise = "14.6"  # electrons
im_pixscale = "3.92"  # arcsec/pixel
mag_zp = 20.4  # ADU/Sec
extinction_coeff = 0.058


#-------------------Load the image & weightmap filenames into arrays-------------------
#Grab all the images that need to be source extracted
im_files = np.sort(glob.glob(sci_im_dir+"/*.fit") + 
                   glob.glob(sci_im_dir+"/*.fits") +
                   glob.glob(sci_im_dir+"/*.FITS"))
# wt_files = np.sort(glob.glob(wt_im_dir+"/*.fits"))

#-------------------Now run source extractor on all the images!-------------------
for ii in range(len(im_files)):
    image = im_files[ii]
    # weight = wt_files[ii]
    # air_mass = fits.open(image)[0].header["AIRMASS"] 

    # se_cat_file = image.replace(sci_im_dir,se_output_dir).replace(".fits","_se.ldac") #use an output catalogue name based on the image filename
    se_cat_file = image.replace(sci_im_dir,se_output_dir).replace(".fit","_se.ldac") #use an output catalogue name based on the image filename
    if checkplot == "y":
        background_img = image.replace(sci_im_dir,se_output_dir).replace(".fit","_back.fit")

    if checkplot == "n":
        background_img = image.replace(sci_im_dir,se_output_dir).replace(".fit","_aper.fit")

    #-------------------Set up the command to run Source-Extractor-------------------
    #Set up the call to source extractor. Use the default values, except where listed below.
    command = [se_command, image,
                "-c", se_config_file,
                "-CATALOG_NAME", se_cat_file,
                "-CATALOG_TYPE", "FITS_LDAC",
                "-PARAMETERS_NAME", se_param_file,
                # "-WEIGHT_IMAGE", weight,
                "-FILTER_NAME", se_config_dir+"/default.conv",
                "-SATUR_LEVEL", "65535",
                # "-MAG_ZEROPOINT", str(mag_zp+(extinction_coeff*air_mass)),
                "-GAIN", im_gain,
                "-PIXEL_SCALE", im_pixscale,
                "-CHECKIMAGE_TYPE", CHECKIMAGE_TYPE,
                "-CHECKIMAGE_NAME", background_img]
    print(f"\n{command}\n")
    subprocess.call(command) #Run source extractor!