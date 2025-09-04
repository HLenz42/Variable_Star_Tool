import subprocess
import os
import glob
import pandas as pd

# cd ffmpeg-*-static
# export PATH="$PWD:$PATH"
# export PATH="$HOME/ffmpeg-*-static:$PATH"
# ffmpeg -version

target = "Wasp12b"
date = "2025Jan27"

# target = "M92"
# date = "2021Jun02"

telescope = "kgmt"
# telescope = "arct"
# telescope = "cmt"

# prefix = "/media/localadmin/KGMTTOSHIBA/variable_star_tool"
prefix = "/media/hlenz/KGMTTOSHIBA/variable_star_tool"

# Main Dir
path_2_python_code = prefix + "/python_code"
path_2_flat = prefix + f"/flats/flat_{telescope}"
path_2_se_config = prefix + f"/se_config/se_config_{telescope}"
path_2_se_wcs_config = prefix + f"/se_config/se_config_wcs_{telescope}"
path_2_scamp_config = prefix + f"/scamp_config/scamp_config_{telescope}"
path_2_swarp_config = prefix + f"/swarp_config/swarp_config_{telescope}"



# --------- Project Specific Dir -----------
# Data
path_2_images              = prefix + f"/{target}_{date}/1_data/images"
path_2_darks               = prefix + f"/{target}_{date}/1_data/darks"
path_2_bias                = prefix + f"/{target}_{date}/1_data/bias"
path_2_processed           = prefix + f"/{target}_{date}/1_data/processed"
path_2_processed_bias      = prefix + f"/{target}_{date}/1_data/master_bias"
path_2_processed_dark      = prefix + f"/{target}_{date}/1_data/master_dark"
path_2_processed_rough_wcs = prefix + f"/{target}_{date}/1_data/processed_rough_wcs"
path_2_processed_wcs       = prefix + f"/{target}_{date}/1_data/processed_wcs"
path_2_background_imgs     = prefix + f"/{target}_{date}/1_data/processed_wcs_background_sub"

# Photometry Pipline
path_2_se_output          = prefix + f"/{target}_{date}/2_photometry_pipline/se_output"
path_2_se_output_v1       = prefix + f"/{target}_{date}/2_photometry_pipline/se_output_backsub"
path_2_se_wcs_output      = prefix + f"/{target}_{date}/2_photometry_pipline/se_output_wcs"
path_2_scamp_output       = prefix + f"/{target}_{date}/2_photometry_pipline/scamp_output"
path_2_se_psf_run1_config = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/se_r1_config"
path_2_se_psf_run1_output = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/se_r1_output"
path_2_psfex_config       = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/psf_config"
path_2_psfex_output       = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/psf_output"
path_2_se_psf_run2_config = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/se_r2_config"
path_2_se_psf_run2_output = prefix + f"/{target}_{date}/2_photometry_pipline/psf_phot/se_r2_output"

# CSV Pipline
path_2_combined_star_cats      = prefix + f"/{target}_{date}/3_csv_pipline/combined_star_cats"
path_2_combined_star_only_cats = prefix + f"/{target}_{date}/3_csv_pipline/combined_star_filtered_cats"
path_2_star_tracked_cat        = prefix + f"/{target}_{date}/3_csv_pipline/star_tracked_cats"
path_2_star_tracked_simbad_cat = prefix + f"/{target}_{date}/3_csv_pipline/star_tracked_simbad_cats"
path_2_simbad_star_class_plots = prefix + f"/{target}_{date}/3_csv_pipline/star_tracked_simbad_cats/simbad_star_class_plots"

# Diff Phot Pipline
path_2_light_curves_wo_diffphot = prefix + f"/{target}_{date}/4_diff_phot_pipline/light_curves_before_diffphot"
path_2_csv_variables_defined    = prefix + f"/{target}_{date}/4_diff_phot_pipline/identify_variables"
path_2_csv_diff_phot            = prefix + f"/{target}_{date}/4_diff_phot_pipline"
path_2_light_curves_w_diffphot  = prefix + f"/{target}_{date}/4_diff_phot_pipline/light_curves_after_diffphot"
path_2_indiv_star_study         = prefix + f"/{target}_{date}/4_diff_phot_pipline/individual_star_index"

os.makedirs(path_2_python_code, exist_ok=True)
os.makedirs(path_2_images, exist_ok=True)
os.makedirs(path_2_swarp_config, exist_ok=True)
os.makedirs(path_2_darks, exist_ok=True)
os.makedirs(path_2_bias, exist_ok=True)
os.makedirs(path_2_processed, exist_ok=True)
os.makedirs(path_2_processed_bias, exist_ok=True)
os.makedirs(path_2_processed_dark, exist_ok=True)
os.makedirs(path_2_processed_wcs, exist_ok=True)
# os.makedirs(path_2_processed_rough_wcs, exist_ok=True)
os.makedirs(path_2_se_output, exist_ok=True)
os.makedirs(path_2_scamp_output, exist_ok=True)
os.makedirs(path_2_se_wcs_config, exist_ok=True)
os.makedirs(path_2_se_wcs_output, exist_ok=True)
os.makedirs(path_2_combined_star_cats, exist_ok=True)
os.makedirs(path_2_combined_star_only_cats, exist_ok=True)
os.makedirs(path_2_star_tracked_cat, exist_ok=True)
os.makedirs(path_2_star_tracked_simbad_cat, exist_ok=True)
os.makedirs(path_2_simbad_star_class_plots, exist_ok=True)
os.makedirs(path_2_light_curves_wo_diffphot, exist_ok=True)
os.makedirs(path_2_csv_variables_defined, exist_ok=True)
os.makedirs(path_2_csv_diff_phot, exist_ok=True)
os.makedirs(path_2_light_curves_w_diffphot, exist_ok=True)
os.makedirs(path_2_se_psf_run1_config, exist_ok=True)
os.makedirs(path_2_se_psf_run1_output, exist_ok=True)
os.makedirs(path_2_psfex_config, exist_ok=True)
os.makedirs(path_2_psfex_output, exist_ok=True)
os.makedirs(path_2_se_psf_run2_config, exist_ok=True)
os.makedirs(path_2_se_psf_run2_output, exist_ok=True)
os.makedirs(path_2_indiv_star_study, exist_ok=True)
# os.makedirs(, exist_ok=True)


print()
print("======================= RAO Variable Star Detection Tool =======================" )

print("""
    *   *            *       *     *        *    *   *   *            *       *      
      *         *         *    *         *            *     *      *     *      
 *           *        *      *      *        *    *             *              *
      *           *        *     *       //              *              *
             *                          //            *
                               ___o |==// 
                              /\  \/  //|\    
                             / /        | \    
                             ` `        '  '  

""")
print("============== Written by Heather Lenz (heather.lenz@ucalgary.ca) ==============" +"\n")

print(">>> Program Overview\n")
print(">>> Section 1: Image Pre-Processing")
print("(1)  Check prefix")
print("(2)  Create master bias")
print("(3)  Create master dark")
print("(4)  Apply dark, bias, flat to images\n")
print(">>> Section 2: Photometry & Star Filtering")
print("(5)  Get images WCS (astronomy.net/Scamp)")
print("(6)  Do photometry on images using Source-Extractor or PSFEx")
print("(7)  Filter stars out of Source-Extractor catalogues")
print("(8)  Track stars though images")
print("(9)  Create star tracking check plots")
print("(10)  Plot a star timelaps to check tracking")
print("(11) Get star classifications from SIMBAD")
print("(12) Plot SIMBAD star classifications over images\n")
print(">>> Section 3: Variable Star Analysis")
print("(13) Plot light curves before differential photometry")
print("(14) Perform differential photometry")
print("     (a) Check for variable stars defined by SIMBAD & Generate analysis log")
print("     (b) Specify known variable and non-variable star indices")
print("     (c) 1st run of differential photometry")
print("     (d) Find new variables based on average variability")
print("     (e) 2nd run of differential photometry")
print("(15) Plot light curves after differential photometry")
print("(16) Plot star profiles for specific stars")
print("(17) Plot star timelaps")
print("(18) Save one star csv")
print("(19) Fit exoplanet curve!!!")
print()
print("(0) Print manual\n")

# Initial Question to select which step you would like to start from

# Initialize skip flags
skip_2_step1 = False
skip_2_step2 = False
skip_2_step3 = False
skip_2_step4 = False
skip_2_step5 = False
skip_2_step6 = False
skip_2_step7 = False
skip_2_step8 = False
skip_2_step9 = False
skip_2_step10 = False
skip_2_step11 = False
skip_2_step12 = False
skip_2_step13 = False
skip_2_step14 = False
skip_2_step15 = False
skip_2_step16 = False
skip_2_step17 = False
skip_2_step18 = False
skip_2_step19 = False
skip_2_step20 = False


# Prompt user to select a step
while True:
    print("================================================================================")
    print(">>> Which step would you like to start at? (NOTE: If this is your first run, you must start at 1)")
    
    try:
        choice1 = int(input("Answer: "))
        
        if choice1 == 1:
            print("\n>>> Going to step 1 ...\n")
            skip_2_step1 = True
            break

        elif choice1 == 2:
            print("\n>>> Going to step 2 ...\n")
            skip_2_step2 = True
            break

        elif choice1 == 3:
            print("\n>>> Going to step 3 ...\n")
            skip_2_step3 = True
            break

        elif choice1 == 4:
            print("\n>>> Going to step 4 ...\n")
            skip_2_step4 = True
            break

        elif choice1 == 5:
            print("\n>>> Going to step 5 ...\n")
            skip_2_step5 = True
            break

        elif choice1 == 6:
            print("\n>>> Going to step 6 ...\n")
            skip_2_step6 = True
            break

        elif choice1 == 7:
            print("\n>>> Going to step 7 ...\n")
            skip_2_step7 = True
            break

        elif choice1 == 8:
            print("\n>>> Going to step 8 ...\n")
            skip_2_step8 = True
            break

        elif choice1 == 9:
            print("\n>>> Going to step 9 ...\n")
            skip_2_step9 = True
            break

        elif choice1 == 10:
            print("\n>>> Going to step 10 ...\n")
            skip_2_step10 = True
            break

        elif choice1 == 11:
            print("\n>>> Going to step 11 ...\n")
            skip_2_step11 = True
            break

        elif choice1 == 12:
            print("\n>>> Going to step 12 ...\n")
            skip_2_step12 = True
            break

        elif choice1 == 13:
            print("\n>>> Going to step 13 ...\n")
            skip_2_step13 = True
            break

        elif choice1 == 14:
            print("\n>>> Going to step 14 ...\n")
            skip_2_step14 = True
            break

        elif choice1 == 15:
            print("\n>>> Going to step 15 ...\n")
            skip_2_step15 = True
            break

        elif choice1 == 16:
            print("\n>>> Going to step 16 ...\n")
            skip_2_step16 = True
            break

        elif choice1 == 17:
            print("\n>>> Going to step 17 ...\n")
            skip_2_step17 = True
            break               

        elif choice1 == 18:
            print("\n>>> Going to step 18 ...\n")
            skip_2_step18 = True
            break

        elif choice1 == 19:
            print("\n>>> Going to step 19 ...\n")
            skip_2_step19 = True
            break

        elif choice1 == 0:
            print("\n>>> Printing 'variable_star_tool_manual.txt' ...\n")
            print("================================================================================")
            try:
                with open("/home/hlenz/Projects/variable_star_tool/variable_star_tool_manual.txt", "r") as file:
                    print(file.read())
                    print()
            except FileNotFoundError:
                print(">>> Error: Manual file not found. Please check the file path.\n")
            continue  # Return to prompt after printing manual
        else:
            print("\n>>> Error: Please select a valid step (0-8)\n")

    except ValueError:
        print("\n>>> Error: Please enter a valid integer (0-8)\n")
		

# Ensure the prefex is correct 
if skip_2_step1 == True:
	while True:
		print("================================================================================" +"\n")	

		print(f">>> Is this prefix correct?          =  {prefix}")
		print(f">>> Is this the current target?      =  {target}")
		print(f">>> Is this the date for the target? =  {date}")
		print("(1) Yes")
		print("(2) No" + "\n")

		# from astropy.io import fits

		# with fits.open(path_2_processed + "/Wasp12b_20-001b_bdf.fit") as hdul:
		#     print(f"Number of HDUs: {len(hdul)}")
		#     for i, hdu in enumerate(hdul):
		#         print(f"HDU {i}: {type(hdu)}, shape={getattr(hdu.data, 'shape', None)}")

		# with fits.open(path_2_images + "/Wasp12b_20-001b.fit") as hdul:
		#     print(f"Number of HDUs: {len(hdul)}")
		#     for i, hdu in enumerate(hdul):
		#         print(f"HDU {i}: {type(hdu)}, shape={getattr(hdu.data, 'shape', None)}")

		choice2 = int(input("Answer: "))


		if choice2 == 1:
			print()
			print(">>> Awesome!")
			print(">>> Going to the next step...")
			print()
			skip_2_step2 = True
			break

		elif choice2 == 2:
			print()
			print(">>> Open 'run_variable_star_tool.py' in a text editor and edit the first line to the correct prefix")
			print(">>> Then run the code again!")
			print(">>> Exiting...")
			print()
			exit()

		else:
			print("Error: Please select either 1,2" + "\n")

# Process images with dark and flat feilds (if they have them)

# Make Master Bias
if skip_2_step2 == True:
	while True:
		print("================================================================================" +"\n")	
		print("STEP 2:" +"\n")

		print(f">>> Do you want to create a master BIAS? ")
		print(f">>> Bias fits used from the fallowing folder  -->  {path_2_bias}")
		print("(1) Yes, create a new master bias!")
		print("(2) No, reuse the master bias from a previouse run")
		print("(3) No, reuse the master bias from a previouse run & Open the master bias")
		print("(4) Exit" + "\n")

		choice3 = int(input("Answer: "))

		if choice3 == 1:
			print()
			print(">>> Runing 'make_master_bias.py' ..." + "\n")
			
			master_bias_name = target + "_" + date + "_" + "master_bias.fit"

			subprocess.run(["python3","make_master_bias.py", path_2_bias, path_2_processed_bias, master_bias_name])

			subprocess.run(["ds9", f"{path_2_processed_bias}/{master_bias_name}"])

			print()
			print(f">>> Success! Created a master bias --> {path_2_processed_bias}/{master_bias_name}" + "\n")
			skip_2_step3 = True
			break

		elif choice3 == 2:
			print()
			print(">>> Reusing previouse master bais")
			print(">>> Going to the next step...")
			print()
			skip_2_step3 = True
			break

		elif choice3 == 3:
			print()
			print(">>> Reusing previouse master bais")
			print(">>> Opening previouse master bais")

			master_bias_name = target + "_" + date + "_" + "master_bias.fit"

			subprocess.run(["ds9", f"{path_2_processed_bias}/{master_bias_name}"])

			print(">>> Going to the next step...")
			print()
			skip_2_step3 = True
			break

		elif choice3 == 4:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3,4" + "\n")

# Make Master dark
if skip_2_step3 == True:
	while True:
		print("================================================================================" +"\n")
		print("STEP 3: \n")
		
		print(f">>> Do you want to create a master DARK ")
		print(f">>> Dark fits used from the fallowing folder  -->  {path_2_darks}")
		print("(1) Yes, create a new master dark!")
		print("(2) No, reuse the master dark from a previouse run")
		print("(3) No, reuse the master dark from a previouse run & Open the master dark")
		print("(4) Exit" + "\n")

		choice3 = int(input("Answer: "))

		if choice3 == 1:
			print()
			print(">>> Runing 'make_master_dark.py' ..." + "\n")
			
			master_bias_name = target + "_" + date + "_" + "master_bias.fit"
			master_dark_name = target + "_" + date + "_" + "master_dark.fit"

			master_bias_path = f"{path_2_processed_bias}/{master_bias_name}"

			from astropy.io import fits

			# Check exposure times of first dark and first image
			dark_files = sorted(glob.glob(f"{path_2_darks}/*.fit") + 
			                    glob.glob(f"{path_2_darks}/*.fits") +
			                    glob.glob(f"{path_2_darks}/*.FITS"))
			image_files = sorted(glob.glob(f"{path_2_images}/*.fit") + 
			                     glob.glob(f"{path_2_images}/*.fits") +
			                     glob.glob(f"{path_2_images}/*.FITS"))

			if not dark_files or not image_files:
			    print(">>> Error: No FITS files found in darks or image folder.")
			    exit()

			dark_exptime = fits.getheader(dark_files[0])['EXPTIME']
			image_exptime = fits.getheader(image_files[0])['EXPTIME']

			if dark_exptime != image_exptime:
			    print()
			    print(f">>> Warning! Exposure time mismatch:")
			    print(f">>>   Dark frame EXPTIME:  {dark_exptime} s")
			    print(f">>>   Image frame EXPTIME: {image_exptime} s")
			    print()
			    print(">>> Are you sure you wish to continue? The dark will be scaled to match the image exposure.")
			    print("(1) Yes, continue and scale the dark")
			    print("(2) No, exit so I can change the darks" + "\n")
			    
			    scale_choice = int(input("Answer: "))

			    if scale_choice == 2:
			        print(">>> Exiting. Please replace the darks and run again.")
			        exit()

			# Add both exposure times as environment variables for use in make_master_dark.py
			os.environ["DARK_EXPTIME"] = str(dark_exptime)
			os.environ["IMAGE_EXPTIME"] = str(image_exptime)

			# Run the script
			subprocess.run(["python3", "make_master_dark.py", 
							path_2_darks, 
							master_bias_path, 
							path_2_processed_dark, 
							master_dark_name])

			subprocess.run(["ds9", f"{path_2_processed_dark}/{master_dark_name}"])

			print()
			print(f">>> Success! Created a master dark --> {path_2_processed_dark}/{master_dark_name}" + "\n")

			skip_2_step4 = True
			break

		elif choice3 == 2:
			print()
			print(">>> Reusing previouse master bais")
			print(">>> Going to the next step...")
			print()
			skip_2_step4 = True
			break


		elif choice3 == 3:
			print()
			print(">>> Reusing previouse master dark")
			print(">>> Opening previouse master dark")

			master_dark_name = target + "_" + date + "_" + "master_dark.fit"

			subprocess.run(["ds9", f"{path_2_processed_dark}/{master_dark_name}"])

			print(">>> Going to the next step...")
			print()
			skip_2_step4 = True
			break

		elif choice3 == 4:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3,4" + "\n")


# Apply Calibration (dark,bias,flat)
if skip_2_step4 == True:
	while True:
		print("================================================================================" +"\n")	
		print("STEP 4: \n")

		print(f">>> Do you want to apply the master bias, master dark, and flat to the images?")
		print("(1) Yes!")
		print("(2) No, reuse the processed images from a previouse run")
		print("(3) Exit" + "\n")

		choice3 = int(input("Answer: "))

		if choice3 == 1:
			print()
			print(">>> Runing 'apply_calibrations.py' ..." + "\n")
			
			master_bias_name = target + "_" + date + "_" + "master_bias.fit"
			master_dark_name = target + "_" + date + "_" + "master_dark.fit"

			master_bias_path = f"{path_2_processed_bias}/{master_bias_name}"
			master_dark_path = f"{path_2_processed_dark}/{master_dark_name}"

			subprocess.run(["python3","apply_calibrations.py", path_2_images, master_bias_path, master_dark_path, path_2_flat, path_2_processed])
			skip_2_step5 = True
			break

		elif choice3 == 2:
			print()
			print(">>> Reusing previouse pre-processed images")
			print(">>> Going to the next step...")
			print()
			skip_2_step5 = True
			break

		elif choice3 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")


used_astrometry = False

# Get WCS for each of the pre-processed images
if skip_2_step5 == True:
	while True:
		print("================================================================================" +"\n")	
		print("STEP 5: \n")

		print(f">>> The next step is to get the wcs information from astronomy.net. Ready to continue?")
		print("(1) Yes! --> Use astronomy.net")
		print("(2) Yes! --> Use Scamp (Beta)")
		print("(3) No, reuse the WCS images from a previouse run")
		print("(4) Exit" + "\n")

		choice4 = int(input("Answer: "))

		if choice4 == 1:

			while True:
				print(f">>> Pick how you would like to upload the image to astronomy.net")
				print("(1) Upload the whole image (CMT/ARCT)")
				print("(2) Split each image up into 16 quadrents then upload each quadrent to astronomy.net (KGMT) \n")
				print("\nNote: a good WCS is needed to track stars though the images\n")

				choice4a = int(input("Answer: "))
				print()

				if choice4a == 1:
					print()
					print(">>> Running 'get_wcs_wo_splitting.py' ...")

					subprocess.run(["python3", "get_wcs_wo_splitting.py",  # uploads the whole image to astrometry
									path_2_processed, 
									path_2_processed_wcs])

					used_astrometry = True

					print("\n >>> Successfuly solved the images with astronomy.net!")
					print(">>> Going to the next step... \n")
					skip_2_step6 = True
					break

				if choice4a == 2:
					print(">>> Running 'get_wcs.py' ...") # split into 16 quadrants then upload to astrometry

					subprocess.run(["python3", "get_wcs.py",
									path_2_processed, 
									path_2_processed_wcs])

					used_astrometry = True

					print("\n >>> Successfuly solved the images with astronomy.net!") 
					print(">>> Going to the next step... \n")
					break

				else: 
					print("Error: Please select either 1 or 2" + "\n")
			break

		if choice4 == 2:
			print()
			print(">>> First we must get a rough wcs so that scamp will work")
			print("(1) Add rough wcs")
			print("(2) Reuse previouse rough wcs")
			print()
			reuse_rough_wcs = False
			reuse_rough_wcs_choice = int(input("Answer: "))

			if reuse_rough_wcs_choice == 1:
				reuse_rough_wcs = False

			elif reuse_rough_wcs_choice == 2:
				reuse_rough_wcs = True

			if reuse_rough_wcs == False:
				print()
				print(">>> Running 'get_rough_wcs.py' ...")

				subprocess.run(["python3", "get_rough_wcs.py", 
								 path_2_processed])

				print("Completed getting the rough wcs")

			print()
			print(">>> Running 'get_wcs_scamp.py' ...")

			subprocess.run(["python3", "get_wcs_scamp.py", 
							path_2_processed,				# 1 path to processed images
							path_2_processed_wcs, 			# 2 path to output wcs images
							path_2_scamp_config, 			# 3 path to scamp config dir
							path_2_scamp_output,			# 4 path to scamp output dir
							path_2_se_wcs_config,			# 5 path to sextractor config (for wcs)
							path_2_se_wcs_output, 			# 6 path to sextractor output (for wcs) (used for scamp)
							path_2_swarp_config])			# 7 path to swarp config

			print("\n >>> Successfuly solved the images with Scamp!")
			print(">>> Going to the next step... \n")
			skip_2_step6 = True
			used_astrometry = False
			break

		elif choice4 == 3:
			print()
			print(">>> Reusing WCS processed images")
			print(">>> Going to the next step...")
			print()
			skip_2_step6 = True
			break


		elif choice4 == 4:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")


# Use source extractor to get photometry on each of the sources
if skip_2_step6 == True:
	while True:
		print("================================================================================" +"\n")	
		print("STEP 6: \n")

		print(f">>> The next step is to do photometry on the WCS'ed images. Pick your choice.")
		print("(1) Aperture photometry (Source-Extractor)")
		print("(2) Aperture photometry with backgound subtraction! (Source-Extractor 2x)")
		print("(3) Aperture photometry (Source-Extractor) & PSF photometry (PSFEx)! (Beta)")
		print("(4) No, reuse the photometry from a previouse run")
		print("(5) Exit" + "\n")

		choice5 = int(input("Answer: "))

		if choice5 == 1:
			print()
			print(">>> Running 'source-extract.py' ...")
			config_name = "/kgmt.config"
			checkplot = "n"
			subprocess.run(["python3", "source-extract.py", 
							path_2_processed_wcs, 
							path_2_se_config, 
							path_2_se_output,
							config_name,
							checkplot,
							path_2_background_imgs])

			print("\n >>> Successfuly did photometry on the images!")
			print(">>> Going to the next step... \n")
			skip_2_step7 = True
			break

		elif choice5 == 2:
			os.makedirs(path_2_se_output_v1, exist_ok=True)
			os.makedirs(path_2_background_imgs, exist_ok=True)
			config_name = "/kgmt_back.config"
			checkplot = "y"
			subprocess.run(["python3", "source-extract.py", 
							path_2_processed_wcs, 
							path_2_se_config, 
							path_2_se_output_v1,
							config_name,
							checkplot,
							path_2_background_imgs])

			subprocess.run(["python3", "sub_se_back.py",
							path_2_processed_wcs,
							path_2_se_output_v1,
							path_2_background_imgs])

			config_name = "/kgmt.config"
			checkplot = "n"
			subprocess.run(["python3", "source-extract.py", 
							path_2_background_imgs, 
							path_2_se_config, 
							path_2_se_output,
							config_name,
							checkplot,
			 				path_2_background_imgs])

			skip_2_step7 = True
			break

		elif choice5 == 3:
			print()
			print(">>> psf_phot.py ...")

			subprocess.run(["python3", "psf_phot.py",
							path_2_processed_wcs,
							path_2_se_psf_run1_config,
							path_2_se_psf_run1_output,
							path_2_psfex_config,
							path_2_psfex_output,
							path_2_se_psf_run2_config,
							path_2_se_psf_run2_output])

		elif choice5 == 4:
			print()
			print(">>> Reusing photometricaly processed image catalogues")
			print(">>> Going to the next step...")
			print()
			skip_2_step7 = True
			break


		elif choice5 == 5:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3,4,5" + "\n")


# combine quadrents .ldac to csv to make one big csv for each image. 
# if skip_2_step6 == True and used_astrometry == True:
if skip_2_step7 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> We need to turn the .ldac files to .csv")
		print("(1) My images are in 16 quadrents (Astrometry)")
		print("(2) My images are in one pice (Scamp)")
		print()
		choice6 = int(input("Answer: "))

		if choice6 == 1:
			print(">>> Combining .ldac Source-Extractor catalogues from each quadrent of the original image into one big catalogue! ")
			print(">>> Running 'combine_ldac_quadrants_2_csv.py' ...")

			subprocess.run(["python3", "combine_ldac_quadrants_2_csv.py",
							path_2_se_output,
							path_2_combined_star_cats])
			skip_2_step7 = True
			break

		if choice6 == 2:
			print(">>> Changing the catalogues from .ldac to .csv")
			print(">>> Running 'ldac_2_csv.py")

			subprocess.run(["python3", "ldac_2_csv.py",
							path_2_se_output,
							path_2_combined_star_cats])


			skip_2_step7 = True
			break

		else:
			print("Error: Please select either 1,2" + "\n")

# Filter stars out of sources from each image
if skip_2_step7 == True:
	while True:
		print("================================================================================" +"\n")	
		print(f">>> The next step is to filter stars out of the Source-Extractor catalogues. Ready to continue?")
		print("(1) Yes!")
		print("(2) No, reuse the star filtering from a previouse run")
		print("(3) Exit" + "\n")

		choice7 = int(input("Answer: "))

		#If using aper
		# mag_faint = 10.4
		# mag_bright = 7

		#If using auto
		# mag_faint = 9.5
		# mag_bright = 5.5

		if choice7 == 1:
			
			while True:
				while True:
					print(">>> Set a value for the faintest magnidude (MAG_AUTO) to be considered a star.")
					print("    Suggested (KGMT) 9.5")
					mag_faint = input("mag_faint = ")
					print()

					print(">>> Set a value for the brightest magnidude (MAG_AUTO) to be considered a star.")
					print("    Suggested (KGMT) 5.5")
					mag_bright = input("mag_bright = ")
					print()

					print(">>> Set a value for the maximum number of flags to be considered a star.")
					print("    Suggested (KGMT) 3")
					flags = input("flags = ")
					print()

					print(">>> Set a value for the number to multipy the median flux radius by to set the min radiuse to be considered a star.")
					print("    Suggested (KGMT) 0.55")
					min_flux_radiuse = input("min_flux_radiuse = ")
					print()

					print(">>> Set a value for the number to multipy the median flux radius by to set the max radiuse to be considered a star.")
					print("    Suggested (KGMT) 1.45")
					max_flux_radiuse = input("max_flux_radiuse = ")
					print()

					print(">>> Set a value for the number minimum elongation to be considered a star.")
					print("    1 = circle, >1 = elliptical")
					print("    Suggested (KGMT) 1")
					min_elongation = input("min_elongation = ")
					print()

					print(">>> Do you want to use these values:")
					print(f"mag_faint        = {mag_faint}")
					print(f"mag_bright       = {mag_bright}")
					print(f"max_flags        = {flags}")
					print(f"min_flux_radiuse = {min_flux_radiuse}")
					print(f"max_flux_radiuse = {max_flux_radiuse}")
					print(f"min_elongation   = {min_elongation}")
					print()
					print("(1) Continue, I am happy with these values")
					print("(2) Change values\n")

					choice7a = int(input("Answer: "))

					if choice7a == 1:
						print(">>> Using set values ... ")
						break

					if choice7a == 2:
						print(">>> Going back to change values ...")


				print()
				print(">>> Running 'star_select.py' ...")

				subprocess.run(["python3", "star_select.py",
							path_2_combined_star_cats,
							path_2_combined_star_only_cats,
							str(mag_faint),
							str(mag_bright),
							str(flags),
							str(min_flux_radiuse),
							str(max_flux_radiuse),
							str(min_elongation),
							])

				print("\n>>> Successfuly did photometry on the images!\n")

				print(">>> Check the star_cat_plots to see if you are happy with the star selection")
				print("(1) Continue to next step, im happy with this selection")
				print("(2) Redo star selection")

				choice7c = int(input("Answer: "))

				if choice7c == 1:
					print(">>> Going to next step ...")
					break

				if choice7c == 2:
					print(">>> Going back to star selection ...")

			print()
			while True:
				print(">>> Running 'plot_vignet_profiles.py' ...\n")

				print(">>> Use star filtered catalogues or the unfiltered catalogues?")
				print("(1) Star filtered")
				print("(2) Unfiltered sources")

				choice7b = int(input("Answer: "))

				if choice7b == 1:
					print("\n Using star filtered catalogues\n")
					star_filtered = "filtered"
					subprocess.run(["python3", "plot_vignet_profiles.py", 
									path_2_combined_star_only_cats, 
									star_filtered])

				if choice7b == 2:
					print("\n Using unfiltered catalogues\n")
					star_filtered = "unfiltered"
					subprocess.run(["python3", "plot_vignet_profiles.py", 
									path_2_combined_star_cats, 
									star_filtered])

				print(">>> Would you like to plot more vignet profiles?")
				print("(1) Yes")
				print("(2) No\n")

				choice7a = int(input("Answer: "))

				if choice7a == 1:
					print(">>> Another!")

				elif choice7a == 2:
					print(">>> Not creating other plots")
					break

			print(">>> Going to the next step... \n")
			skip_2_step8 = True
			break

		elif choice7 == 2:
			print()
			print(">>> Reusing photometricaly processed image catalogues")
			print(">>> Going to the next step...")
			print()
			skip_2_step8 = True
			break


		elif choice7 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")



# Fallow stars along in the cataloges to make a catalogue of that star over time
if skip_2_step8 == True:
	while True:
		print("================================================================================" +"\n")	
		print(f">>> The next step is track stars thoughout the images!\n")
		print(">>> What search radiuse would you like to use to find star matches?")
		print(">>> NOTE: This will effect how many stars we can track\n")
		print("Suggested for the KGMT = 3 arcsec")
		
		search_radiuse = int(input("Search Radiuse (arcsec): "))

		print()

		print(">>> Run new tracking script with magnidudes and errors checked or just track stars based on location?")
		print("(1) Filter just based on location")
		print("(2) Run new tracking script")		
		print()

		choice6a = int(input("Answer: "))

		print()
		print(">>> Track stars though star filtered catalogues or on non filtered catalogues?")
		print("(1) Star filtered catalogues")
		print("(2) Unfiltered catalogues")
		print()

		choice6b = int(input("Answer: "))
		print()
		print(">>> Through what is the minimum percentage you would like to track though the images?")
		print("    (ex: 100 for 100%, 80 for 80%, etc)\n")

		min_track_percent = int(input("minimum tracking percentage = "))
		print()
		
		if choice6a == 1:
			print(">>> Running 'track_stars.py' ...")

			if choice6b == 1:

				print(">>> Using star filtered catalogues")
				subprocess.run(["python3", "track_stars.py", 
							path_2_combined_star_only_cats, 
							path_2_star_tracked_cat, 
							str(search_radiuse), 
							path_2_processed,
							str(min_track_percent)])

			if choice6b == 2:
				print(">>> Using unfiltered catalogues")
				subprocess.run(["python3", "track_stars.py", 
								path_2_combined_star_cats, 
								path_2_star_tracked_cat, 
								str(search_radiuse), 
								path_2_processed,
								str(min_track_percent)])

		if choice6a == 2:
			print(">>> Running 'new_track_stars.py' ...")

			if choice6b == 1:
				print(">>> Using star filtered catalogues")
				subprocess.run(["python3", "new_track_stars.py", 
							path_2_combined_star_only_cats, 
							path_2_star_tracked_cat, 
							str(search_radiuse), 
							path_2_processed,
							str(min_track_percent)])

			if choice6b == 2:
				print(">>> Using unfiltered catalogues")
				subprocess.run(["python3", "new_track_stars.py", 
								path_2_combined_star_cats, 
								path_2_star_tracked_cat, 
								str(search_radiuse), 
								path_2_processed,
								str(min_track_percent)])
		
		print("\n>>> Would you like to run this again with a diffrent search radiuse?")
		print("(1) No, I am happy with this result")
		print("(2) Yes, re-run tracking script")
		print("(3) Exit\n")

		choice6 = input("Answer: ")


		if choice6 == "2":
			print()
			print(">>> Going back")
			print()
		
		elif choice6 == "1":
			print()
			print(">>> Going to the next step...")
			print()
			skip_2_step9 = True
			break


		elif choice6 == '3':
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")


if skip_2_step9 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to create star tracking check plots?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice15a = int(input("Answer: "))
		print()

		if choice15a == 1:
			print(">>> Running 'plot_star_index_on_imgs.py'") 
			path_2_tracking_check_plots = path_2_star_tracked_cat + "/tracking_check_plots/star_index_flagged_imgs"
			os.makedirs(path_2_tracking_check_plots, exist_ok=True)
			subprocess.run(["python3", "plot_star_index_on_imgs.py", 
							path_2_star_tracked_cat, 
							path_2_processed_wcs, 
							path_2_tracking_check_plots, 
							target, 
							date])

			print(">>> Running 'plot_tracking_percent.py'") 
			path_2_tracking_check_plots = path_2_star_tracked_cat + "/tracking_check_plots/tracking_percent"
			os.makedirs(path_2_tracking_check_plots, exist_ok=True)
			subprocess.run(["python3", "plot_tracking_percent.py", 
							path_2_star_tracked_cat, 
							path_2_processed_wcs, 
							path_2_tracking_check_plots, 
							target, 
							date])

			print(">>> Running 'plot_when_tracked.py'") 
			path_2_tracking_check_missed_plots = path_2_star_tracked_cat + "/tracking_check_plots/tracking_map"
			os.makedirs(path_2_tracking_check_missed_plots, exist_ok=True)
			subprocess.run(["python3", "plot_when_tracked.py", 
							path_2_star_tracked_cat, 
							path_2_tracking_check_missed_plots, 
							target, 
							date])

			print(">>> Successfuly created tracking star check plots")
			print(">>> Going to the next step...")
			skip_2_step10 = True
			break

		elif choice15a == 2:
			print()
			print(">>> Not making any check plots :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step10 = True
			break


		elif choice15a == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

if skip_2_step10 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to create a star timelaps to check star tracking?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice15 = int(input("Answer: "))
		print()

		if choice15 == 1:

			while True:
				print("================================================================================" +"\n")
				print(">>> Which star would you like to create a star timelaps?\n")

				star_index = int(input("Star index: "))
				print()

				print(">>> Running 'star_track_video.py'")
				
				subprocess.run(["python3", "star_track_video.py",
								path_2_star_tracked_cat,  # data dir path to tracked stars csv
								path_2_star_tracked_cat + f"/star_timelaps", # path to output vid
								target, # target id
								date,   # date
								str(star_index), # desired star index 
								path_2_processed_wcs]) # fits folder

				print(f">>> Created timelaps for star {star_index}\n")
				print(f">>> Would you like to create another star timelaps?")
				print("(1) Yes!")
				print("(2) No :(")

				print()
				create_another = int(input("Answer: "))
				print()

				if create_another == 1:
					print(">>>  Another!")

				if create_another == 2:
					print(">>> Sad :(")
					skip_2_step11 = True
					break
				

		elif choice15 == 2:
			print()
			print(">>> Not making any videos :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step11 = True
			break


		elif choice15 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()


# Get star classifications from simbad
if skip_2_step11 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> The next step is to find the star classifications from SIMBAD. Ready to continue?")
		print("(1) Yes!")
		print("(2) No, reuse the star classifications from a previouse run")
		print("(3) No, skip this step set every star to UNCLASSIFIED")
		print("(4) Exit" + "\n")

		choice9 = int(input("Answer: "))


		if choice9 == 1:
			print()

			print(">>> What search radiuse would you like to use to find star in Simbad?")
			print("    (Distance in arcsec from our catalogue to Simbad's data)")
			print("     Suggested: 8-10\n")
			radiuse_arcsec = int(input("Search Radiuse (arcsec) = "))

			print(">>> Running 'find_SIMBAD_class.py' ...")

			subprocess.run(["python3", "find_SIMBAD_class.py", 
							 path_2_star_tracked_cat, 
							 path_2_star_tracked_simbad_cat, 
							 str(radiuse_arcsec)])

			print("\n>>> Successfuly found star classifications with SIMBAD!")

			print("\n============ Answer Key For SIMBAD Types ============\n")
			print("EB* =  Eclipsing Binary")
			print("*   =  Star")
			print("cC* =  Classical Cepheid Variable Star")
			print("PM* =  High Proper Motion Star")
			print("WD* =  White Dwarf")
			print("LM* =  Low-Mass Star")
			print("WD* =  White Dwarf Star")
			print("WD? =  White Dwarf Candidate")
			print("RR* =  RR Lyrae Variable")
			print("SB* =  Spectroscopic Binary")
			print("GrG =  Gravitationally Lensed Galaxy")
			print("Em* =  Emission-line Star")
			print("RS* =  Radio Star")
			print("**  =  Double Star")
			print("HB* =  Horizontal Branch Star")
			print("RG* =  Red Giant Star")
			print("dS* =  Dwarf Star")
			print("Sy1 =  Seyfert 1 Galaxy")
			print("Sy2 =  Seyfert 2 Galaxy")
			print("G   =  Galaxy")
			print("GiG =  Giant Galaxy")
			print("Pe* =  Peculiar Star")
			print("QSO =  Quasar")
			print("C*? =  Carbon Star Candidate")
			print("ds* =  Dwarf Spheroidal Galaxy")
			print("PaG =  Pair of Galaxies")
			print("BY* =  BY Draconis Variable")
			print("LP* =  Long Period Variable")
			print("LP? =  Long Period Variable Candidate")
			print("LYN =  Low-Ionization Nuclear Emission-line Region")
			print("V*  =  Variable Star")
			print("Ro* =  Rotating Variable Star")
			print("Er* =  Eruptive Variable")
			print("Y*? =  Young Stellar Object Candidate")
			print("HS? =  High-mass Star Candidate")
			print("\n=====================================================\n")


			print(">>> Going to the next step... \n")
			skip_2_step12 = True
			break

		elif choice9 == 2:
			print()
			print(">>> Reusing star classifications with SIMBAD")
			print(">>> Going to the next step...")
			print()
			skip_2_step12 = True
			break

		elif choice9 == 3:
			print()
			print(">>> Running 'skip_SIMBAD_class.py' ...")

			subprocess.run(["python3", "skip_SIMBAD_class.py", 
							 path_2_star_tracked_cat, 
							 path_2_star_tracked_simbad_cat, 
							 ])

			print(">>> Going to the next step... \n")
			skip_2_step12 = True
			break

		elif choice9 == 4:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")


# full_star_data = pd.read_csv(path_2_star_tracked_simbad_cat + "/tracked_stars_with_simbad.csv")
# unique_star_indices = full_star_data['star_index'].unique()
# print(f"Found {len(unique_star_indices)} unique stars")
# first_observations = full_star_data.drop_duplicates(subset=['star_index'])
# first_observations.set_index('star_index', inplace=True)

# Plot the quadrents with sorted stars and star types from simbad
if skip_2_step12 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to plot the individual image quandrents with an overlay of SIMBAD classified stars?")
		print(">>> For just the first image")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice10 = int(input("Answer: "))

		if choice10 == 1:
			print()
			print(">>> Running 'plot_simbad_star_class.py' ...")

			subprocess.run(["python3", "plot_simbad_star_class.py", 
							path_2_star_tracked_simbad_cat, 
							path_2_processed_wcs, 
							path_2_simbad_star_class_plots, 
							target, 
							date])

			print(">>> Successfuly created simbad star classification plots")
			print(">>> Going to the next step...")
			skip_2_step13 = True
			break

		elif choice10 == 2:
			print()
			print(">>> Not making any plots :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step13 = True
			break


		elif choice10 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")



if skip_2_step13 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to plot the stars magnidude over time?")
		print(">>> This shows stars trends before differential photometry")
		print(">>> This will generate multi plots for UNCLASSIFIED, Variable, and Misc star groups")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice11 = int(input("Answer: "))

		if choice11 == 1:
			print()
			print(">>> Running 'plot_lightcurves_wo_diffphot.py' ...")

			subprocess.run(["python3", "plot_lightcurves_wo_diffphot.py",
							 path_2_star_tracked_simbad_cat + "/tracked_stars_with_simbad.csv", 
							 path_2_light_curves_wo_diffphot, 
							 target, 
							 date])

			print(">>> Successfuly created simbad star classification plots")
			print(">>> Going to the next step...")
			skip_2_step14 = True
			break

		elif choice11== 2:
			print()
			print(">>> Not making any plots :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step14 = True
			break


		elif choice11 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")

if skip_2_step14 == True:	
	diff_phot_step = 1
	print("\n================================================================================")
	print("                         Begining Variable Star Analysis                        ")
	print("================================================================================" +"\n")	
	print(">>> Differential Photometry Section Steps")
	print("(1) Check for variable stars defined by SIMBAD & Generate analysis log")
	print("(2) Specify known variable and non-variable star indices")
	print("(3) 1st run of differential photometry")
	print("(4) Find new variables based on average variability")
	print("(5) 2nd run of differential photometry\n")

	print(">>> Which step would you like to start at? (NOTE: If this is your first run, you must start at 1)\n")
	diff_phot_step = int(input("Answer: "))	
	print()
	if diff_phot_step == 1:
		print(">>> Adding new column to csv to define variable & non-variable stars")
		print("\n >>> Running 'variable_star_cheak.py' ...")

		subprocess.run(["python3", "variable_star_cheak.py",
						path_2_star_tracked_simbad_cat, 
						path_2_csv_variables_defined, 
						target, 
						date])

		print("\n>>> Successfuly added the variability flag for known variable stars from SIMBAD")

		print("\n>>> Running 'generate_analysis_log.py")

		subprocess.run(["python3", "generate_analysis_log.py",
						path_2_csv_variables_defined,	# Input csv path
						path_2_csv_variables_defined,	# Ouput log path
						target,
						date])

		print("\n>>> Running 'clean_mag_errors.py' ... ")
		subprocess.run(["python3", "clean_mag_errors.py", 
						path_2_csv_variables_defined, 
						target, 
						date])

		diff_phot_step = 2

	if diff_phot_step == 2:
		print("\n >>> Running 'manually_update_variability.py' ...")

		subprocess.run(["python3", "manually_update_variability.py",
						path_2_csv_variables_defined,
						path_2_csv_variables_defined,
						target,
						date])
		print(">>> Done manually updating csv ...")
		diff_phot_step = 3

	if diff_phot_step == 3:
		print("\n>>> Now for the fun part :)")
		print(">>> It is differential photometry time!")
		print()
		print("\n>>> Running 'diff_phot.py' ... ")
		csv_name = f"tracked_stars_with_variability_flag_{target}_{date}.csv"
		print(f">>> Using csv: {csv_name}")
		subprocess.run(["python3", "diff_phot.py", 
						path_2_csv_variables_defined, 
						path_2_csv_diff_phot, 
						target, 
						date,
						csv_name])

		print("\n>>> Successfuly did differential photometry!")

		print("\n>>> Running 'clean_diff_phot_errors.py' ... ")
		subprocess.run(["python3", "clean_diff_phot_errors.py", 
						path_2_csv_diff_phot, 
						target, 
						date])
		diff_phot_step = 4

	if diff_phot_step == 4:
		print("\n>>> Running 'find_new_variables.py")
		subprocess.run(["python3", "find_new_variables.py", 
						path_2_csv_diff_phot,
						path_2_csv_variables_defined,
						target,
						date])

		print("\n>>> Successfuly found new variables!")
		diff_phot_step = 5

	if diff_phot_step == 5:
		print("\n>>> Running 'diff_phot.py' ... ")

		csv_name = f"tracked_stars_with_updated_variability_flag_{target}_{date}.csv"
		print(f">>> Using csv: {csv_name}")
		subprocess.run(["python3", "diff_phot.py", 
						path_2_csv_variables_defined, 
						path_2_csv_diff_phot, 
						target, 
						date,
						csv_name])

		print("\n>>> Running 'clean_diff_phot_errors.py' ... ")
		subprocess.run(["python3", "clean_diff_phot_errors.py", 
						path_2_csv_diff_phot, 
						target, 
						date])

		print("\n>>> Successfuly did differential photometry!")


		print("\n>>> Going to the next step ...")

	skip_2_step15 = True

if skip_2_step15 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to plot the stars magnidude over time?")
		print(">>> This shows stars trends after differential photometry")
		print(">>> This will generate multi plots for UNCLASSIFIED, Variable, and Misc star groups")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice13 = int(input("Answer: "))

		if choice13 == 1:
			print()
			while True:
			    print(">>> Use weighted or unweighted differential magnitudes to find new variables?")
			    print("(1) unweighted")
			    print("(2) weighted\n")

			    which_mag = int(input("Answer: "))
			    print()

			    if which_mag == 1:
			        DIFF_MAG_COLUMN = 'differential_mag'
			        DIFF_MAG_ERROR_COLUMN = 'differential_mag_err'
			        w = "uw"
			        break
			    if which_mag == 2:
			        DIFF_MAG_COLUMN = 'weighted_differential_mag'
			        DIFF_MAG_ERROR_COLUMN = 'weighted_differential_mag_err'
			        w = "w"
			        break
			    else:
			        print("Error: Please select either 1, 2" + "\n")

			# Display selected magnitude columns
			print(">>> Magnitude option selection")
			print(f"    DIFF_MAG_COLUMN: {DIFF_MAG_COLUMN}")
			print(f"    DIFF_MAG_ERROR_COLUMN: {DIFF_MAG_ERROR_COLUMN}\n")

			print(">>> Running 'plot_lightcurves_w_diffphot.py' ...")

			subprocess.run(["python3", "plot_lightcurves_w_diffphot.py", 
							path_2_csv_diff_phot, 
							path_2_light_curves_w_diffphot, 
							target, 
							date,
							DIFF_MAG_COLUMN,
							DIFF_MAG_ERROR_COLUMN,
							w])

			print(">>> Running 'plot_lightcurves_w_diffphot_same_scale.py' ...")

			subprocess.run(["python3", "plot_lightcurves_w_diffphot_same_scale.py", 
							path_2_csv_diff_phot, 
							path_2_light_curves_w_diffphot, 
							target, 
							date,
							DIFF_MAG_COLUMN,
							DIFF_MAG_ERROR_COLUMN,
							w])
			
			print(">>> Running 'plot_both_lightcurves.py'")
			subprocess.run(["python3", "plot_both_lightcurves.py", 
							path_2_csv_diff_phot, 
							path_2_light_curves_w_diffphot, 
							target, 
							date,
							DIFF_MAG_COLUMN,
							DIFF_MAG_ERROR_COLUMN,
							w])


			print(">>> Successfuly created light curves witht the differential photometry data")
			
			while True:
				print(">>> Would you like to plot a light curve multi-plot starting at a specific star index?")
				print("(1) Yes!")
				print("(2) No\n")
				choice13a = int(input("Answer: "))
				print()
				while True:
				    print(">>> Use weighted or unweighted differential magnitudes to find new variables?")
				    print("(1) unweighted")
				    print("(2) weighted\n")

				    which_mag = int(input("Answer: "))
				    print()

				    if which_mag == 1:
				        DIFF_MAG_COLUMN = 'differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'differential_mag_err'
				        w = "uw"
				        break
				    if which_mag == 2:
				        DIFF_MAG_COLUMN = 'weighted_differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'weighted_differential_mag_err'
				        w = "w"
				        break
				    else:
				        print("Error: Please select either 1, 2" + "\n")
				if choice13a == 1:
					print(">>> Running 'plot_lightcurves_w_diffphot_same_scale_specified_index.py'")
					subprocess.run(["python3", "plot_lightcurves_w_diffphot_same_scale_specified_index.py", 
									path_2_csv_diff_phot, 
									path_2_light_curves_w_diffphot, 
									target, 
									date,
									DIFF_MAG_COLUMN,
									DIFF_MAG_ERROR_COLUMN,
									w])

				if choice13a == 2:
					print("Not makeing any more plots ...")
					break

			print(">>> Going to the next step...")
			skip_2_step16 = True
			break

		elif choice13== 2:
			print()
			print(">>> Not making any plots :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step16 = True
			break


		elif choice13 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

		else:
			print("Error: Please select either 1,2,3" + "\n")


if skip_2_step16 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to plot a star profile?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice14 = int(input("Answer: "))
		print()

		if choice14 == 1:

			while True:
				print("================================================================================" +"\n")
				print(">>> Which star would you like to create a profile for?\n")

				star_index = int(input("Star index: "))
				print()
				while True:
				    print(">>> Use weighted or unweighted differential magnitudes to find new variables?")
				    print("(1) unweighted")
				    print("(2) weighted\n")

				    which_mag = int(input("Answer: "))
				    print()

				    if which_mag == 1:
				        DIFF_MAG_COLUMN = 'differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'differential_mag_err'
				        w = "uw"
				        break
				    if which_mag == 2:
				        DIFF_MAG_COLUMN = 'weighted_differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'weighted_differential_mag_err'
				        w = "w"
				        break
				    else:
				        print("Error: Please select either 1, 2" + "\n")

				# Display selected magnitude columns
				print(">>> Magnitude option selection")
				print(f"    DIFF_MAG_COLUMN: {DIFF_MAG_COLUMN}")
				print(f"    DIFF_MAG_ERROR_COLUMN: {DIFF_MAG_ERROR_COLUMN}\n")

				print(">>> Running 'star_profile.py'")
				
				subprocess.run(["python3", "star_profile.py", 
								path_2_csv_diff_phot, 
								path_2_light_curves_w_diffphot + f"/star_profiles", 
								target, 
								date, 
								str(star_index), 
								path_2_processed_wcs,
								DIFF_MAG_COLUMN,
								DIFF_MAG_ERROR_COLUMN,
								w])

				print(f">>> Created profile for star {star_index}\n")
				print(f">>> Would you like to create another star profile?")
				print("(1) Yes!")
				print("(2) No :(")

				print()
				create_another = int(input("Answer: "))
				print()

				if create_another == 1:
					print(">>>  Another!")

				if create_another == 2:
					print(">>> Sad :(")
					skip_2_step17 = True
					break
				

		elif choice14 == 2:
			print()
			print(">>> Not making any plots :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step17 = True
			break


		elif choice14 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

if skip_2_step17 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to create a star timelaps?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice15 = int(input("Answer: "))
		print()

		if choice15 == 1:

			while True:
				print("================================================================================" +"\n")
				print(">>> Which star would you like to create a star timelaps?\n")

				star_index = int(input("Star index: "))
				print()

				while True:
				    print(">>> Use weighted or unweighted differential magnitudes to find new variables?")
				    print("(1) unweighted")
				    print("(2) weighted\n")

				    which_mag = int(input("Answer: "))
				    print()

				    if which_mag == 1:
				        DIFF_MAG_COLUMN = 'differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'differential_mag_err'
				        w = "uw"
				        break
				    if which_mag == 2:
				        DIFF_MAG_COLUMN = 'weighted_differential_mag'
				        DIFF_MAG_ERROR_COLUMN = 'weighted_differential_mag_err'
				        w = "w"
				        break
				    else:
				        print("Error: Please select either 1, 2" + "\n")

				# Display selected magnitude columns
				print(">>> Magnitude option selection")
				print(f"    DIFF_MAG_COLUMN: {DIFF_MAG_COLUMN}")
				print(f"    DIFF_MAG_ERROR_COLUMN: {DIFF_MAG_ERROR_COLUMN}\n")
				print(">>> Running 'star_video.py'")
				
				subprocess.run(["python3", "star_video.py", 
								path_2_csv_diff_phot, 
								path_2_light_curves_w_diffphot + f"/star_timelaps", 
								target, 
								date, 
								str(star_index), 
								path_2_processed_wcs,
								DIFF_MAG_COLUMN,
								DIFF_MAG_ERROR_COLUMN,
								w])

				print(f">>> Created timelaps for star {star_index}\n")
				print(f">>> Would you like to create another star timelaps?")
				print("(1) Yes!")
				print("(2) No :(")

				print()
				create_another = int(input("Answer: "))
				print()

				if create_another == 1:
					print(">>>  Another!")

				if create_another == 2:
					print(">>> Sad :(")
					skip_2_step18 = True
					break
				

		elif choice15 == 2:
			print()
			print(">>> Not making any videos :(")
			print(">>> Going to the next step...")
			print()
			skip_2_step18 = True
			break


		elif choice15 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

if skip_2_step18 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to save a star index to its own csv?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice18 = int(input("Answer: "))
		print()

		if choice18 == 1:

			while True:
				print("================================================================================" +"\n")
				print(">>> Which star would you like to save a star csv?\n")

				star_index = int(input("Star index: "))
				print()

				print(">>> Running 'save_star_index_csv.py'")

				# Saves one star index that has been tracked though the fits
				subprocess.run(["python3", "save_star_index_csv.py", 
								path_2_csv_diff_phot, 
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index), 
								])

				print(f">>> Created csv for star {star_index}\n")

				# Plot with error colour bar the curve and the image numbers
				print("\n>>> Running 'plot_single_star_csv_check.py'")
				version = "1"
				subprocess.run(["python3", "plot_single_star_csv_check.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"differential_mag",
								"differential_mag_err",
								"uw"])

				subprocess.run(["python3", "plot_single_star_csv_check.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"weighted_differential_mag",
								"weighted_differential_mag_err",
								"w"])

				# Plot a nice light cuve of the star alone
				print("\n>>> Running 'plot_single_star_csv.py'")
				subprocess.run(["python3", "plot_single_star_csv.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"differential_mag",
								"differential_mag_err",
								"uw" 
								])

				subprocess.run(["python3", "plot_single_star_csv.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"weighted_differential_mag",
								"weighted_differential_mag_err",
								"w" 
								])

				# Remove user specified rows from the csv of bad points
				print("\n>>> Running 'remove_points.py'")
				subprocess.run(["python3", "remove_points.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index), 
								])

				# Plot with error colour bar the curve and the image numbers after removing points				
				print("\n>>> Running 'plot_single_star_csv_check.py'")
				version = "2"
				subprocess.run(["python3", "plot_single_star_csv_check.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"differential_mag",
								"differential_mag_err",
								"uw"])

				subprocess.run(["python3", "plot_single_star_csv_check.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"weighted_differential_mag",
								"weighted_differential_mag_err",
								"w"])

				# Plot a nice light cuve of the star alone after removing points
				print("\n>>> Running 'plot_single_star_csv.py'")
				subprocess.run(["python3", "plot_single_star_csv.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"differential_mag",
								"differential_mag_err",
								"uw" 
								])

				subprocess.run(["python3", "plot_single_star_csv.py",  
								path_2_indiv_star_study + f"/star_{star_index}", 
								target, 
								date, 
								str(star_index),
								version,
								"weighted_differential_mag",
								"weighted_differential_mag_err",
								"w" 
								])
				print(f">>> Would you like to save another star index csv?")
				print("(1) Yes!")
				print("(2) No :(")

				print()
				create_another = int(input("Answer: "))
				print()

				if create_another == 1:
					print(">>>  Another!")

				if create_another == 2:
					print(">>> Sad :(")
					skip_2_step19 = True
					break
				

		elif choice18 == 2:
			print()
			print(">>> Going to the next step...")
			print()
			skip_2_step19 = True
			break


		elif choice18 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()

if skip_2_step19 == True:
	while True:
		print("================================================================================" +"\n")	
		print(">>> (Optional) would you like to fit an exoplanet curve to a star index?")
		print("(1) Yes!")
		print("(2) No, continue")
		print("(3) Exit" + "\n")

		choice19 = int(input("Answer: "))
		print()

		if choice19 == 1:

			while True:
				print("================================================================================" +"\n")
				print(">>> Which star would you like to fit an exoplanet curve to?\n")

				star_index = int(input("Star index: "))
				print()

				print(">>> Running 'fit_exoplanet.py'")
				
				subprocess.run(["python3", "fit_exoplanet.py", 
								])

				print(f">>> Created csv for star {star_index}\n")
				print(f">>> Would you like to fit another star index?")
				print("(1) Yes!")
				print("(2) No :(")

				print()
				create_another = int(input("Answer: "))
				print()

				if create_another == 1:
					print(">>>  Another!")

				if create_another == 2:
					print(">>> Sad :(")
					skip_2_step20 = True
					break
				

		elif choice19 == 2:
			print()
			print(">>> Going to the next step...")
			print()
			skip_2_step20 = True
			break


		elif choice19 == 3:
			print()
			print(">>> Goodbye!")
			print()
			exit()
