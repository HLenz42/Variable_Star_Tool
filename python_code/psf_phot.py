import os
import subprocess

# === USER-SPECIFIED INPUT ===
catalog_path = "/home/heather/Projects/PHYS599/python_code/psf_phot/NGC7332_megacam_u_combined_image_regionABCDEFGH.cat"
output_dir = "/home/heather/Projects/PHYS599/python_code/psf_phot"
psfex_config = "/home/heather/Projects/PHYS599/python_code/psf_phot/default.psfex"
# =============================

def run_psfex(catalog_path, psfex_config, output_dir):
    """Run PSFEx on the SExtractor catalog to generate the PSF model."""
    if not os.path.exists(catalog_path):
        raise FileNotFoundError(f"Catalog file not found at {catalog_path}")
    
    cmd = [
        'psfex', catalog_path,
        '-c', psfex_config,
        '-OUTCAT_NAME', os.path.join(output_dir, 'psfex_output.cat'),
        '-PSF_DIR', output_dir
    ]
    print(f"Running PSFEx on {catalog_path}...")
    try:
        subprocess.run(cmd, check=True)
        print(f"PSFEx output saved in {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"PSFEx failed: {e}")
        raise

def main():
    print(f"Checking catalog at {catalog_path}...")
    if os.path.exists(catalog_path):
        print(f"Catalog file size: {os.path.getsize(catalog_path)} bytes")
        print("Catalog header and first 3 data rows:")
        with open(catalog_path, 'r') as f:
            data_rows_printed = 0
            for i, line in enumerate(f):
                if line.startswith('#'):
                    print(line.strip())
                else:
                    if data_rows_printed < 3:
                        print(line.strip())
                        values = line.split()
                        if len(values) >= 9:
                            vignet_values = len(values[8:])  # VIGNET is column 9
                            print(f"VIGNET column values: {vignet_values}")
                        data_rows_printed += 1
                    else:
                        break
    else:
        raise FileNotFoundError(f"Catalog file not found at {catalog_path}")

    run_psfex(catalog_path, psfex_config, output_dir)
    print(f"\nPSF model saved in: {output_dir}")

if __name__ == "__main__":
    main()

def run_sextractor(image_path, sex_config, catalog_path):
    """Run Source Extractor on a FITS image to generate a catalog."""
    cmd = [
        'sex', image_path,
        '-c', sex_config,
        '-CATALOG_NAME', catalog_path
    ]
    print(f"Running SExtractor on {image_path}...")
    try:
        subprocess.run(cmd, check=True)
        if not os.path.exists(catalog_path):
            raise FileNotFoundError(f"SExtractor did not create catalog at {catalog_path}")
        os.chmod(catalog_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
        print(f"Catalog created at {catalog_path}")
        print(f"Catalog file size: {os.path.getsize(catalog_path)} bytes")
        print(f"Catalog permissions: {oct(os.stat(catalog_path).st_mode)[-3:]}")
        # Read FITS_LDAC catalog
        try:
            with fits.open(catalog_path) as hdul:
                # FITS_LDAC typically has the catalog in HDU 2
                if len(hdul) < 2:
                    raise ValueError("FITS_LDAC catalog does not contain expected data HDU")
                catalog = hdul[1].data
                columns = hdul[1].columns.names
                print("Catalog columns:", columns)
                # Verify required columns
                required = ['X_IMAGE', 'Y_IMAGE', 'FLUX_RADIUS', 'FLUX_AUTO', 'FLUXERR_AUTO', 
                           'MAG_AUTO', 'MAGERR_AUTO', 'FLAGS', 'VIGNET']
                missing = [col for col in required if col not in columns]
                if missing:
                    print(f"Warning: Missing columns: {missing}")
                # Check VIGNET shape
                if 'VIGNET' in columns:
                    vignet_shape = catalog['VIGNET'].shape
                    print(f"VIGNET data shape: {vignet_shape}")
                    if len(vignet_shape) == 3 and vignet_shape[1:] == (35, 35):
                        print("VIGNET size is correct: 35x35")
                    else:
                        print(f"Warning: VIGNET size is {vignet_shape[1:]}")
                # Print first few rows
                print("First 3 catalog rows (selected columns):")
                for i in range(min(3, len(catalog))):
                    row = {col: catalog[col][i] for col in columns if col != 'VIGNET'}
                    print(row)
        except Exception as e:
            print(f"Error reading FITS catalog: {e}")
    except subprocess.CalledProcessError as e:
        print(f"SExtractor failed: {e}")
        raise
    except FileNotFoundError as e:
        print(e)
        raise

def main():
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(input_image))[0]
    catalog_name = f"{base_name}.cat"
    catalog_path = os.path.join(output_dir, catalog_name)
    run_sextractor(input_image, sex_config, catalog_path)
    print(f"\nCatalog saved at: {catalog_path}")

if __name__ == "__main__":
    main()