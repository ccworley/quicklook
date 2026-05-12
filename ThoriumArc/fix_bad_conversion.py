# Open the converted fits file, which has negative fluxes, and replace it with the original image from the raw data folder, which has the correct positive fluxes. This is necessary for the quicklook reduction to work properly, as it relies on the original image to perform the reduction steps. The converted fits file is only used for the initial inspection of the data, and is not used for the actual reduction process.

import os
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

UINT32_SIGN_OFFSET = 2**31

# Function to extract target names and filenames from a directory of FITS files
def get_files_types(raw_dir):
    files = []
    types = []
    orig_fname = []
    poor_conversion_dir = os.path.join(raw_dir, 'Poor_Conversion')
# For all files in raw_dir, open and extract target info
    for filename in os.listdir(poor_conversion_dir):
        #print(f"Processing file: {filename}")
        if filename.endswith('.fit') and '_converted' not in filename and '_test' not in filename and '_probe' not in filename:
            with fits.open(os.path.join(poor_conversion_dir, filename)) as hdul:
                header = hdul[0].header
                data = hdul[0].data.astype(np.int64)
                types
                data_min = np.min(data)
                data_max = np.max(data)
                print(data_min, data_max)
                new_header = header.copy()
                orig_bzero = new_header.get('BZERO')
                print(f"{filename} original BZERO: {orig_bzero}")

                if data_min >= 0:
                    print(f"{filename} does not have negative fluxes. Skipping conversion.")
                    continue

                corrected_data = data + UINT32_SIGN_OFFSET
                print(
                    f"{filename} corrected flux range: {np.min(corrected_data)} {np.max(corrected_data)}"
                )

                for key in ('BZERO', 'BSCALE'):
                    if key in new_header:
                        del new_header[key]

                newfilename = filename.replace('.fit', '_converted.fit')
                new_path = os.path.join(raw_dir, newfilename)
                new_hdu = fits.PrimaryHDU(data=corrected_data.astype(np.uint32), header=new_header)
                new_hdu.writeto(new_path, overwrite=True)
                print(f"Converted file saved as: {os.path.join(raw_dir, newfilename)}")

                with fits.open(new_path) as hdul:
                    chk_header = hdul[0].header
                    print(f"{filename} saved BZERO: {chk_header.get('BZERO')}")
                    chk_data = hdul[0].data.copy()
                    print(np.min(chk_data), np.max(chk_data))


    return files,types

# --------------------------------------------------------------
# Main
# --------------------------------------------------------------

# Get CCD32398 Thorium Arc files
ccd32398_raw_dir = '//home/astronomer/QuickLook_Data/Raw_Data/20260511/'
ccd32398_files, ccd32398_types = get_files_types(ccd32398_raw_dir)
ccd32398_thor_files = [f for f, t in zip(ccd32398_files, ccd32398_types) if 'THORIUM' in t.upper()]
print(f"CCD32398 - Thorium Arc Files: {ccd32398_thor_files}")

