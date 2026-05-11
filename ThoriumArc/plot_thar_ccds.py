# Compare the thorium arc exposures from CCD486 and CCD32398 (prior to use of IR filter)
import os
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

# Function to extract target names and filenames from a directory of FITS files
def get_files_types(raw_dir):
    files = []
    types = []
# For all files in raw_dir, open and extract target info
    for filename in os.listdir(raw_dir):
        print(f"Processing file: {filename}")
        if filename.endswith('.fit'):
            print('here')
            with fits.open(os.path.join(raw_dir, filename)) as hdul:
                header = hdul[0].header
                target_name = header.get('OBJECT', 'Unknown')
                #print(f"CCD4 - Target: {target_name}, Filename: {filename}")
                files.append(filename)
                types.append(target_name)
    return files,types

# --------------------------------------------------------------
# Main
# --------------------------------------------------------------

# Get CCD486 Thorium Arc files
ccd486_raw_dir = '/home/astronomer/QuickLook_Data/Raw_Data/20260112'
ccd486_files, ccd486_types = get_files_types(ccd486_raw_dir)
ccd486_thor_files = [f for f, t in zip(ccd486_files, ccd486_types) if 'THORIUM' in t.upper()]
print(f"CCD486 - Thorium Arc Files: {ccd486_thor_files}")

# Get CCD32398 Thorium Arc files
#ccd32398_raw_dir = '//home/astronomer/QuickLook_Data/Raw_Data/20260501/'
ccd32398_raw_dir = '//home/astronomer/QuickLook_Data/Raw_Data/20260511/'
ccd32398_files, ccd32398_types = get_files_types(ccd32398_raw_dir)
ccd32398_thor_files = [f for f, t in zip(ccd32398_files, ccd32398_types) if 'THORIUM' in t.upper()]
print(f"CCD32398 - Thorium Arc Files: {ccd32398_thor_files}")

# Plot first Thorium Arc from CCD486 and from CCD32398 side by side
print(f"CCD486 - First Thorium Arc File: {ccd486_thor_files[0]}")
print(f"CCD32398 - First Thorium Arc File: {ccd32398_thor_files[0]}")


# Load the first Thorium Arc file from each CCD
with fits.open(os.path.join(ccd486_raw_dir, ccd486_thor_files[1])) as hdul:
    thar_ccd486_bkg = np.median(hdul[0].data)
    print("CCD486 - Median (Bkg): ", thar_ccd486_bkg)
    thar_ccd486 = hdul[0].data# - np.median(hdul[0].data)  # Subtract bias level to make it more comparable to CCD32398 (which has bias already subtracted)

with fits.open(os.path.join(ccd32398_raw_dir, ccd32398_thor_files[1])) as hdul:
    thar_ccd32398_bkg = np.median(hdul[0].data)
    print("CCD32398 - Median (Bkg): ", thar_ccd32398_bkg)
    thar_ccd32398 = hdul[0].data# - np.median(hdul[0].data)  # Subtract bias level to make it more comparable to CCD486 (which has bias already subtracted)
    
print(thar_ccd486.shape)  # Note: CCD486 (4096,4130) is (y,x) in image plot, x are pixels corresponding to order location, y are pixels indicating wavelength location
print(thar_ccd32398.shape) # Note: CCD32398 (4108,4096) is (y,x) in image plot, y are pixels corresponding to order location, x are pixels indicating wavelength location

thar_ccd486_flat = []
thar_ccd32398_flat = []

for i in range(thar_ccd486.shape[1]):  # Want to sample along order, therefore in x (CCD486 is rotated)
    maxpix_col = np.nanmean(thar_ccd486[:, i])
    thar_ccd486_flat.append(maxpix_col)

for i in range(thar_ccd32398.shape[0]): # Want to sample along order, therefore in y (CCD32398 is not rotated)
    maxpix_col = np.nanmean(thar_ccd32398[i, :])
    thar_ccd32398_flat.append(maxpix_col)


print(len(thar_ccd486_flat), type(thar_ccd486_flat))
print(len(thar_ccd32398_flat), type(thar_ccd32398_flat))

uppc = 99
lopc = 10

do_plt1 = True #False #

if do_plt1:
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))

    # Plot histogram for CCD486
    im1 = ax1.bar(range(len(thar_ccd486_flat)), thar_ccd486_flat)
    ax1.plot([0,thar_ccd486.shape[1]], [thar_ccd486_bkg, thar_ccd486_bkg], 'r--', label='95th Percentiles')
    ax1.set_title('CCD486 - Thorium Argon Histogram')
    ax1.set_xlabel('Pixel Value (N={})'.format(len(thar_ccd486_flat)))
    ax1.set_ylabel('Mean Column Counts - Along Order')
    ax1.set_yscale('log')  # Set y-axis to logarithmic scale for better visibility of low-frequency bins

    # Plot histogram for CCD32398
    im2 = ax2.bar(range(len(thar_ccd32398_flat)), thar_ccd32398_flat)
    ax2.plot([0,thar_ccd32398.shape[0]], [thar_ccd32398_bkg, thar_ccd32398_bkg], 'r--', label='95th Percentiles')
    ax2.set_title('CCD32398 - Thorium Argon Histogram')
    ax2.set_xlabel('Pixel Value (N={})'.format(len(thar_ccd32398_flat)))
    ax2.set_ylabel('Mean Column Counts - Along Order')
    ax2.set_yscale('log')  # Set y-axis to logarithmic scale for better visibility of low-frequency bins

    # Plot Thorium Arc from CCD486
    im3 = ax3.imshow(thar_ccd486, vmin=np.percentile(thar_ccd486, lopc), vmax=np.percentile(thar_ccd486, uppc))
    ax3.set_title('CCD486 - Thorium Argon Lamp')
    xpixel_range = np.arange(thar_ccd486.shape[1])
    ypixel_range = np.arange(thar_ccd486.shape[0])
    ax3.set_xlabel('X Pixel (N={})'.format(thar_ccd486.shape[1]))
    ax3.set_ylabel('Y Pixel (N={})'.format(thar_ccd486.shape[0]))
    plt.colorbar(im3, ax=ax3)

    # Plot Thorium Arc from CCD32398
    im4 = ax4.imshow(thar_ccd32398, vmin=np.percentile(thar_ccd32398, lopc), vmax=np.percentile(thar_ccd32398, uppc))
    ax4.set_title('CCD32398 - Thorium Argon Lamp')
    ax4.set_xlabel('X Pixel (N={})'.format(thar_ccd32398.shape[1]))
    ax4.set_ylabel('Y Pixel (N={})'.format(thar_ccd32398.shape[0]))
    plt.colorbar(im4, ax=ax4)

    plt.tight_layout()
    plt.show()


# Zoom in on the brightest end of the orders and overlap to see if can line up the pixel columns
min_pix = 0
max_pix = 4110

fig, ((ax1, ax2)) = plt.subplots(2, 1, figsize=(16, 12))   #, (ax3, ax4)

# Plot histogram for CCD486
im1 = ax1.bar(range(len(thar_ccd486_flat)), thar_ccd486_flat)
ax1.plot([0,thar_ccd486.shape[1]], [thar_ccd486_bkg, thar_ccd486_bkg], 'r--', label='95th Percentiles')
ax1.set_title('CCD486 - Thorium Argon Histogram')
ax1.set_xlabel('Pixel Value (N={})'.format(len(thar_ccd486_flat)))
ax1.set_ylabel('Mean Column Counts - Along Order')
ax1.set_yscale('log')  # Set y-axis to logarithmic scale for better visibility of low-frequency bins
ax1.set_xlim(min_pix, max_pix)  # Zoom in on the brightest end of the orders

# Plot histogram for CCD32398
im2 = ax2.bar(range(len(thar_ccd32398_flat)), thar_ccd32398_flat)
ax2.plot([0,thar_ccd32398.shape[0]], [thar_ccd32398_bkg, thar_ccd32398_bkg], 'r--', label='95th Percentiles')
ax2.set_title('CCD32398 - Thorium Argon Histogram')
ax2.set_xlabel('Pixel Value (N={})'.format(len(thar_ccd32398_flat)))
ax2.set_ylabel('Mean Column Counts - Along Order')
ax2.set_yscale('log')  # Set y-axis to logarithmic scale for better visibility of low-frequency bins
ax2.set_xlim(min_pix, max_pix)  # Zoom in on the brightest end of the orders

plt.tight_layout()
plt.show()