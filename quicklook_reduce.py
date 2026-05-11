# SImple Reduction of HERCULES spectra
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import scipy as sp
import matplotlib.colors as colors
from scipy.optimize import curve_fit


# Function to extract target names and filenames from a directory of FITS files
def get_files_types(raw_dir):
    files = []
    types = []
# For all files in raw_dir, open and extract target info
    for filename in os.listdir(raw_dir):
        print(f"Processing file: {filename}")
        if filename.endswith('.fit'):
            with fits.open(os.path.join(raw_dir, filename)) as hdul:
                header = hdul[0].header
                target_name = header.get('OBJECT', 'Unknown')
                print(f"Target: {target_name}, Filename: {filename}")
                files.append(filename)
                types.append(target_name)
    return files,types



def plot_image(imdata, upc, lpc, type):
    plt.imshow(imdata,vmin=np.percentile(imdata, lpc), vmax=np.percentile(imdata, upc)) #,norm=colors.LogNorm(vmin=1000, vmax=4000))
    plt.title(type)
    plt.colorbar()
    plt.show()

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean)**2) / (2 * stddev**2))

def extract_order_peaks(white, background_level, nlevel, horzpix):
    vertical_slice = white[:, horzpix]
    vertical_slice_pixvector = np.arange(len(vertical_slice))
    order_pixel_values, order_info = sp.signal.find_peaks(vertical_slice, height=background_level * nlevel, distance=10)
    print('Order pixel values:', order_pixel_values)

    if len(order_pixel_values) == 0:
        raise RuntimeError('No order peaks found in the vertical slice. Try lowering nlevel.')
    
    return vertical_slice,order_pixel_values

def calculate_centroid(vertical_slice, centre_peakpix):
    centroid_half_window = 8
    ymin = max(0, centre_peakpix - centroid_half_window)
    ymax = min(vertical_slice.shape[0], centre_peakpix + centroid_half_window + 1)
    centroid_pixels = np.arange(ymin, ymax)
    centroid_flux = vertical_slice[centroid_pixels] - np.median(vertical_slice)
    centroid_flux = np.clip(centroid_flux, 0, None)

    if np.sum(centroid_flux) > 0:
        centre_ordpix = int(np.round(np.sum(centroid_pixels * centroid_flux) / np.sum(centroid_flux)))
    else:
        centre_ordpix = centre_peakpix

    print('Centre order pixel value:', centre_ordpix)
    
    return centroid_half_window,centre_ordpix

def calculate_fwhm(gaussian, white, horzpix, centre_ordpix):
    iprof = np.arange(centre_ordpix-10, centre_ordpix+10)
    ordprof_slice = white[iprof, horzpix]

    # Find full width at half maximum of the order profile, and overplot the nearest peak pixel and the refined centroid pixel.
    # Initial guess [amplitude, mean, stddev] helps convergence
    initial_guess = [max(ordprof_slice), np.mean(iprof), 1.0]
    popt, _ = curve_fit(gaussian, iprof, ordprof_slice, p0=initial_guess)
    amplitude_fit, mean_fit, sigma_fit = popt
    fwhm = np.floor(2 * np.sqrt(2 * np.log(2)) * abs(sigma_fit))
    #print(f"Fitted Sigma: {sigma_fit:.4f}")
    #print(f"Calculated FWHM: {fwhm:.4f}")
    
    return iprof,ordprof_slice,mean_fit,fwhm

def calculate_flux_ppix(ord_vpixvector, ord_hpixvector,ord_flux_ppsum, ord_flux_ppvar, ord_flux_ppstd, horzpix, centre_ordpix, iprof, ordprof_slice, mean_fit, fwhm, location):
    
    flux_fwhm = ordprof_slice[(iprof >= mean_fit - fwhm/2) & (iprof <= mean_fit + fwhm/2)]
    flux_sum_fwhm = np.sum(flux_fwhm)
    flux_variance_fwhm = np.var(flux_fwhm)
    flux_std_fwhm = np.std(flux_fwhm)
    #print('Flux across FWHM:', flux_fwhm)
    #print('Sum of flux across FWHM:', flux_sum_fwhm)
    #print('Variance of flux across FWHM:', flux_variance_fwhm)
    #print('Standard deviation of flux across FWHM:', flux_std_fwhm)

    if location =='right':
        print('Calculating flux for right side of order trace')
        ord_vpixvector.append(centre_ordpix)
        ord_hpixvector.append(horzpix)
        ord_flux_ppsum.append(flux_sum_fwhm)
        ord_flux_ppvar.append(flux_variance_fwhm)
        ord_flux_ppstd.append(flux_std_fwhm)
    elif location == 'left':
        print('Calculating flux for left side of order trace')
        ord_vpixvector.insert(0,centre_ordpix)
        ord_hpixvector.insert(0,horzpix)
        ord_flux_ppsum.insert(0,flux_sum_fwhm)
        ord_flux_ppvar.insert(0,flux_variance_fwhm)
        ord_flux_ppstd.insert(0,flux_std_fwhm)

# ----------------------------------------------------------------
# Main
#folder = '/home/astronomer/QuickLook_Data/Raw_Data/20260501/'
folder = 'test_data/20260501/'

files,types = get_files_types(folder)
print(types)

flat_files = [f for f, t in zip(files, types) if 'WHITE' in t.upper()]
thar_files = [f for f, t in zip(files, types) if 'THORIUM' in t.upper()]
dark_files = [f for f, t in zip(files, types) if 'DARK' in t.upper()]
star_files = [f for f, t in zip(files, types) if 'ALPHA' in t.upper()]
#print(star_files)

hdul = fits.open(folder + flat_files[0])
white = hdul[0].data 
hdul.close()

hdul = fits.open(folder + star_files[0])
star = hdul[0].data 
hdul.close()

hdul = fits.open(folder + thar_files[0])
thar = hdul[0].data 
hdul.close()

hdul = fits.open(folder + dark_files[0])
dark = hdul[0].data 
hdul.close()

print(folder)
print(flat_files[0], thar_files[0], dark_files[0],  star_files[0])


upc = 99
lpc = 10
#plot_image(white, upc, lpc, 'White Light Flat')
# plot_image(star, upc, lpc, 'Star')
# plot_image(thar, upc, lpc, 'Thorium Argon Lamp')
# plot_image(dark, upc, lpc, 'Dark')

background_level = np.median(white)
nlevel = 1.05
cen_hpix = 2250
max_fwhm = 20

# Initial horizontal pixel to slice through for order tracing. This is just a starting point, will need to loop through all columns and fit a polynomial to the order trace.
horzpix = cen_hpix

#-----------------------------------------------------------
# Trace order to the right of the initial horizontal pixel, then loop back to the left of the initial horizontal pixel, then fit a polynomial to the order trace. Save the vertical pixel locations of the order trace and the flux across the FWHM of the profile at each horizontal pixel, which will be used to weight the order trace fitting in the next step.
ord_vpixvector = []
ord_hpixvector = []
ord_flux_ppsum = []
ord_flux_ppvar = []
ord_flux_ppstd = []

# Find vertical slice and order peak locations
vertical_slice, order_pixel_values = extract_order_peaks(white, background_level, nlevel, horzpix)

# Pick the peak nearest the detector center, then refine to an intensity-weighted centroid.
detector_midpix = vertical_slice.shape[0] / 2.0
centre_peakpix = int(order_pixel_values[np.argmin(np.abs(order_pixel_values - detector_midpix))])

#Center the profile on the centroid of the order, which will be used to track the order position across the detector.
centroid_half_window, centre_ordpix = calculate_centroid(vertical_slice, centre_peakpix)

# Calculate the FWHM of the order profile, which will be used to weight the order trace fitting in the next step.
iprof, ordprof_slice, mean_fit, fwhm = calculate_fwhm(gaussian, white, horzpix, centre_ordpix)

# Overplot the 'y' pixel locations as lines on the column slice
fig = plt.figure(figsize=(12, 12))
plt.subplot(2, 1, 1)
# Overplot the order trace on the Dark-Subtract Raw image.
plt.imshow(white,vmin=1000,vmax=2000, aspect='auto')
for i in range(len(order_pixel_values)):
    plt.plot(horzpix, order_pixel_values[i], '*k', markersize=10)
plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('While Region with Locations of Order Peaks')

plt.subplot(2, 1, 2)
plt.plot(iprof,ordprof_slice)
plt.axvline(centre_peakpix, color='tab:orange', linestyle='--', linewidth=1.5)
plt.axvline(centre_ordpix, color='tab:red', linestyle='-', linewidth=1.5)
plt.plot(centre_peakpix, white[centre_peakpix, horzpix], 'o', color='tab:orange', markersize=6)
plt.plot(centre_ordpix, white[centre_ordpix, horzpix], 'o', color='tab:red', markersize=6)
plt.title('Order Profile Slice')
plt.xlabel('Y pixels of Column Slice')
plt.ylabel('Counts')
plt.text(centre_ordpix-10, ordprof_slice[10], 'Centre Order Pixel: VertColPix ' + str(centre_ordpix) + ' HorColPIx ' + str(horzpix), color='red', fontsize=12, ha='left', va='bottom')
plt.plot([mean_fit-fwhm/2,mean_fit+fwhm/2], [np.max(ordprof_slice)/2,np.max(ordprof_slice)/2], '--', color='tab:green', markersize=10, label='FWHM')
plt.legend(['Order profile', 'Nearest peak pixel', 'Refined centroid pixel'])
plt.show()

# Saved the flux across the FWHM of the profile and the variance of the flux across the FWHM, which will be used to weight the order trace fitting in the next step.
calculate_flux_ppix(ord_vpixvector, ord_hpixvector,ord_flux_ppsum, ord_flux_ppvar, ord_flux_ppstd, horzpix, centre_ordpix, iprof, ordprof_slice, mean_fit, fwhm,'right')

while horzpix < white.shape[1]-4:
    # Now got to the next horizontal pixel and repeat the process, then fit a polynomial to the order trace.
    horzpix = horzpix + 5
    print('Horizontal pixel:', horzpix, 'Max pixel:', white.shape[1])


    # Find vertical slice and order peak locations
    vertical_slice, order_pixel_values = extract_order_peaks(white, background_level, nlevel, horzpix)

    # Now find the peak nearest the previous order position, then refine to an intensity-weighted centroid.
    peakpixdist = np.abs(order_pixel_values - centre_ordpix)
    nearest_peakpix = int(order_pixel_values[np.argmin(peakpixdist)])
    #print('Nearest peak pixel value:', nearest_peakpix)
    #print('Previous centre order pixel value:', centre_ordpix)
    
    prev_center_ordpix = centre_ordpix
    
    #Center the profile on the centroid of the order, which will be used to track the order position across the detector.
    centroid_half_window, centre_ordpix = calculate_centroid(vertical_slice, nearest_peakpix)


    try:
        # Calculate the FWHM of the order profile, which will be used to weight the order trace fitting in the next step.
        iprof, ordprof_slice, mean_fit, fwhm = calculate_fwhm(gaussian, white, horzpix, centre_ordpix)
        # Saved the flux across the FWHM of the profile and the variance of the flux across the FWHM, which will be used to weight the order trace fitting in the next step.
        print(fwhm, horzpix, centre_ordpix)
        if fwhm < max_fwhm:
            calculate_flux_ppix(ord_vpixvector, ord_hpixvector,ord_flux_ppsum, ord_flux_ppvar, ord_flux_ppstd, horzpix, centre_ordpix, iprof, ordprof_slice, mean_fit, fwhm, 'right')
    except RuntimeError:
        print(f"Error occurred while calculating FWHM for horizontal pixel {horzpix}")
        continue
    
print('Horizontal pixel values:', ord_hpixvector)
print('Vertical pixel values:', ord_vpixvector)
print('Flux across FWHM:', ord_flux_ppsum)
print('Variance of flux across FWHM:', ord_flux_ppvar)
print('Standard deviation of flux across FWHM:', ord_flux_ppstd)

print(ord_hpixvector[0],ord_vpixvector[0])
# Overplot the 'y' pixel locations as lines on the column slice
fig = plt.figure(figsize=(12, 6))
# Overplot the order trace on the Dark-Subtract Raw image.
#plt.imshow(white[np.min(ord_vpixvector)-10:np.max(ord_vpixvector)+10,:],vmin=1000,vmax=2000, aspect='auto')
plt.imshow(white,vmin=np.percentile(white, 5),vmax=np.percentile(white, 95), aspect='auto')
for i in range(len(ord_vpixvector)):
    plt.plot(ord_hpixvector[i], ord_vpixvector[i], '*k', markersize=10)
plt.plot(ord_hpixvector[0],ord_vpixvector[0],  '*r', markersize=10)
#plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('While Region with Traced Central order')
#plt.ylim([np.max(ord_vpixvector)+10,np.min(ord_vpixvector)-10])
plt.colorbar()
plt.show()


#-----------------------------------------------------------
# Trace order to the LEFT of the initial horizontal pixel, then loop back to the left of the initial horizontal pixel, then fit a polynomial to the order trace. Save the vertical pixel locations of the order trace and the flux across the FWHM of the profile at each horizontal pixel, which will be used to weight the order trace fitting in the next step.
horzpix = cen_hpix

# Find vertical slice and order peak locations
vertical_slice, order_pixel_values = extract_order_peaks(white, background_level, nlevel, horzpix)

# Pick the peak nearest the detector center, then refine to an intensity-weighted centroid.
detector_midpix = vertical_slice.shape[0] / 2.0
centre_peakpix = int(order_pixel_values[np.argmin(np.abs(order_pixel_values - detector_midpix))])

#Center the profile on the centroid of the order, which will be used to track the order position across the detector.
centroid_half_window, centre_ordpix = calculate_centroid(vertical_slice, centre_peakpix)

# Calculate the FWHM of the order profile, which will be used to weight the order trace fitting in the next step.
iprof, ordprof_slice, mean_fit, fwhm = calculate_fwhm(gaussian, white, horzpix, centre_ordpix)

# Overplot the 'y' pixel locations as lines on the column slice
fig = plt.figure(figsize=(12, 12))
plt.subplot(2, 1, 1)
# Overplot the order trace on the Dark-Subtract Raw image.
plt.imshow(white,vmin=1000,vmax=2000, aspect='auto')
for i in range(len(order_pixel_values)):
    plt.plot(horzpix, order_pixel_values[i], '*k', markersize=10)
plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('While Region with Locations of Order Peaks')

plt.subplot(2, 1, 2)
plt.plot(iprof,ordprof_slice)
plt.axvline(centre_peakpix, color='tab:orange', linestyle='--', linewidth=1.5)
plt.axvline(centre_ordpix, color='tab:red', linestyle='-', linewidth=1.5)
plt.plot(centre_peakpix, white[centre_peakpix, horzpix], 'o', color='tab:orange', markersize=6)
plt.plot(centre_ordpix, white[centre_ordpix, horzpix], 'o', color='tab:red', markersize=6)
plt.title('Order Profile Slice')
plt.xlabel('Y pixels of Column Slice')
plt.ylabel('Counts')
plt.text(centre_ordpix-10, ordprof_slice[10], 'Centre Order Pixel: VertColPix ' + str(centre_ordpix) + ' HorColPIx ' + str(horzpix), color='red', fontsize=12, ha='left', va='bottom')
plt.plot([mean_fit-fwhm/2,mean_fit+fwhm/2], [np.max(ordprof_slice)/2,np.max(ordprof_slice)/2], '--', color='tab:green', markersize=10, label='FWHM')
plt.legend(['Order profile', 'Nearest peak pixel', 'Refined centroid pixel'])
plt.show()

# Saved the flux across the FWHM of the profile and the variance of the flux across the FWHM, which will be used to weight the order trace fitting in the next step.
#calculate_flux_ppix(ord_vpixvector, ord_hpixvector,ord_flux_ppsum, ord_flux_ppvar, ord_flux_ppstd, horzpix, centre_ordpix, iprof, ordprof_slice, mean_fit, fwhm, 'duplicate')

while horzpix > 4:
    # Now got to the next horizontal pixel and repeat the process, then fit a polynomial to the order trace.
    horzpix = horzpix - 5
    print('Horizontal pixel:', horzpix, 'Max pixel:', white.shape[1])


    # Find vertical slice and order peak locations
    vertical_slice, order_pixel_values = extract_order_peaks(white, background_level, nlevel, horzpix)

    # Now find the peak nearest the previous order position, then refine to an intensity-weighted centroid.
    peakpixdist = np.abs(order_pixel_values - centre_ordpix)
    nearest_peakpix = int(order_pixel_values[np.argmin(peakpixdist)])
    #print('Nearest peak pixel value:', nearest_peakpix)
    #print('Previous centre order pixel value:', centre_ordpix)
    
    prev_center_ordpix = centre_ordpix
    
    #Center the profile on the centroid of the order, which will be used to track the order position across the detector.
    centroid_half_window, centre_ordpix = calculate_centroid(vertical_slice, nearest_peakpix)


    try:
        # Calculate the FWHM of the order profile, which will be used to weight the order trace fitting in the next step.
        iprof, ordprof_slice, mean_fit, fwhm = calculate_fwhm(gaussian, white, horzpix, centre_ordpix)
        
        # plt.plot(iprof,ordprof_slice)
        # plt.axvline(nearest_peakpix, color='tab:orange', linestyle='--', linewidth=1.5)
        # plt.axvline(centre_ordpix, color='tab:red', linestyle='-', linewidth=1.5)
        # plt.plot(nearest_peakpix, white[nearest_peakpix, horzpix], 'o', color='tab:orange', markersize=6)
        # plt.plot(centre_ordpix, white[centre_ordpix, horzpix], 'o', color='tab:red', markersize=6)
        # plt.title('Order Profile Slice')
        # plt.xlabel('Y pixels of Column Slice')
        # plt.ylabel('Counts')
        # plt.text(centre_ordpix-10, ordprof_slice[10], 'Centre Order Pixel: VertColPix ' + str(centre_ordpix) + ' HorColPIx ' + str(horzpix), color='red', fontsize=12, ha='left', va='bottom')
        # plt.plot([mean_fit-fwhm/2,mean_fit+fwhm/2], [np.max(ordprof_slice)/2,np.max(ordprof_slice)/2], '--', color='tab:green', markersize=10, label='FWHM')
        # plt.legend(['Order profile', 'Nearest peak pixel', 'Refined centroid pixel'])
        # plt.show()
                
        
        # Saved the flux across the FWHM of the profile and the variance of the flux across the FWHM, which will be used to weight the order trace fitting in the next step.
        print(fwhm, horzpix, centre_ordpix)
        if fwhm < max_fwhm:
            calculate_flux_ppix(ord_vpixvector, ord_hpixvector,ord_flux_ppsum, ord_flux_ppvar, ord_flux_ppstd, horzpix, centre_ordpix, iprof, ordprof_slice, mean_fit, fwhm, 'left')
    except RuntimeError:
        print(f"Error occurred while calculating FWHM for horizontal pixel {horzpix}")
        continue
    
    
print('Horizontal pixel values:', ord_hpixvector)
print('Vertical pixel values:', ord_vpixvector)
print('Flux across FWHM:', ord_flux_ppsum)
print('Variance of flux across FWHM:', ord_flux_ppvar)
print('Standard deviation of flux across FWHM:', ord_flux_ppstd)

print(ord_hpixvector[0],ord_vpixvector[0])
# Overplot the 'y' pixel locations as lines on the column slice
fig = plt.figure(figsize=(12, 12))
plt.subplot(131)
# Overplot the order trace on the Dark-Subtract Raw image.
#plt.imshow(white[np.min(ord_vpixvector)-10:np.max(ord_vpixvector)+10,:],vmin=1000,vmax=2000, aspect='auto')
plt.imshow(white,vmin=np.percentile(white, 5),vmax=np.percentile(white, 95), aspect='auto')
for i in range(len(ord_vpixvector)):
    plt.plot(ord_hpixvector[i], ord_vpixvector[i], '*k', markersize=10)
plt.plot(ord_hpixvector[0],ord_vpixvector[0],  '*r', markersize=10)
#plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('White Region with Traced Central order')
#plt.ylim([np.max(ord_vpixvector)+10,np.min(ord_vpixvector)-10])
plt.colorbar()

plt.subplot(132)
# Overplot the order trace on the Dark-Subtract Raw image.
plt.imshow(thar,vmin=np.percentile(thar, 5),vmax=np.percentile(thar, 95), aspect='auto')
for i in range(len(ord_vpixvector)):
    plt.plot(ord_hpixvector[i], ord_vpixvector[i], '*k', markersize=10)
plt.plot(ord_hpixvector[0],ord_vpixvector[0],  '*r', markersize=10)
#plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('Thar Region with Traced Central order')
#plt.ylim([np.max(ord_vpixvector)+10,np.min(ord_vpixvector)-10])
plt.colorbar()

plt.subplot(133)
# Overplot the order trace on the Dark-Subtract Raw image.
plt.imshow(star,vmin=np.percentile(star, 5),vmax=np.percentile(star, 95), aspect='auto')
for i in range(len(ord_vpixvector)):
    plt.plot(ord_hpixvector[i], ord_vpixvector[i], '*k', markersize=10)
plt.plot(ord_hpixvector[0],ord_vpixvector[0],  '*r', markersize=10)
#plt.plot(horzpix, centre_ordpix, '*r', markersize=10)
plt.title('Star Region with Traced Central order')
#plt.ylim([np.max(ord_vpixvector)+10,np.min(ord_vpixvector)-10])
plt.colorbar()
plt.show()


