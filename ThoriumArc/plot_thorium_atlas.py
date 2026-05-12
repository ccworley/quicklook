# Read in thorium files (HERCULES@pyreduce)
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.io import loadmat
import h5py


# pyreduce files
herc_fits = 'hercules_all.thar_master.fits'
herc_npz = 'hercules_all.ord_default.npz'

# Load the first Thorium Arc file from each CCD
with fits.open(herc_fits) as hdul:
    herc_thar = hdul[0].data
    print("Hercs fits table keys: ", list(hdul[0].header.keys()))
    
with np.load(herc_npz) as data:
    # get variable names
    variable_names = data.files
    print("Hercs npz variable names: ", variable_names)
    orders = data['orders']
    columns = data['column_range']
    
print(orders)
print(columns)

# # Won't load --- Define the directory containing the thorium files
# directory = '/home/astronomer/QuickLook_Data/Raw_Data/20260112/'

# meg_thar = 'ThAr_info_blue_500.mat'

# with h5py.File(directory + meg_thar, 'r') as f:
#     # Get all variable names
#     variable_names = list(f.keys())
    
#     print("Variables in file:", variable_names)




# # Define the directory containing the thorium files
# directory = 'ThoriumArc/'

# meg_thar = 'megara_thardat_o.mat'

# with h5py.File(directory + meg_thar, 'r') as f:
#     # Get all variable names
#     variable_names = list(f.keys())
    
#     print("Variables in file:", variable_names)
#     # Access a variable
#     order = f['order'][:]
#     wv_air = f['air'][:]
#     wv_vac = f['vac'][:]
#     uc = f['uc'][:]
#     ul = f['ul'][:]
#     ur = f['ur'][:]
#     vc = f['vc'][:]
#     vl = f['vl'][:]
#     vr = f['vr'][:]
#     wavnum = f['wavnum'][:]
#     spec = f['species'][:]


# print(order)
# print(wv_air)
# print(wv_vac)
# print(uc)
# print(ul)
# print(ur)
# print(vc)
# print(vl)
# print(vr)
# print(wavnum)
# print(spec) 

# plt.figure(figsize=(10, 6))
# plt.scatter(wv_air, order, label='Air Wavelength', color='blue')
# plt.scatter(wv_vac, order, label='Vacuum Wavelength', color='orange')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Order')
# plt.title('Wavelength vs Order for HERCULES@pyreduce')
# plt.legend()
# plt.grid()
# plt.show()

# ordlist = list(set(order.tolist()))
# print(ordlist)


