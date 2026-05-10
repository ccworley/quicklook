# Read in thorium files (HERCULES@pyreduce)

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.io import loadmat
import h5py

# Define the directory containing the thorium files
directory = 'ThoriumArc/'

meg_thar = 'megara_thardat_o.mat'

with h5py.File('megara_thardat_o.mat', 'r') as f:
    # Get all variable names
    variable_names = list(f.keys())
    
    print("Variables in file:", variable_names)
    # Access a variable
    order = f['order'][:]
    wv_air = f['air'][:]
    wv_vac = f['vac'][:]
    uc = f['uc'][:]
    ul = f['ul'][:]
    ur = f['ur'][:]
    vc = f['vc'][:]
    vl = f['vl'][:]
    vr = f['vr'][:]
    wavnum = f['wavnum'][:]
    spec = f['species'][:]


print(order)
print(wv_air)
print(wv_vac)
print(uc)
print(ul)
print(ur)
print(vc)
print(vl)
print(vr)
print(wavnum)
print(spec) 

plt.figure(figsize=(10, 6))
plt.scatter(wv_air, order, label='Air Wavelength', color='blue')
plt.scatter(wv_vac, order, label='Vacuum Wavelength', color='orange')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Order')
plt.title('Wavelength vs Order for HERCULES@pyreduce')
plt.legend()
plt.grid()
plt.show()

ordlist = list(set(order.tolist()))
print(ordlist)