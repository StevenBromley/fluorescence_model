# -*- coding: utf-8 -*-
"""
Running Fluorescence Model for Ni I
SJB
05-26-2021
This uses version 13_2 of the model
"""
from fluor_v13_2 import *
import matplotlib.pyplot as plt
#
cwd = os.getcwd()
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
ritz_col = 1 #required user input
aval_col = 3 #required
lower_col = 5 #required
upper_col = 6 #required
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
#Import the Kurucz spectra provided; units W / m^2 per nm:
kurucz_flux = np.genfromtxt('kurucz_150nm-81um.txt',delimiter ='\t',dtype=float,skip_header=1) 
#Calculate the integrated fluxes; this takes into account a doppler-broadedning line profile at temp t_comet by default
fluxes = fluxes_with_profiles(kurucz_flux,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist)
model = fluorescence_spectra(fluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col)  
lines = model[1]
#The lines comes in 3 columns: col0 air wavelength, col1 vacuum wavelength, col2 intensity (normalized such that max=1 unless specified otherwise)
#%%
#plotting the spectra:
plt.clf()
plt.figure(figsize=(16,8))
plt.rcParams.update({'font.size': 22})
plt.stem(lines[:,0], lines[:,2], label = 'Ni I Fluorescence Model')
plt.xlabel("Wavelength (nm)")
plt.grid()
plt.ylabel("Intensity (arb. units)")
plt.xlim(200,600)

np.savetxt('v13_2_ni0_synthetic.txt',lines,delimiter ='\t')
