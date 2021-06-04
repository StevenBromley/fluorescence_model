# -*- coding: utf-8 -*-
"""
Running Fluorescence Model for Ni I
SJB
05-26-2021
This uses version 13_2 of the model
"""
from fluor_v13_3 import *
#
cwd = os.getcwd()
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
ritz_col = 1 #required user input
lower_col = 5 #required
upper_col = 6 #required
aval_col = 3 #required
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
#Import the Kurucz spectra provided; units W / m^2 per nm:

rad_uv = np.genfromtxt('kurucz_vac_comp150-200.txt', delimiter='\t')
rad_vis = np.genfromtxt('kurucz_air_comp200-300_meas300-1000.txt', delimiter='\t')
rad_ir = np.genfromtxt('kurucz_vac_comp1-81um.txt', delimiter='\t')

#Calculate the integrated fluxes; this takes into account a doppler-broadedning line profile at temp t_comet by default
nifluxes = fluxes_with_profiles(rad_uv,rad_vis,rad_ir,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist)
nimodel = fluorescence_spectra(nifluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col)  
nilines = nimodel[1]
#The lines comes in 3 columns: col0 air wavelength, col1 vacuum wavelength, col2 intensity (normalized such that max=1 unless specified otherwise)
#%%
#plotting the spectra:
plt.clf()
plt.figure(figsize=(16,8))
plt.rcParams.update({'font.size': 22})
plt.stem(nilines[:,0], nilines[:,2], label = 'Ni I Fluorescence Model')
plt.xlabel("Wavelength (nm)")
plt.grid()
plt.ylabel("Intensity (arb. units)")
plt.xlim(200,600)

np.savetxt('v13_3_ni0_synthetic.txt',nilines,delimiter ='\t')
