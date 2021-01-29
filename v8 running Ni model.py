# -*- coding: utf-8 -*-
"""
Running Fluorescence Model for Ni I
SJB
01-28-2021
This uses version 7 of the model
"""
from fluor_v8 import *
"""
fluor_v7 contains a number of dependencies and other useful packages:
    numpy (np), matplotlib.pyplot (plt), os, math
    Numerical integratin is handled internally using a simple rectangle code I have written. 
"""
cwd = os.getcwd()
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
species_string = 'Ni_I' #required for outputting plot and datafile at end.
ritz_col = 1 #required user input
lower_col = 5 #required
upper_col = 6 #required
aval_col = 3 #required
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
m_species= 9.746267511e-26 #mass of Nickel in kg
t_comet = 100 # Kelvin; temperature used for doppler profile.
#Import the Kurucz spectra provided; units W / m^2 per nm:
kurucz_flux = np.genfromtxt('kurucz_150nm-81um.txt',delimiter ='\t',dtype=float,skip_header=1) 
#Calculate the integrated fluxes; this takes into account a doppler-broadedning line profile at temp t_comet by default
fluxes = fluxes_with_profiles(kurucz_flux,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist,t_comet, m_species)
model = fluorescence_spectra(fluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col)  
lines = model[1]
#The lines comes in 3 columns: col0 air wavelength, col1 vacuum wavelength, col2 intensity (normalized such that max=1 unless specified otherwise)

#plotting the spectra:
plt.clf()
plt.stem(lines[:,1], lines[:,2], label = 'Synthetic {:} spectra'.format(species_string))
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity (arb. units)")
plt.xlim(0,1000)

np.savetxt('Ni_I_synthetic.txt',lines,delimiter='\t')
