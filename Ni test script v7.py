# -*- coding: utf-8 -*-
"""
Ni test script for fluor_v7
01-27-2021
SJB
"""
from fluor_v7 import *

#%%
cwd = os.getcwd()
species_string = 'Ni_I' #required for outputting plot and datafile at end.
ritz_col = 1 #required user input
lower_col = 5 #required
upper_col = 6 #required
aval_col = 3 #required
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
Rsun = 6.957e8 # meters
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
m_species= 9.746267511e-26 #mass of Nickel in kg
t_comet = 100 # Kelvin; temperature used for doppler profile.
#Import the Kurucz spectra provided; units W / m^2 per nm:
kurucz_flux = np.genfromtxt('kurucz_150nm-81um.txt',delimiter ='\t',dtype=float,skip_header=1) 
#
fluxes = fluxes_with_profiles(kurucz_flux,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist,t_comet, m_species)
model = fluorescence_spectra(fluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col)  
pops = model[0]
lines = model[1]
#
plt.clf()
plt.stem(lines[:,0], lines[:,2], label = 'Synthetic Ni I Spectra')
plt.grid(True)
plt.legend()
plt.xlim(0,1000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (arb. units)')