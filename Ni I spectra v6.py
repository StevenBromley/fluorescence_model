# -*- coding: utf-8 -*-
"""
SJB
script for running fluorescence model and plotting synthetic spectra
01-13-2021
#Import the model and its dependencies:
"""
from fluor_v6 import *
"""
    fluorescence_model() function explanatio
    This function calculates fluorescence spectra from NIST data. The function returns a tuple; tuple elements are:
    Element 0: level populations, indexed the same as "raw_levs"
    Element 1: line wavelengths and intensities normalized to maximum intensity.
            
                                                USER INPUTS:            
    fluxes: Integrated flux in W/m^2; calculated with "grab_integrated_fluxes" function
    obj_vel: object heliocentric velocity in m/s; positive is defined AWAY from sun
    solar_dist: distance in meters from sun
    raw_lines: string array of NIST line data. Wavelengths must be in nm
    raw_levs: string array of NIST level data. Ensure both lines and levels have no headers. Energies must be in cm^-1
    ritz_col: column for ritz wavelengths (in VACUUM)
    lower_col: column of lower energy in lines array. Note 0 indexing, and energies must be in cm^-1
    upper_col: column of upper energy. Same limitations as 'lower_col'
    aval_col: column for A values in lines array.
    temp: optional; temperature of blackbody; assumed 5777K. Only used when radfield does not have flux at required wavelengths
    supp: optional; assumed True. Supplements flux array w/ blackbody for missing wavelengths. If unused, missing
    transition wavelengths in radfield are assigned rates of 0 for all radiative processes    
    renorm: optional; assumed True. Normalizes line intensities such that the maximum intensity = 1.       
"""
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
"""
For Kurucz solar spectrum provided, Col0: Wavelength [nm], Col1 Flux
The Kurucz MEASURED flux list is already in W / m^2 per nm
The Kurucz COMPUTED line lists have flux in units of ergs/cm**2/s/ster/nm; we need to convert to W / m^2 per nm
We want W/m**2/s/sr/nm. Note: 1 erg/cm**2 = 0.001 J/m**2
Therefore, we mutiply by 0.001 =  1e-3, and scale to solar distance and integrate over sky:
"""
#%%
#If using the COMPUTED solar spectra, keep this cell uncommented
radfield = np.genfromtxt('kurucz_solar_spectra.txt',delimiter ='\t',dtype=float) 
radfield1 = radfield.copy()
#the kurucz flux is in units of ergs / cm^2 per ster per nm; we convert to 
#J/m^2 per nm and scale it to 
radfield1[:,1] = np.pi * 1e-3 * (Rsun/solar_dist)**2 * radfield1[:,1] #Now in W/m^2 per nm
fluxes_comp = grab_integrated_fluxes(radfield1,raw_lines,lower_col,upper_col,obj_vel,solar_dist)
spectra1 = fluorescence_spectra(fluxes_comp,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp=5777,supp=True,renorm=False)
model_choice = spectra1[1] #the tuple element '1' is our synthetic spectra; tuple element 0 are the populations
#%%
#If using the MEASURED kurucz solar spectra, keep this cell uncommented
#W/m^2 /nm MEASURED kurucz spectrum from 300 to 1000 nm:
#This is at 1 AU by default; if scaled as (Rsun/solar_dist)**2, you get ~0.02 W/m^2 at 1AU (NOT ~1370 as expected)
radfield_meas = np.genfromtxt('kurucz_measured_300_1000_processed.txt',delimiter='\t') #2 col format
#Scale to solar distance:
radfield_meas[:,1] = (1.496e+11/solar_dist)**2 * radfield_meas[:,1]
fluxes_meas = grab_integrated_fluxes(radfield_meas,raw_lines,lower_col,upper_col,obj_vel,solar_dist)
spectra2 = fluorescence_spectra(fluxes_meas,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp=5777,supp=True,renorm=False)
model_choice = spectra2[1]
#%%
#If you are interested in seeing the difference in line intensities from computed vs measured solar spectra, uncomment
#the following:
measured = spectra1[1]
computed = spectra2[1] 
plt.clf()
plt.scatter(measured[:,2],computed[:,2])
plt.xlim(0.1 *measured[0,2], measured[-1,2]*1.1)
plt.ylim(0.1 *computed[0,2], computed[-1,2] * 1.1)
plt.grid(True)
plt.xlabel('Ni I lines using MEASURED Solar Spectrum + BLackbody')
plt.ylabel('Ni I lines using COMPUTED Solar Spectrum + Blackbody')
plt.rcParams.update({'font.size': 22})
#%%
#Normalize max intensity to 1:
max_int = np.amax(model_choice[:,2])
model_choice[:,2] = model_choice[:,2] / max_int
#%
plt.clf()
plt.figure(figsize=(16,12))
plt.rcParams.update({'font.size': 18})
plt.ylabel('Intensity [arb. units]')
plt.xlabel('Wavelength [nm]')
plt.stem(model_choice[:,0],model_choice[:,2], label = '{:} synthetic spectra'.format(species_string))
plt.legend()
plt.grid(True)
plt.xlim(100,500)
#Save a figure + datafile 
plt.savefig('{:}_spectra.pdf'.format(species_string))
np.savetxt('{:}_model_data.txt'.format(species_string),model_choice,delimiter ='\t')





