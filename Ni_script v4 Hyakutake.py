# -*- coding: utf-8 -*-
"""
SJB
script for running fluorescence model and plotting synthetic spectra
01-07-2021

To use this code, do the following:
1. Download lines from NIST in tab delimited format
    -- Only grab the lines with A values
    -- Remove any header lines (NOTE: They may also be at the lower end of the line list, as transitions above 2000nm
    are observed in vacuum)
2. Open the lines file in excel and re-save; this process removes the extra "" around all data
3. Do the same for levels file. 
    -- LEVEL ENERGIES SHOULD BE IN INVERSE CENTIMETERS.
4. Download Kurucz high-res solar spectrum, or use whichever you prefer. The code is setup so that the radfield passed
to the fluorescence model is in W/m^3 (W/m^2 per wavelength (meters)) 
5. Adjust all parameters in the "User Parameters" section as necessary
6. Run the full script. 

-------------------------
The tuple returned from the fluorescence_spectra function contains two elements:
    (1) level populations, and
    (2) line wavelengths and intensities in 3 column format. Col0 wavelengths in air, Col1 wavelengths in vacuum, col2 intensity (arb)

"""
#Import the model and its dependencies:
from fluor_v4 import *
"""
    fluorescence_model() function:
                USER INPUTS:           
    radfield: flux per wavelength interval array; 2 col (wavelength, flux) in W/m^2 per nm
    obj_vel: object heliocentric velocity in m/s; positive is defined AWAY from sun
    solar_dist: distance in meters from sun
    raw_lines: string array of NIST line data. Wavelengths must be in nm
    raw_levs: string array of NIST level data. Ensure both lines and levels have no headers. Energies must be in cm^-1
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
hyak0 = np.genfromtxt('offset_0_arcsec.tab',usecols = (0,1))
#This cell of code can be run (if in Spyder) to verify that column indexes are set correctly before running model further down
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
radfield = np.loadtxt('kurucz_solar_spectra.txt',delimiter ='\t',dtype=float) 
#For Kurucz solar spectrum provided, Col0: Wavelength [nm], Col1 Flux [ergs/cm**2/s/ster/nm] at solar surface.
"""
The Kurucz line lists have flux in units of ergs/cm**2/s/ster/nm
We want W/m**2/s/sr/meter. Note: 1 erg/cm**2 = 0.001 J/m**2
We will avoid going to frequency form as the conversion is not so straightforward, and easier to work in wavelength.
Therefore, we mutiply by 0.001 =  1e-3 
"""
radfield1 = radfield.copy()
radfield1[:,1] = 4 * np.pi * 1e-3 * (Rsun/solar_dist)**2 * radfield1[:,1] #Now in W/m^2 per nm
#%%
synthetic_spectra = fluorescence_spectra(radfield1,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col)
pops = synthetic_spectra[0]
model_lines = synthetic_spectra[1]
#%%
hy_max = np.amax(hyak0[:,1])
hyak0[:,1] = hyak0[:,1]/hy_max

#Let's try using the 345.86 line to normalize
norm_line = 345.86 #nm
norm_model_idx = find_nearest_index(model_lines[:,0],norm_line)
norm_comet_idx = find_nearest_index(hyak0[:,0]/10,norm_line)
norm_int_model = model_lines[norm_model_idx,2]
norm_int_comet = hyak0[norm_comet_idx,1]
#%
model_lines[:,2] = model_lines[:,2] * norm_int_comet / norm_int_model
plt.clf()
plt.figure(figsize=(16,12))
plt.rcParams.update({'font.size': 18})
plt.ylabel('Intensity [arb. units]')
plt.xlabel('Wavelength [nm]')
plt.stem(model_lines[:,0],model_lines[:,2], label = '{:} synthetic spectra'.format(species_string))
plt.plot(hyak0[:,0]/10,hyak0[:,1],color='green',label='Hyakutake 0 offset')
plt.legend()
plt.grid(True)
plt.xlim(100,500)
plt.savefig('{:}_spectra_Hyak.pdf'.format(species_string))
np.savetxt('{:}_data_Hyak.txt'.format(species_string),model_lines,delimiter ='\t')





