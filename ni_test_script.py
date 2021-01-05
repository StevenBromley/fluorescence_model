# -*- coding: utf-8 -*-
"""
SJB
script for running fluorescence model and plotting synthetic spectra
01-04-2021

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

The code is organized as follows:
  
1. Load in data / parameters for model run.
2. Process the lines/levels files, i.e. match lines to the level indexing scheme.
3. Generate all stimulated emission / absorption rates from A values
4. Populate LHS and RHS of matrix equation.
5. Solve matrix equation for populations.
6. Calculate line intensities.
7. Save outputs in two ways: (1) a figure showing normalized intensity vs wavelength, and a 2col text file containing
the wavelength vs intensity data

The tuple returned from the fluorescence_spectra function contains two elements:
    (1) level populations, normalized to the total population, and
    (2) line wavelengths and intensities in 2 column fo

"""
#Import the model and its dependencies:
from fluor_model_func_v1 import *

"""
                USER INPUTS
ritz_col: column (note: python is 0 indexed) of ritz wavelengths in NIST lines file
lower_col: column for lower level ENERGY (in inverse cm)
upper_col: column for upper level ENERGY (in inverse cm)
aval_col: column for Einstein A values
obj_vel: comet/object heliocentric velocity in meters for person (-36.7 x 1e3 [m/s] for Hyakutake)
raw_lines: load in NIST lines data as strings
raw_levels: load in NIST levels data as strings. 
solar_dist: distance from sun in meters; assumed 
radfield: high resolution Kurucz solar spectrum
Supp: if 'True', supplements the imported solar spectra by using a blackbody intensity value if a given lines wavelength
is outside the range of values in the imported solar spectrum
                Inputs unlikely to change:
        temp: blackbody temp. for supplementing solar spectrum.
lower_cutoff: lower wavelength considered in model
upper_cutoff: highest wavelength (in nm) to consider in model
"""
#%%
#This cell of code can be run (if in Spyder) to verify that column indexes are set correctly before running model further down
cwd = os.getcwd()
species_string = 'Ni_I' #required for outputting plot and datafile at end.
ritz_col = 1 #required user inpu
lower_col = 5
upper_col = 6
aval_col = 3
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
temp = 5777 #Kelvin
lower_cutoff = 0
upper_cutoff = 2000
Rsun = 6.957e8 # meters
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
radfield = np.loadtxt('kurucz_solar_spectra.txt',delimiter ='\t',dtype=float) 
#For Kurucz solar spectrum provided, Col0: Wavelength [nm], Col1 Flux [ergs/cm**2/s/ster/nm] at solar surface.
supp = True #set to True if using blackbody to estimate solar intensity at lines outside the range of the 'radfield' spectrum
"""
The Kurucz line lists have flux in units of ergs/cm**2/s/ster/nm
We want W/m**2/s/sr/meter. Note: 1 erg/cm**2 = 0.001 J/m**2
We will avoid going to frequency form as the conversion is not so straightforward, and easier to work in wavelength.
Therefore, we mutiply by 0.001 * (1/1e-9) = 1e-3 * 1e9 = 1e6 to get to W/m^3/sr/
#EDIT: Corrected sign of exponent to positive on 01-05-2021
We will also take into account the distance to our object and the solid angle:
"""

radfield1 = radfield.copy()

radfield1[:,1] = 4 * np.pi * 1e6 * (Rsun/solar_dist)**2 * radfield1[:,1]
#The end result is W/m^3 (W/m**2 per wavelength in meters)
#%%
synthetic_spectra = fluorescence_spectra(radfield1,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp,lower_cutoff,upper_cutoff,supp)
synthetic_spectra2 = fluorescence_spectra(radfield2,obj_vel2,solar_dist2,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp,lower_cutoff,upper_cutoff,supp)

#%%
pops_normalized = synthetic_spectra[0]
model_lines = synthetic_spectra[1]
#
plt.figure(figsize=(16,12))
plt.rcParams.update({'font.size': 18})

plt.ylabel('Intensity [arb. units]')
plt.xlabel('Wavelength [nm]')
plt.stem(model_lines[:,0],model_lines[:,1], label = '{:} synthetic spectra'.format(species_string))
plt.legend()
plt.grid(True)
plt.xlim(model_lines[0,0] - 5, 1000)
plt.savefig('{:}_spectra.pdf'.format(species_string))
np.savetxt('{:}_data.txt'.format(species_string),model_lines,delimiter ='\t')

