#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SJB
Aug 9, 2023

Sample File showing how to run `error iterations' using Ni I data from the NIST ASD
"""
from fluor_dict_v12122023 import *
from molecular_utils_florpy import *

#Output file for the computed g-factors:
gfac_file = 'Ni0'
#Output Format:
#   Wavelength (air nm (rest))   |  Wavelength (vacuum nm (rest))    |   g-factor (J/s/mol)
element_string = 'nickel'
#Lines and Level input files:
levs_file = 'ni0_levs_txt.txt'
lines_file = 'ni0_lines_txt.txt'
# The ECO(TR) and ECO+ models can be run by simply swapping the filenames in variables levs_file and lines_file to their respective versions.
# Note: It is advised NOT to mix versions.
#Define orbital parameters:
solar_dist = 1 #0.66 #Astronomical Units
obj_vel = 0 #-1.7 #km/s; internally converted to m/s when needed.
geo_vel = 0 #km/s
geo_dist = 1 #geocentric distance; 1AU
orbital_conds = [solar_dist,obj_vel, geo_dist,geo_vel]
#Unique identifier string for this model with these orbital conditions:
orbit_id = '1au_rest'#indicates orbital conditions.
#Radiation field file; 2column format: Wavelength (vac nm) | flux (W / m^2 / nm)
radiation_field = '0-14000_vac.txt'
#Parameter for controlling line profile for all calcs:
prof = 'delta'
#If one wants to use a doppler line profile, the following two parameters must be passed to load_orbit_data:
temp_comet = 279 * (1/np.sqrt(solar_dist)) #Kelvin; comet temperature for doppler profiles
emitter_mass = 4.65e-26 #Mass of CO in kg
"""
FlorPy calculations are broken into several sequential steps. They are as follows:
    build_model() : Define a dict structure to hold the model
    add_element_to_model() : Define a dictionary entry for your unique run(s)
    load_nist_data() : Load the lines/levels files into the dictionary
    load_orbit_data(): This loads in details about the heliocentric distance / velocity (see above) into the dict structure
    define_rad_fields() : Define & load the radiation field to be used into the dict structure
    generate_fluxes_general() : A generalized function that computes the pumping rates for all included lines:
    calculate_pops_and_gfactors() : The workhorse of FlorPy. This function will generate the relevant matrices, invert them,
    and compute the level populations &  g-factors
    
    Sample ways to grab this data from the dict structure are provided.
    Sample documentation is in preparation.
"""

#First, we build FlorPy and run the equilibrium version:
flor = build_model()
flor = add_element_to_model(flor, element_string)
flor = load_nist_data(flor,element_string,lines_file, levs_file)
flor = define_rad_fields(flor,element_string,rad_source='user',radfiles=[radiation_field], vac_or_air = ['vac'])
flor = load_orbit_data(flor,element_string,orbit_id,orbital_conds,t_comet = temp_comet, m_species = emitter_mass)
flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = 'user', profiles=prof)
flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)

#%%
#Define a number of iterations. We pass many of the same 
number_of_err_iterations = 100
flor = error_iterations(flor,element_string,orbit_id, bbsupp=False, bbtemp = 5777, lower_cutoff = 1e30, upper_cutoff=1e30, aval_min = 0, bbflux_min = 0, permit_cascade=True, rad_choice = 'user', profiles = 'doppler', err_id = 'default_setting', num_iterations=number_of_err_iterations,process=True)
#%%
error_dict = flor[element_string][orbit_id]['error_calcs']['default_setting']
error_array = error_dict['processed_iteration_output']
#The contents of 'error_array' are:
#   Col0    Col1     Col2     Col3     Col4    Col5    Col6
#   Wave (air nm)   Wave (vac nm)   G-factor (j/s/mol)      Minimum g-factor    Max g-factor    Stdev
#%%
#   We can prepare a stick plot with error bars:
vscaler = 1e21 #Used to remove "e-21" from y-axis to clean up the plot:    
plt.clf()
plt.figure(figsize=(10,10))
plt.stem(error_array[:,0], error_array[:,2] * vscaler, markerfmt = 'bo',label = 'Synthetic Ni I spectra')
plt.errorbar(error_array[:,0],error_array[:,2] * vscaler, yerr=error_array[:,-1] * vscaler ,fmt='o', color = 'red',markersize=2, capsize=2)
plt.legend()
plt.title('Ni I Fluorescence Spectrum with approximate errorbars based on {:} iterations'.format(number_of_err_iterations), fontsize=15)
plt.xlim(290,400)
plt.grid(True)
plt.ylabel(r'G-factor ($10^{-21}$ J/s/mol)',fontsize=15)
plt.xlabel('Wavelength (air nm)', fontsize=15)
plt.savefig('sample_ni0_errorbar_plot.pdf', dpi=200)
