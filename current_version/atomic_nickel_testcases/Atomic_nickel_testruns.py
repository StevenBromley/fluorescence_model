# -*- coding: utf-8 -*-
"""
SJB

Example file for running FlorPy with atomic Nickel

Updated June 2024
"""
import os
cwd = os.getcwd()
os.chdir('..')
from fluor_dict_v05152024 import *
element_string = 'ni0'
#Two input files: lines file (tab-delimited) from NIST, and tab-delimited levels from NIST.
os.chdir('../input_files/atomic_nickel/')
lines_file = 'ni0_lines_txt.txt'
levs_file = 'ni0_levs_txt.txt'
#Define orbital parameters:
solar_dist = 1.0                                        #Astronomical Units
obj_vel = 0                                             #km/s; internally converted to m/s when needed.
geo_vel = 0                                             #km/s;
geo_dist = 1                                            #geocentric distance; 1AU; not used at present.
temp_comet = 279                                        #Kelvin; comet temperature for doppler profiles.
mass = 9.746267511e-26                                  #Atomic mass of Ni in kg. Used for doppler profile. 
orbital_conds = [solar_dist,obj_vel, geo_dist,geo_vel]
orbit_id = 'test'                                       #unique identifier for the calculations in the dictionary structure
#####################################################
#                                                   #
#                 RUNNING THE MODEL                 #
#                                                   #
#####################################################
#In sequence, we run:
flor = build_model()
flor = add_element_to_model(flor, element_string)
flor = load_nist_data(flor,element_string,lines_file, levs_file)
#We use the newer 0-14 micron spectrum:
radiation_field = '0-14000_vac.txt'
os.chdir('../../solar_spectra/')
flor = define_rad_fields(flor,element_string,rad_source='user',radfiles=[radiation_field], vac_or_air = ['vac'])
flor = load_orbit_data(flor,element_string,orbit_id,orbital_conds,t_comet = temp_comet, m_species = mass)
flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = 'user', profiles='delta', bbsupp=False, lower_cutoff = 1e30, upper_cutoff=1e30, aval_min=0,permit_cascade=True)
#The 'bbsupp' true indicates that a blackbody should be used to approximate wavelengths outside the radiation field, i.e. above 2730 for the 'default' file
# 'upper_cutoff' sets the cutoff value for the higher allowed energy level. Default values are shown such that *all* atomic data is included for this example.
#For alternative radiation field options, see documentation.
#Two line profiles are available: "delta" (dirac delta), and "doppler" broadened.
flor = calculate_pops_and_gfactors(flor,element_string,orbit_id,renorm=False)
lines = flor[element_string][orbit_id]['outputs']['gfactors']

#The 'lines' array contains 3 columns of data: 
#Col0 is wavelength in air nm, Col1 is wavelength in vacuum nm, col3 is gfactor in J/s/particle at the emitting source
os.chdir(cwd)
plt.clf()
plt.figure(figsize = (8,8))
plt.stem(lines[:,0], lines[:,2])
plt.ylabel(r'g-factor (J s$^{-1}$ particle$^{-1}$)')
plt.xlabel('Air Wavelength (nm)')
plt.grid(True)
plt.xlim(0,1000)
plt.savefig('ni0_gfactors.pdf')

#%%
#####################################################
#                                                   #
#      RUNNING THE MONTE-CARLO ERROR ESTIMATION     #
#                                                   #
#####################################################
#
#  One can use a Monte-Carlo method (see Bromley et al PSJ 2021) to estimate the 
# 'uncertainty' in the g-factors driven by uncertainties in the transition rates.   
# For the sake of allowing this to run in a reasonable time, let's only do 100 iterations as a demonstration

#Define a number of iterations:
number_of_err_iterations = 100
#   Note: error_iterations can take the same optional arguments as generate_fluxes_general, so we leave off the default values.
#   error_iterations is effectively a combination of generate_fluxes_general() and calculate_pops_and_gfactors(), with some additional overhead for storing the iterations.
flor = error_iterations(flor,element_string,orbit_id, rad_choice = 'user', profiles = 'delta', err_id = 'default_setting', num_iterations=number_of_err_iterations,process=True)
error_dict = flor[element_string][orbit_id]['error_calcs']['default_setting']
error_array = error_dict['processed_iteration_output']
#The contents of 'error_array' are:
#   Col0    Col1     Col2     Col3     Col4    Col5    Col6
#   Wave (air nm)   Wave (vac nm)   G-factor (j/s/mol)      Minimum g-factor    Max g-factor    Stdev
#%
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
