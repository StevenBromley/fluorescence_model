#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SJB
Aug 9, 2023

Sample File showing some functionality of FlorPy
"""
from fluor_dict_v12122023 import *
from molecular_utils_florpy import *

#Output file for the computed g-factors:
gfac_file = 'MAH_testrun'
#Output Format:
#   Wavelength (air nm (rest))   |  Wavelength (vacuum nm (rest))    |   g-factor (J/s/mol)
element_string = 'co+'
#Lines and Level input files:
levs_file = 'Base(MAH)_energy_levels.txt'
lines_file = 'Base(MAH)_transitions.txt'
# The ECO(TR) and ECO+ models can be run by simply swapping the filenames in variables levs_file and lines_file to their respective versions.
# Note: It is advised NOT to mix of energy level and transitions files. 
#Define orbital parameters:
solar_dist = 1 #0.66 #Astronomical Units
obj_vel = 0 #-1.7 #km/s; internally converted to m/s when needed.
geo_vel = 0 #km/s
geo_dist = 1 #geocentric distance; 1AU
orbital_conds = [solar_dist,obj_vel, geo_dist,geo_vel]
#Unique identifier string for this model with these orbital conditions:
orbit_id = 'co+'#indicates orbital conditions.
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
flor = build_model()
flor = add_element_to_model(flor, element_string)
flor = load_nist_data(flor,element_string,lines_file, levs_file)
flor = define_rad_fields(flor,element_string,rad_source='user',radfiles=[radiation_field], vac_or_air = ['vac'])
flor = load_orbit_data(flor,element_string,orbit_id,orbital_conds,t_comet = temp_comet, m_species = emitter_mass)
flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = 'user', profiles=prof)
flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)

#We grab and plot some outputs from the dict structure:
level_populations = flor[element_string][orbit_id]['outputs']['pops']
gfacs = flor[element_string][orbit_id]['outputs']['gfactors']
#Save the g-factors to a file for other processing if desired:
np.savetxt('{:}_gfacs.txt'.format(gfac_file),gfacs,delimiter='\t')

#%%
#We generate a sample stick plot of the spectrum from 200 - 700 nm:
lower_x = 200
upper_x = 700
lower_ind =  find_nearest_index(gfacs[:,0], lower_x)
upper_ind = find_nearest_index(gfacs[:,0], upper_x)
relevant_gfacs = gfacs[lower_ind:upper_ind,:]
plt.clf()
plt.figure(figsize=(10,10))
plt.stem(relevant_gfacs[:,0], relevant_gfacs[:,2], markerfmt=' ')
plt.xlim(200,700)
plt.grid(True)
plt.ylabel('G-factor (J/s/mol)',fontsize=15)
plt.xlabel('Wavelength (air nm)', fontsize=15)
plt.savefig('Sample stick spectrum.pdf',dpi=200)

#%%
#   The provided library "molecular_utils" also contains some sample functions to show how to parse the florpy output if molecules
#   are being calculated.
#   We first index the transitions into the various bands by running separate_vib_bands().
#   A second function computes the band luminosities of all of the bands in the model.
#   Sample syntax for grabbing this information from the dictionary is provided below, and the band luminosities are saved to a file:
#%%
flor = separate_vib_bands(flor,element_string,orbit_id)
#%%
flor = auto_gen_band_lum(flor,element_string,orbit_id)
#%%
x_x_bandlums = flor[element_string][orbit_id]['outputs']['band_gfacs']['X_X']
pi_x_bandlums = flor[element_string][orbit_id]['outputs']['band_gfacs']['Pi_X']
b_x_bandlums = flor[element_string][orbit_id]['outputs']['band_gfacs']['B_X']
b_pi_bandlums = flor[element_string][orbit_id]['outputs']['band_gfacs']['B_Pi']
np.savetxt('{:}_X_X_allbandlum.txt'.format(gfac_file), x_x_bandlums,delimiter='\t',fmt='%s')
np.savetxt('{:}_Pi_X_allbandlum.txt'.format(gfac_file),pi_x_bandlums,delimiter='\t',fmt='%s')
np.savetxt('{:}_B_X_allbandlum.txt'.format(gfac_file),b_x_bandlums,delimiter='\t',fmt='%s')
np.savetxt('{:}_B_Pi_allbandlum.txt'.format(gfac_file),b_pi_bandlums,delimiter='\t',fmt='%s')






