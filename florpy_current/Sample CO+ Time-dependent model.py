#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SJB
Aug 9, 2023

Sample File showing some of the time-dependent functionality of FlorPy
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
# Note: It is advised NOT to mix versions.
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

    To run FlorPy in the time-dependent mode, many of the functions are identical. Only one function
    is different when running time-dependently, which is the matrix solver + gfactor calculation. The sytnax below
    will calculate the level populations and g-factors time-dependently.
    For ease of operation, I recommend running the equilibrium version first so that defining
    some of the necessary arrays is easier.

"""
#First, we build FlorPy and run the equilibrium version:
flor = build_model()
flor = add_element_to_model(flor, element_string)
flor = load_nist_data(flor,element_string,lines_file, levs_file)
flor = define_rad_fields(flor,element_string,rad_source='user',radfiles=[radiation_field], vac_or_air = ['vac'])
flor = load_orbit_data(flor,element_string,orbit_id,orbital_conds,t_comet = temp_comet, m_species = emitter_mass)
flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = 'user', profiles=prof)
flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)
#Second, we grab the population array and use it to define the initial conditions for the time-dependent run:
equil_pops = flor[element_string][orbit_id]['outputs']['pops']
gfacs = flor[element_string][orbit_id]['outputs']['gfactors']
init_popdef = np.zeros((len(equil_pops[:]),1))
init_popdef[0] = 1
#We make a unique identifier for this time-dependent run and a grid of timesteps:
td_id = 'sample_timedependent_run'
time_grid = np.geomspace(1e-2,1e3,500)
#And run the time-dependent model:
flor = calculate_pops_and_gfactors_time_dep(flor,element_string,orbit_id,td_id,time_grid,init_conds=init_popdef,time_integrated=False)
#%%
"""
Just as in the equilibrium case (see the sample spectrum file), we can auto-compute the band luminosities and grab them for saving if desired:
"""
flor = separate_vib_bands(flor,element_string,orbit_id)
flor = auto_gen_td_bandlums(flor,element_string,orbit_id,td_id)
tdbands = flor[element_string][orbit_id][td_id]['band_gfacs']
#
xx_bands = tdbands['X_X']
ct_bands = tdbands['Pi_X']
fn_bands = tdbands['B_Pi']

#%%
#   A sample plot of a level population versus time. We'll do the ground state. For reference, the time-dependent populations array has levels in energy energy (1 level per row), and
#   the number of columns is equal to the number of timesteps defined in time_grid above.
td_pops = flor[element_string][orbit_id][td_id]['time_dep_pops']
#   The td_pops array contains the time-dependent level populations. Each population is stored in a row, with time increasing as one moves to further columns
#   Ex: td_pops[0,:] -> Population of level 0 (ground) as a function of time
#       td_pops[1,:] -> population of level 1 (the first level above ground) as a function of time, and so on.
plt.clf()
plt.figure(figsize=(10,10))
plt.xscale('log')
plt.plot(time_grid[:], td_pops[0,:] / td_pops[0,-1])
plt.xlabel('Time (s)', fontsize=15)
plt.ylabel('Level population w.r.t. last time step', fontsize=15)
plt.grid()
plt.savefig('Sample_td_population.pdf',dpi=200)
#%%
# We can also show a time-dependent g-factor. We'll take a random transition somewhere near 400 nm:
td_gfacs = flor[element_string][orbit_id][td_id]['time_dep_gfactors']
#   Just like the populations, the td_gfacs array is indexed such that each row is a transition, and the columns designate the timesteps, with the last column as the last timestep.
index = 1672 #index of a random transition
plt.clf()
plt.figure(figsize=(10,10))
plt.xscale('log')
plt.plot(time_grid[:], td_gfacs[index,:])
plt.xlabel('Time (s)', fontsize=15)
plt.ylabel('Time-dependent g-factor (J/s/mol)', fontsize=15)
plt.grid()
plt.savefig('Sample_td_gfactor.pdf',dpi=200)
#%%
#
#   One could also generate a spectrum from the g-factors at a few different times:

#We grab the wavelength information:
lower_x = 200
upper_x = 700
lower_ind =  find_nearest_index(gfacs[:,0], lower_x)
upper_ind = find_nearest_index(gfacs[:,0], upper_x)    
wavelengths = gfacs[:,0]

#The plotting:
#We search the time grid for the time closest to 0.01 seconds:
time_index = find_nearest_index(time_grid, 0.01) 
plt.clf()
plt.subplots(nrows=3,ncols=1,figsize=(10,10))
plt.subplot(311)
plt.stem(wavelengths[lower_ind:upper_ind], td_gfacs[lower_ind:upper_ind,time_index], label = 't = {:} s'.format(time_grid[time_index]), markerfmt = ' ')
plt.legend()
plt.grid()
#
#We search the time grid for the time closest to 20 seconds:
time_index = find_nearest_index(time_grid, 20)
plt.subplot(312)
plt.stem(wavelengths[lower_ind:upper_ind], td_gfacs[lower_ind:upper_ind,time_index], label = 't = {:} s'.format(time_grid[time_index]), markerfmt = ' ')
plt.legend()
plt.grid()
plt.ylabel('Time-dependent g-factor (J/s/mol)', fontsize=15)
#
#We search the time grid for the time closest to 600 seconds:
time_index = find_nearest_index(time_grid, 600)
plt.subplot(313)
plt.stem(wavelengths[lower_ind:upper_ind], td_gfacs[lower_ind:upper_ind,time_index], label = 't = {:} s'.format(time_grid[time_index]), markerfmt = ' ')
plt.legend()
plt.grid()
plt.xlabel('Wavelength (air nm)', fontsize=15)
plt.subplots_adjust(hspace=0.5)

plt.savefig('Sample_plot_timedep_spectrum.pdf',dpi=200)
