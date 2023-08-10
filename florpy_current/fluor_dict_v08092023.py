# -*- coding: utf-8 -*-
"""
Version 08.09.2023 of "FlorPy"

See Bromley et al PSJ 2021 and Bromley et al 2023 for details

This code is intended to be used for computing fluorescence efficiencies of atoms or molecules

See documentation for details

Distributed under a GNU GPL license

"""
import numpy as np #The OG
import matplotlib.pyplot as plt
import time
import math
import sys
import scipy
import os
from random import random #Used for generating new A values according to NIST uncertainty rating
from scipy.integrate import quad
from scipy import interpolate
from scipy.integrate import solve_ivp
from scipy.optimize import nnls

# Some useful constants (not global; just a reference for you.) All required units are defined
#inside the functions where they are used
Rsun = 6.957e8 # meters
h = 6.62607004e-34
hbar = 1.054571817e-34
k = 1.38064852e-23
c = 299792458
cm_to_J = 1.986491855e-23

def build_model():
    #06-11-2021 SJB Initiates Dict structure for flor model
    flor = {}
    return(flor)

def add_element_to_model(flor,element_string):
    """
    SJB
    06-11-2021
    This function takes in the flmod dict and adds to it a dict for a new element
    i.e. flor['fe0'] for neutral 0. The key of the dict is added in as a string;
    Note there is no requirement on the naming convention; if I add in silicon but name it iron,
    nothing breaks. The element_string is best if YOU understand what each entry into flor contains.
    
    Only one entry with a given key can exist. 
    
    """
    if (element_string not in flor):
        flor[element_string] = {}
    else:
        print('Element key "{:}" already exists in dictionary! Ignoring addition request'.format(element_string))
    return(flor)

def load_nist_data(flor,element_string,lines_file, levels_file):
    """
    SJB
    06-11-2021
    This function reads in NIST lines and levels files downloaded by the user from NIST ASD and adds the data to the model.
    """
    load_start = time.time()
    print('Loading Atomic Data...')
    raw_lines_un = np.genfromtxt(lines_file, delimiter='\t', dtype = str)
    raw_levs_un = np.genfromtxt(levels_file,delimiter ='\t',skip_header=1,usecols=(0,1,2,3),dtype=str)
    raw_lines = np.copy(raw_lines_un)
    raw_levs = np.copy(raw_levs_un)
    #We strip all of the extra "" characters:
    for i in range(0,len(raw_lines[:,0])):
        for j in range(0,len(raw_lines[0,:])):
            raw_lines[i,j] = raw_lines[i,j].replace('"','')
    #Grab column headers needed:        
    for i in range(0,len(raw_lines[0,:])):
        if ('ritz' in raw_lines[0,i]):
            ritz_col = i
        if ('Aki' in raw_lines[0,i]):
            aval_col = i
        if ('Ei' in raw_lines[0,i]):
            lower_col = i
        if ('Ek' in raw_lines[0,i]):
            upper_col = i 
        if ('Acc' in raw_lines[0,i]):
            uncert_col = i 
        #J value checks added on Sep 23, 2022 to handle degenerate in energy levels:
        if ('J_i' in raw_lines[0,i]):
            lowerj_col = i
        if ('J_k' in raw_lines[0,i]):
            upperj_col = i
    #We need to know the number of header rows:
    header_rows = np.zeros((0,1),dtype=int)
    for i in range(0,len(raw_lines[:,0])):
        if ('Ei' in raw_lines[i,lower_col]):
            temp = np.zeros((1,1))
            temp[0,0] = i
            header_rows = np.append(header_rows,temp,axis=0)
    #We swap the order of the array,i.e. descending order:
    header_rows = header_rows[::-1]
    if (len(header_rows) > 0):
        for i in range(len(header_rows[:])):
            raw_lines = np.delete(raw_lines,int(header_rows[i]),0)
    #Save the "raw_lines" data:
    flor[element_string]['lines_data_str'] = raw_lines
    #We convert all necessary values to float. This will also remove troublesome characters:
    lines = np.zeros((len(raw_lines[:,0]),len(raw_lines[0,:])),dtype=float)
    #Save data in usable form and convert half-integer J's to floats where necessary in the lines and levels arrays:
    for i in range(0,len(lines[:,0])):
        for j in range(0,len(lines[0,:])):
            try:
                lines[i,j] = convert_to_float(raw_lines[i,j])
            except ValueError:
                pass
    ####
    ####    Following added on 08-10-2021:
    ####
    #For the float data, we convert the accuracy rating over to a float:
    uncert_grabbed = grab_uncertainty(raw_lines,uncert_col)
    lines[:,uncert_col] = np.reshape(uncert_grabbed,len(uncert_grabbed))
    #Note: The reshape is necessary to re-cast in appropriate format
    flor[element_string]['lines_data_float'] = lines
    ####
    ####    End of 08-10-2021 addition
    ####
    #Save the column indices we'll need:
    flor[element_string]['column_indices'] = {'ritz_col' : ritz_col,
                                              'aval_col' : aval_col,
                                              'lower_col' : lower_col,
                                              'upper_col' : upper_col,
                                              'uncert_col' : uncert_col,
                                              'lowerj_col' : lowerj_col,
                                              'upperj_col' : upperj_col}
    #%
    #Now we process the levels. First, remove extra "" characters:
    for i in range(0,len(raw_levs[:,0])):
        for j in range(0,len(raw_levs[0,:])):
            raw_levs[i,j] = raw_levs[i,j].replace('"','')
    levs = np.zeros((len(raw_levs[:,0]),len(raw_levs[0,:])),dtype=float)
    bad_levels_out = np.empty((0)) #Returns nothing for a good run; for a bad run, will contain the probable 'bad' levels found through model iterations.
    for i in range(0,len(levs[:,0])):
        for j in range(0,len(levs[0,:])):
            try:
                levs[i,j] = convert_to_float(raw_levs[i,j])
            except ValueError:
                pass
    
    flor[element_string]['levs_str'] = raw_levs
    #Added sort to following line on 06-22-2021
    flor[element_string]['levs_float'] = levs[np.argsort(levs[:,3])]
    flor[element_string]['bad_levels'] = bad_levels_out #Saved for later
    load_end = time.time()
    print('Lines and levels loaded for model "{:}" in {:} seconds'.format(element_string, round(load_end - load_start,2)))
    return(flor)

def define_rad_fields(flor,element_string,rad_source,bbtemp = 5777, radfiles=[], vac_or_air = []):
    """
    SJB
    06-11-2021
    This function defines the radiation field(s) used. 
    
    IF rad_source = 'default', the provided file "0-2730_vac.txt' is used. This file is constructed from:
        UV:
          0 - 168 nm: the FUV model '13x1' of Fontentla et al 2014 is used, 
          168 - 200.06nm (vac): we assumed the solar spectrum collected during solar minimum from the SOLSTICE instrument with 0.1 nm resolution
          200.06 - 202nm: the Hall & Anderson (1991) data set, matched at the edges, and
          202 - 2730 nm: the combined & calibrated spectrum from Coddington et at 2021
         All wavelengths for the default field are in vacuum.
            
    A blackbody choice is provided with a default temp of  5777K (~solar). User is free to set any temp they wish.
    This option can be used as the entire source of excitation, or only for "supplementing" the radiation field (discussed later)
    
    IF rad_source = 'user', the user must provide one of two forms:
        Form 1:
            User provides one big 2-col array of : col0 wavelength in VACUUM (nm), col1 flux (W / m^2 / nm) at 1 AU
        Form 2:
            User provides multiple files. For example, if the user provided 2 files called spec1 and spec2:
                spec1: 0 - 300 nm (vac), and 
                spec2: 500 - 2000 nm (air), 
            The inputs 'radfiles' and vac_or_air would look like:
                radfiles = ['spec1.txt', 'spec2.txt']
                vac_or_air = ['vac', 'air']
            For points outside these ranges, the user can either set rates to 0 (not recommended; leads to potentially bad solutions) by setting bbsupp = False later,
            or use a blackbody to estimate fluxes at wavelengths outside those bounds; for the above case, this would be for lines between 300 - vac_to_air(500),
            and above 2000. Using a blackbody to estimate fluxes would be setting the variable bbsupp = True later.
    If radfield = 'blackbody', a blackbody of temperature "bbtemp" is used for all pumping rate calculations.

    """    
    define_rad_start = time.time()
    if ('rad_choices' not in flor[element_string]):
        flor[element_string]['rad_choices'] = {}
    if (rad_source == 'default'):
        print('Default Radiation fields selected. Loading data...')
        flor[element_string]['rad_choices']['default'] = {}
        #Load the 'default' files
        flor[element_string]['rad_choices']['default']['all'] = {}
        flor[element_string]['rad_choices']['default']['all']['medium'] = 'vac'
        flor[element_string]['rad_choices']['default']['all']['lims'] = np.zeros((1,2))
        flor[element_string]['rad_choices']['default']['all']['fluxes'] = np.genfromtxt('0-2730_vac.txt', delimiter='\t')
        flor[element_string]['rad_choices']['default']['all']['lims'][0,0] = flor[element_string]['rad_choices']['default']['all']['fluxes'][0,0]
        flor[element_string]['rad_choices']['default']['all']['lims'][0,1] = flor[element_string]['rad_choices']['default']['all']['fluxes'][-1,0]
        flor[element_string]['rad_choices']['blackbody'] = {}
        flor[element_string]['rad_choices']['blackbody']['temperature_K'] = bbtemp
        define_rad_end = time.time()
        print ('Default solar spectrum (see docs) loaded in {:} seconds'.format(round(define_rad_end - define_rad_start,2)))
    if (rad_source == 'user'):
        print('Custom user-defined radiation inputs selected')
        flor[element_string]['rad_choices']['user'] = {}
        for i in range(0,len(radfiles[:])):            
            radfile = radfiles[i]
            #Define new dict:
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)] = {}
            #Load data into 'fluxes' key:
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['fluxes'] = np.genfromtxt(radfile,delimiter='\t')        
            #define an array for the lower/upper bounds:
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['lims'] = np.zeros((1,2))
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['lims'][0,0] = flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['fluxes'][0,0]
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['lims'][0,1] = flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['fluxes'][-1,0]
            flor[element_string]['rad_choices']['user']['rad{:}'.format(i)]['medium'] = vac_or_air[i] 
        define_rad_end = time.time()
        print ('Custom solar spectrum loaded in {:} seconds'.format(round(define_rad_end - define_rad_start,2)))  
    return(flor)


def load_orbit_data(flor,element_string,orbit_id,orbital_conds, t_comet = 280, m_species = 9.74620109e-26):
    """
    SJB
    Documentation Updated Aug9
    
    This function takes in the florpy dict structure, the element_string and orbit_id identifiers, and the orbital conditions
    the orbital conditions must be saved in a 1-D array such as:
        orbital_conds = [heliocentric_distance, heliocentric_velocity, geocentric distance, geocentric velocity, temperature, mass]
        where distances are in AU, velocities are in km/s (defined as positive = away from the sun/earth)
        
    Optional Parameters "t_comet" and "m_species" are in Kelvin and kg, respectively, and are only required
    if the user wishes to use a doppler line profile. Defaults are 280K and the mass of atomic nickel

    Each set of distance/velocity is assigned a unique id number, and all calculated data will be saved under that ID number (the orbit_id).
    """
    if orbit_id not in flor[element_string]:
        flor[element_string][orbit_id] = {}
    flor[element_string][orbit_id]['orbit_params'] = {}
    flor[element_string][orbit_id]['orbit_params']['helio_dist_au'] = orbital_conds[0]
    flor[element_string][orbit_id]['orbit_params']['helio_vel_kms'] = orbital_conds[1]
    flor[element_string][orbit_id]['orbit_params']['geo_dist_au'] = orbital_conds[2]
    flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] = orbital_conds[3]
    flor[element_string][orbit_id]['orbit_params']['comet_temp'] = t_comet #default 100K
    flor[element_string][orbit_id]['orbit_params']['emitter_mass'] = m_species #default to nickel
    print('Orbital Conditions Loaded:')
    print('Distance from Radiation Source: {:} AU'.format(orbital_conds[0]))
    print('Velocity w.r.t. Radiation Source: {:} km/s'.format(orbital_conds[1]))
    return(flor)


def generate_fluxes_general(flor,element_string,orbit_id,rad_choice='default', bbsupp=False, bbtemp = 5777, profiles='delta', lower_cutoff = 1e30, upper_cutoff=1e30, aval_min = 0, bbflux_min=0, permit_cascade=True):
    """
    SJB
    06-12-2021
    Flux calculation function for the dictionary version of florpy.

    Required Inputs:
    flor: input dict after loading atomic data and radiation fields
    element_string: specifies which element is being calculated.
    orbital_conds: list of orbital parameters; distance in AU, velocity in km/s
    orbit_id: user-provided string to indicate that the dict should have a unique entry based on the orbital conditions
    This can be useful if one wanted to run models with multiple sets of distances/velocities and compare outputs easily
    
    Optional Inputs:
    rad_choice: choice of radiation field. The two options are 'default' and 'user'
    bbsupp: Boolean (True or False). If bbsupp=True, then a blackbody of temperature "bbtemp" (default to 5777K)
    will be used to compute the pumping rate for the lines outside the bounds of the provided radiation field(s).
    Note: This appears to cause some problems if used for lines with with very large wavelengths.
    profiles: Optionals are 'delta' (default) or 'doppler'. These are case-sensitive options.   
     
    Cutoffs:
        lower_cutoff: sets the maximum value of the LOWER level energy in cm-1 for included lines / levels
        upper_cutoff: sets the maximum value of the UPPER level energy in cm-1 for included lines / levels.
        aval_min: sets the minimum transition rate to include in the model
    bbflux_min: set the minimum rate for inclusion in the model when using bbsupp=True
    
    permit_cascade: Boolean (True or False). Default is True. When True, this parameter will include lines in the model whose wavelengths are outside the radiation field,
    but their pumping rates are set to 0. This allows models to be included if they act as weak perturbers or connect levels, but otherwise are not important in the rate equations.
    I recommend keeping this set to True unless you know 100% what you are doing.
    
    """
    print('Starting Calculation of Integrated Fluxes...')
    start = time.time()
    Rsun = 6.957e8 # meters
    #Grab necessary data from dict structure:
    #turned off prep_indexing in favor of line_check on Jan10 2023
    flor = prep_indexing(flor,element_string,orbit_id,rad_choice, bbsupp, lower_cutoff, upper_cutoff, aval_min, permit_cascade)
    raw_lines = flor[element_string][orbit_id]['relevant_lines']
    #Column indexes
    lower_col = flor[element_string]['column_indices']['lower_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    integrated_fluxes = np.zeros((len(raw_lines[:,0]),1),dtype=float)
    #Grab orbital data and convert as necessary to SI:
    solar_dist = 1.496e+11 * flor[element_string][orbit_id]['orbit_params']['helio_dist_au'] #meters
    obj_vel = 1e3 * flor[element_string][orbit_id]['orbit_params']['helio_vel_kms'] #m/s
    dist_scale_1au = (1.496e+11/solar_dist)**2
    #We want to loop over the choice of radiation fields and generate fluxes for each line in that range.
    upper_energies = (raw_lines[:,upper_col])
    lower_energies = (raw_lines[:,lower_col])
    vac_waves = 1e7 / (upper_energies - lower_energies)     
    air_waves = vac_to_air(vac_waves)
    all_waves = np.concatenate((vac_waves.reshape(vac_waves.shape[0],-1),air_waves.reshape(vac_waves.shape[0],-1)),axis=1)
    #
    if (rad_choice != 'blackbody'):
        if (profiles == 'delta'):
            for key in flor[element_string]['rad_choices'][rad_choice]:
                if (flor[element_string]['rad_choices'][rad_choice][key]['medium'] == 'vac'):
                    wave_col = 0
                else:
                    wave_col = 1
                search_idx_iter = 0
                for i in range(0,len(all_waves[:,0])):
                    if (flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,0] < doppler(all_waves[i,wave_col],obj_vel) < flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,1]):
                        #
                        searched_rad_indices = fniiv(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0],doppler(all_waves[i,wave_col],obj_vel),search_idx_iter)
                        #start_idx is passed to next iteration.
                        search_idx_iter = searched_rad_indices[1]
                        #Index for the lines radiation field:
                        rad_idx = searched_rad_indices[0]
                        #FOllowing edited on Sep 26 2022 (removed 1e9 conv.)   
                        rad_at_wav = dist_scale_1au * flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][rad_idx,1] * 1e9 #converts from W/m^2 per nm to W /m^2 per meter (wavelength)
                        #rad_at_wav = dist_scale_1au * flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][rad_idx,1]
                        integrated_fluxes[i] = rad_at_wav
            #If bbsupp is true, then there may be lines for which the imported radiation fields are insufficient. We calculate fluxes for them:
            if (bbsupp == True):
                for i in range(0,len(integrated_fluxes[:])):
                    if (integrated_fluxes[i] == 0):
                        dist_scale  = (Rsun/solar_dist)**2
                        bb = blackbody_wave(doppler(vac_waves[i],obj_vel),temp=bbtemp) #Only need that specific wavelength at temp=5777K (default); assumes nm in
                        rad_at_wav = dist_scale * bb
                        integrated_fluxes[i] = rad_at_wav
        if (profiles == 'doppler'):
            comet_temp = flor[element_string][orbit_id]['orbit_params']['comet_temp']
            m_species = flor[element_string][orbit_id]['orbit_params']['emitter_mass']
            for key in flor[element_string]['rad_choices'][rad_choice]:
                if (flor[element_string]['rad_choices'][rad_choice][key]['medium'] == 'vac'):
                    wave_col = 0
                else:
                    wave_col = 1
                search_idx_iter = 0
                for i in range(0,len(all_waves[:,0])):
                    if (flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,0] < doppler(all_waves[i,wave_col],obj_vel) < flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,1]):
                        profile = doppler_dist(doppler(all_waves[i,wave_col],obj_vel),comet_temp,m_species)
                       # dist_scale_1au = (1.496e+11/solar_dist)**2
                        searched_rad_indices = fniiv(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0],doppler(all_waves[i,wave_col],obj_vel),search_idx_iter)
                        rad_idx = searched_rad_indices[0]
                        search_idx_iter = searched_rad_indices[1]
                        # We want to auto-detect what wavelengths range to pass in.
                        # Two counters for the "mini loops" which are used to find start/end indices automatically
                        # For the lower one, we intentionally offset to the red and search backwards.
                        #For the upper index, we intentionally offset to the blue and search forwards:
                        lower_rad_idx = np.copy(rad_idx) + 5
                        upper_rad_idx = np.copy(rad_idx) - 5
                        #
                        for p in range(0,len(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0])):
                            if (flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][lower_rad_idx,0] > profile[0,0]):
                                lower_rad_idx = lower_rad_idx - 1
                            else:
                                break
                        for p in range(0,len(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0])):
                            if (flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][upper_rad_idx,0] < profile[-1,0]):
                                upper_rad_idx = upper_rad_idx + 1
                            else:
                                break
                        #Buffer with extra 2. Will likely be removed, but effect is minor at best.
                        rad_start_idx = lower_rad_idx - 2
                        rad_end_idx = upper_rad_idx + 2
                        new_rad_grid = generate_rad_grid(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][rad_start_idx:rad_end_idx,:],profile[:,0])
                        #Use fast sci-py module for the integration. ~10x increase in speed over homemade older integrator.
                        rad_at_wav = dist_scale_1au * 1e9 * scipy.integrate.simps(new_rad_grid[:,1] * profile[:,1], profile[:,0])
                        integrated_fluxes[i] = rad_at_wav
            if (bbsupp == True):
                dist_scale  = (Rsun/solar_dist)**2
                for i in range(0,len(integrated_fluxes[:])):
                  if (integrated_fluxes[i] == 0):
                      profile = doppler_dist(doppler(all_waves[i,wave_col],obj_vel),comet_temp,m_species)
                      bb = blackbody_wave(profile[:,0],temp=bbtemp) #Only need the profile wavelengths at temp=5777K (default); assumes nm in
                      rad_at_wav = dist_scale * scipy.integrate.simps(bb * profile[:,1], profile[:,0])
                      integrated_fluxes[i] = rad_at_wav         
    #Handle blackbody in a similar way. Less accurate, but less cumbersome to deal with since we don't have to do any array searching:
    if (rad_choice == 'blackbody'):
        wave_col = 0
        comet_temp = flor[element_string][orbit_id]['orbit_params']['comet_temp']
        m_species = flor[element_string][orbit_id]['orbit_params']['emitter_mass']
        dist_scale  = (Rsun/solar_dist)**2
        if (profiles == 'delta'):
            for i in range(0,len(integrated_fluxes[:])):
                bb = blackbody_wave(doppler(vac_waves[i],obj_vel),temp=bbtemp) 
                rad_at_wav = dist_scale * bb
                integrated_fluxes[i] = rad_at_wav
        if (profiles == 'doppler'):
            for i in range(0,len(integrated_fluxes[:])):
                profile = doppler_dist(doppler(all_waves[i,wave_col],obj_vel),comet_temp,m_species)
                bb = blackbody_wave(profile[:,0],temp=bbtemp) #Only need that specific wavelength at temp=5777K (default); assumes nm in
                #new:
                rad_at_wav = dist_scale * scipy.integrate.simps(bb * profile[:,1], profile[:,0])
                integrated_fluxes[i] = rad_at_wav    
    #We perform a final check to perform a threshold on the very weak fluxes. Default is 0.
    if (bbsupp==True):        
        for i in range(0,len(integrated_fluxes[:])):
               if (integrated_fluxes[i] < bbflux_min):
                   integrated_fluxes[i] = 0
    end = time.time()
    print('Calculation of Integrated Fluxes finished. Elapsed Time {:} s'.format(round(end - start,2)))
    flor[element_string][orbit_id]['fluxes'] = integrated_fluxes
    return(flor)    

def prep_indexing(flor,element_string,orbit_id,rad_choice, bbsupp, lower_en_cutoff, upper_en_cutoff, aval_min,permit_cascade=True):
    """
    06-22-2021 
    SJB
    The purpose of this function is to prepare the lines/levels list for ONLY those that are within the limits of the radiation fields.
    If we include the lines without radiation field components and the blackbody suppelementation is applied, we'll have
    a singular matrix and the solution methods will struggle. 
    This function also applies the supplied cutoffs to the lines/levels included in the model.
    """
    print('Checking raw line list for relevant lines...')
    #First, we make a copy of the TOTAL lines file.
    checklines = np.copy(flor[element_string]['lines_data_float'])
    checklines_str = np.copy(flor[element_string]['lines_data_str'])
    upper_col = flor[element_string]['column_indices']['upper_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    aval_col = flor[element_string]['column_indices']['aval_col']
    #Row indexes of lines to eventually REMOVE from checklines:
    idx_to_remove = np.zeros((len(checklines[:,0]),2),dtype=int)
    #Generate all of the limits:
    lims = np.empty((0,2))
    media = np.empty((0,1), dtype = str)
    obj_vel = flor[element_string][orbit_id]['orbit_params']['helio_vel_kms'] * 1e3 #convert to m/s
    for key in flor[element_string]['rad_choices'][rad_choice]:
        lims_range = flor[element_string]['rad_choices'][rad_choice][key]['lims']
        lims = np.append(lims,lims_range,axis=0)
        medium = np.empty((1,1),dtype=str)
        medium[0,0] = flor[element_string]['rad_choices'][rad_choice][key]['medium']
        media = np.append(media,medium,axis=0)
        for i in range(0,len(checklines[:,0])):
            idx_to_remove[i,0] = i
            #We compute the vacuum wavelength from the energy levels directly, then convert to air as needed:
            wave_vac = doppler(1e7 /(checklines[i,upper_col] - checklines[i,lower_col]),obj_vel) #
            wave_air = doppler(vac_to_air(wave_vac),obj_vel)
            for j in range(0,len(lims[:,0])):
                if ('v' in media[j,0]):
                    wave_to_use = wave_vac
                else:
                    wave_to_use = wave_air
                if ((lims[j,0] < wave_to_use < lims[j,1]) and (checklines[i,upper_col] < upper_en_cutoff) and (checklines[i,lower_col] < lower_en_cutoff) and (checklines[i,aval_col] > aval_min)):  #Check wavelength limits and upper energy cutoff
                    idx_to_remove[i,1] = 1
    #The only lines w/o a '1' to indicate keeping the line are outside the radiation field. We flag those for keeping if bbsupp == True:
    #Added permit_cascade check on Jan 23, 2023
    for i in range(0,len(checklines[:,0])):
        if ( ( (bbsupp == True) or (permit_cascade == True) ) and (checklines[i,upper_col] < upper_en_cutoff) and (checklines[i,lower_col] < lower_en_cutoff) and (checklines[i,aval_col] > aval_min)):
            idx_to_remove[i,1] = 1
    #After this procedure, all lines w/ a valid radfield will have a '1' in their corresponding entry in idx_to_remove
    #We set by descending order of col0:
    idx_to_remove = idx_to_remove[idx_to_remove[:,0].argsort()[::-1]]    
    for i in range(0,len(idx_to_remove[:,0])):
        if (idx_to_remove[i,1] == 0):
            checklines = np.delete(checklines,idx_to_remove[i,0],0)
            checklines_str = np.delete(checklines_str,idx_to_remove[i,0],0)
    #The checklines array with deleted rows is now what we want.
    flor[element_string][orbit_id]['relevant_lines'] = checklines
    flor[element_string][orbit_id]['relevant_lines_str'] = checklines_str
    return(flor)

def doppler(wavelength,velocity):
    """
    SJB
    Function for calculating the doppler shift of a line center.
    Inputs:
    wavelength in nm.
    velocity is in m/s (positive is AWAY from radiation source)
    outputs doppler shifted wavelength.

    2023 edits:
    Swapped signs of velocities as listed below:
       frac = (1 - velocity/c) / (1 + velocity/c) instead of (1 + velocity/c) / (1 - velocity/c)
       to correct for shifting of cometary wavelengths instead of solar spectrum
       See discussion in Schleicher & A'Hearn (1988) Section IIb
       This is much less cumbersome than doppler shifting the entire solar spectrum, and reduces the computational demands substantially
       for the same end result.
    """
    c = 299792458
    rest_wavelength = wavelength #nm
    frac = np.sqrt((1 - velocity/c)/(1 + velocity/c))
    shifted_wavelength = rest_wavelength * frac
    return(shifted_wavelength)

def blackbody_wave(wavelength, temp=5777):
    """
    SJB
    Function for a blackbody radiation distribution.
    Note that this has funky units: input wavelength is nm, and output is in 
    W / m^3.
    Temperature defaults to solar (5777 K)
    Integrates to approximately the correct energy output of the Sun when handled properly.
    """
    wavelength = wavelength * 1e-9
    h = 6.62607004e-34 #planck constant in SI
    c = 299792458 #speed of light in m/s (SI)
    k = 1.38064852e-23 #Boltzmann constant in SI
    exp_factor = h*c / (wavelength * k * temp) #compute item inside exponential factor
    exp_term = 1/(np.exp(exp_factor) - 1) #calculate exponential. Then, calculate BB intensity (W/m^2 per wavelength(m) per ster.):
    intensity = (2 * h * c**2 / (wavelength**5) )* exp_term
    intensity = np.pi * intensity
    #The factor of pi corrects for solid angle; 
    #Above function integrates to 1358 W/m^2 when integrated and projected to 1AU by (Rsun (m) / 1 AU (m))^2    
    return(intensity)

def find_nearest_index(array,value):
    """
    SJB
    Custom function to find the index of the element nearest in value to the desired "value"
    Not very fast for large arrays, but works well enough for small arrays and tasks:
    """
    diff_temp = abs(value-array[0])
    index = 0 
    for i in range(0,len(array[:])):
        diff = abs(array[i] - value)
        if (diff < diff_temp):
            index = i
            diff_temp = diff
    return index

def vac_to_air(wavelength_vac):
    """
    SJB
    IAU Standard Conversion (Morton 2000)
    for converting vacuum wavelengths to air wavelengths.
    
    Input wavelength is in VACUUM nm. 
    """
    wavelength_vac = wavelength_vac * 10 #convert to angstroms
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength = wavelength_vac / n
    wavelength_conv = wavelength / 10 #convert back to nm
    return(wavelength_conv)
#%%
def convert_to_float(frac_str):
    """
    Updated 06-11-2021 by SJB
    This converts a string to float; it is intended to do 2 things:
     1. Remove flag characters from "questionable" levels in NIST level files, and 
     2. Convert half-integer J's, e.g. 3/2, to decimals (1.5)
    Known questionable flags (to date) are: 'a', brackets i.e. '[ ]', '*', and '?'
    These problematic characters show up in some NIST data for atomic species when there are problem(s) with certain energy levels
    in experimental analyses.
    """
    number_string = (frac_str.split('?'))[0]
    try:
        number_string = number_string.replace('[', '')
        number_string = number_string.replace(']', '')
        number_string = number_string.replace('*', '')
        number_string = number_string.replace('a', '')
        return float(number_string)
    except ValueError:        
        top, bottom = number_string.split('/')
        try:
            #If >1:
            leading, top = top.split(' ')
            whole = float(leading)
        except ValueError: #if <1 handle thusly:
            whole = 0
        frac = float(top) / float(bottom)
        if (whole < 0):
            return(whole - frac)
        else:
            return(whole + frac)  

#%%
def grab_uncertainty(lines,col):
    """
    SJB
    01-08-2021. This function is used for converting A value uncertainty ratings in NIST to numbers
    pass in raw_lines array, and 'col' is the col# of the accuracy ratings; using standard input names,
    this is the "uncert_col" variable
    """
    a_uncert = np.zeros((len(lines[:,0]),1),dtype=float) #output array A value uncertainties
    for i in range(0,len(lines[:,0])):
        letter = lines[i,col]
        if (letter == 'AAA'):
            a_uncert[i,0] = 0.3
        if (letter == 'AA'):
            a_uncert[i,0] = 1
        if (letter == 'A+'):
            a_uncert[i,0] = 2
        if (letter == 'A'):
            a_uncert[i,0] = 3
        if (letter == 'B+'):
            a_uncert[i,0] = 7
        if (letter == 'B'):
            a_uncert[i,0] = 10
        if (letter == 'C+'):
            a_uncert[i,0] = 18
        if (letter == 'C'):
            a_uncert[i,0] = 25
        if (letter == 'D+'):
            a_uncert[i,0] = 40
        if (letter == 'D'):
            a_uncert[i,0] = 50
        if (letter == 'E'):
            a_uncert[i,0] = 100
    #Return the equivalent uncertainty rating in terms of a percentage:
    return(a_uncert)
#%%
def generate_new_avals(avals,uncerts):
    """
    Added on 01-08-2021
    Pass in str array of a values, and a 1d array of uncertainties (floats)
    Then, generate a new A value within the range Aval +/- uncert.
    """
    adjusted_avals = np.empty((len(avals[:]),1),dtype='U256')
    for i in range(0,len(avals[:])):
        uncert = uncerts[i]
        aval = float(avals[i])
        low_aval = aval - aval*uncert*0.01 #change from percent to decimal
        high_aval = aval + aval*uncert*0.01 #change from percent to decimal
        adjusted_avals[i] = (str(np.random.uniform(low_aval,high_aval)).replace('[','')).replace(']','')
    return(adjusted_avals)

def generate_new_avals_float(avals,uncerts):
    """
    Added on 08-10-2021
    Generates new A values (transition rates) for array of input array values and uncertainties. 
    Uncertainties are imported in the form of a percent, and new values are chosen within original value +/- uncertainty.
    """
    adjusted_avals = np.empty((len(avals[:]),1),dtype='U256')
    for i in range(0,len(avals[:])):
        uncert = uncerts[i]
        aval = avals[i]
        low_aval = aval - aval*uncert*0.01 #change from percent to decimal
        high_aval = aval + aval*uncert*0.01 #change from percent to decimal
        adjusted_avals[i] = np.random.uniform(low_aval,high_aval)
    return(adjusted_avals)

def error_process(errdat):
    """
    Checked on 02-25-2021.
    This function processes the output of "error_calc" and calculates standard deviations of line intensities following the iterations performed in
    error_calc. This function is intended to directly pass the output of error_calc into this function.
    """
    import numpy as np
    #Pre-allocate the array size; quicker than using np.append to add new rows every time:
    err_out = np.zeros((len(errdat[:,0]),6)) 
    err_out[:,0] = errdat[:,0] #Air wave
    err_out[:,1] = errdat[:,1] #vac Wave 
    err_out[:,2] = errdat[:,2] #Intensity
    for i in range(0,len(errdat[:,0])):
        max_val = np.amax(errdat[i,3:])
        min_val = np.amin(errdat[i,3:])
        stdev = np.std(errdat[i,3:])
        #Save:
        err_out[i,3] = min_val
        err_out[i,4] = max_val
        err_out[i,5] = stdev
    return(err_out)
  
#%%
def fniiv(array,value,start_idx):   
    """
    " Nearest Index Iterative Version"
    Added on 01-27-2021
    This is essentially a "find nearest index" code for a sorted array,
    where you pass the place you stopped to the next iteration to prevent re-searching parts of the array that are already searched.
    """
    index = 0
    for i in range(start_idx,len(array[:])):
        diff = array[i] - value
        if (diff > 0):
            lower_val = array[i-1]
            if (abs(lower_val) < abs(diff)):
                index = i-1
                break
            else:
                index = i
                break
    if (index > 100):
        next_index = i
    else:
        next_index = i-25        
    return(index,next_index)
#%%
def zeros_check(arr):
    """
    SJB
    Added on 03-02-2012.
    Searches for zeros in input array. If present, break and return True. 
    Used to verify that the matrix solver returned a valid solution
    """
    check_bool = False
    for p_check in range(0,len(arr[:])):
        if (arr[p_check] == 0):
           check_bool = True
           break
    return(check_bool)
#%%
def matrix_solve(lhs,rhs):
    """
    Updated on Jan 17, 2023 to handle larger matrices. 
    
    First, we check the condition number of the rate matrix. If it is too large, regular inverse methods may fail. We default to using the pseudo-inverse if the condition is OK, as this
    will automatically using a standard inverse method before attempting pseudo-inverse in the event the regular method fails.
    
    if the condition number is large, a least squares solution to Ax=B is more appropriate. In that case, we want to use a least squares (np.linalg.lstsq specifically) to solve Ax=B
    In cases where the least squares solution runs into trouble, we default to a non-negative least-squares solver with scipy. 
    """
    solve_start_time = time.time()
    pops = np.zeros((len(rhs),1))
    #We need to do this iteratively. For most systems, the solution easily follows from np.linalg.inv.    

    condition_number = np.linalg.cond(lhs)
    print('Matrix Condition Number: {:.2e}'.format(condition_number))
    #Condition cutoff set to 1e30 instead of 1e12.
    if (condition_number < 1e30):
        a_inverse = np.linalg.pinv(lhs)
        pops = np.dot(a_inverse,rhs)    
        if (zeros_check(pops) == False):
            solve_end = time.time()
            print('Matrix solution with Numpy Inverse methods completed in {:} seconds. Checking for negative values...'.format(round(solve_end - solve_start_time,3)))
            neg_check = 0
            #We loop over the populations and check for negative values. The below check will multiply the populations by -1 if a negative value is found. 
            #If the entire solution is negative, this is not a problem; a vector times a scalar is still a vector, and could be a valid solution potentially (need to double check some texts)
            for i in range(0,len(pops[:])):
                if (pops[i] < 0):
                    pops = pops * -1
                    neg_check = neg_check + 1
            if (neg_check > 1):                
                print('Problematic negative populations present in matrix solution.')
            else:
                print('No negative populations persisting in matrix solution. Proceeding...')
        else:
            print('Zeros detected in final level populations! Check atomic data inputs.')
    else:
        print('Large condition number detected for rate matrix. Proceeding with least-squares solver.')
        #In this case, the matrix has a very large condition number and there is some doubt in the standard matrix inverse solution method. Instead, we apply a least squares solver:
        lstsq_solution = np.linalg.lstsq(lhs,rhs,rcond=None)
        pops = lstsq_solution[0]
        if (zeros_check(pops) == False):
            solve_end = time.time()
            print('Least-squares matrix solution completed in {:} seconds. Checking for negative values...'.format(round(solve_end - solve_start_time,3)))
            neg_check = 0
            #We loop over the populations and check for negative values. The below check will multiply the populations by -1 if a negative value is found. 
            #If the entire solution is negative, this is not a problem; a vector times a scalar is still a vector, and could be a valid solution potentially (need to double check some texts)
            for i in range(0,len(pops[:])):
                if (pops[i] < 0):
                    pops = pops * -1
                    neg_check = neg_check + 1
            if (neg_check > 1):                
                print('Problematic negative populations present in matrix solution.')
            else:
                print('No negative populations persisting in matrix solution. Proceeding...')
                
    if (neg_check > 1):
        print('Persistent Negative values in solution. Attempting SciPy Non-Negative Least Squares solver...')
        pops = nnls(lhs,rhs.reshape(-1),maxiter=1e10)[0]
        neg_check = 0
        for i in range(0,len(pops[:])):
            if (pops[i] < 0):
                pops = pops * -1
                neg_check = neg_check + 1
        if (neg_check > 1):                
            print('Problematic negative populations present in matrix solution. Be wary of NNLS solution!')
        else:
            print('Valid solution achieved with Non-Negative Least Squares solver.')            
    return(pops)    
#%%
def generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col,lowerj_col,upperj_col):
    """
    Added 03-02-2021
    SJB
    Used internally to generate matrix "idx" array used to encapsulate all rate data that goes into the LHS
    matrix of our matrix equation.
    
    Input variables are pulled from elsewhere as determined in the initial data reading. 
    
    """    
    c = 299792458
    #First, find the LEVELS for the LINES
    # Array 'idx' is the same length as the lines array and will hold the indexes corresponding to the index of the 
    # lower level (col0) and upper level (col1), among other data. It will also store all of the final rates which are loaded into the rate matrix
    rel_levs = np.zeros((len(levs[:,0]),3))
    for i in range(0,len(rel_levs[:,0])):
        rel_levs[i,0] = i
    #We go through and match transition energies to those in NIST
    #At the same time, we keep track of which levels have transitions w/ A vals and generate the rates (product of einstein coeffs. and rad. field)
    #If a level has a transition, we store a value '1' in the second column of "rel_levs" array
    idx = np.zeros((len(lines[:,0]),7))
    print('Mapping Internal Indexing Scheme & Generating Rates...')
    datagen_start = time.time()
    for i in range(0,len(lines[:,0])):
        #Check we are not on a header line; otherwise, proceed!
        upper_lev = lines[i,upper_col]
        upper_j = lines[i,upperj_col]
        lower_lev = lines[i,lower_col]
        lower_j = lines[i,lowerj_col]
        #find level indices:
        for j in range(0,len(levs[:,0])):
            #For each line, we have  the upper/lower level indices.
            if ((upper_lev == levs[j,3]) and (float(upper_j) == float(levs[j,2]))):
                idx[i,1] = j 
                rel_levs[j,1] = 1
                upper_energy = levs[j,3]
                #print('Upper Level found for line {:}, index {:}'.format(lines[i,0],i))
            if ((lower_lev == levs[j,3]) and (float(lower_j) == float(levs[j,2]))):
                idx[i,0] = j
                rel_levs[j,1] = 1
                lower_energy = levs[j,3]
                #print('Lower Level found for line {:}, index {:}'.format(lines[i,0],i))
        wavelength_vac = 1e7 / (upper_energy - lower_energy)
        #
        wavelength_air = vac_to_air(wavelength_vac)
        g_i = 2 *float(levs[int(idx[i,0]),2]) + 1 #lower level degeneracy
        g_j = 2 *float(levs[int(idx[i,1]),1]) + 1 #upper level degeneracy
        aval = float(lines[i,aval_col]) #spont. emission rate
        #OLD; 
        stim_coeff = (wavelength_vac * 1e-9)**5 * aval * 1/(8 * np.pi * h * c**2)
        #Temporarily "wrong" to test theory:
#        stim_coeff = 2 * np.pi * (wavelength_vac * 1e-9)**5 * aval * 1/(8 * np.pi * h * c**2)
        #NEW; test:
        #stim_coeff = (wavelength_vac * 1e-9)**5 * aval * 1/(8 * np.pi * h * c**2)
        #
        absorp_coeff = stim_coeff * g_j / g_i     #calc'd absorption coeff. from stim. coeff. This is regardless of form of rad. field. 4pi factor added sep 27
        #Calculate rates from integrated fluxes that were pre-calculated:
        rad_at_wav = fluxes[i]
        absorp_rate = absorp_coeff * rad_at_wav
        stim_rate = stim_coeff * rad_at_wav 
        idx[i,2] = wavelength_air
        idx[i,3] = wavelength_vac
        idx[i,4] = aval
        idx[i,5] = absorp_rate
        idx[i,6] = stim_rate
    datagen_end = time.time()
    #Now that we know which levels have transitions, we generate a mapping from original index to new index scheme:
    #We store that in col 2 of the rel_levs array:
    new_idx = 0
    for i in range(0,len(rel_levs[:,0])):
        if (rel_levs[i,1] == 1):
            rel_levs[i,2] = new_idx
            new_idx = new_idx + 1
    print('Total time generating rates from fluxes: {:} s'.format(round(datagen_end - datagen_start,4)))
    print('{:} Transitions Connecting {:} out of {:} possible levels'.format(len(lines),new_idx,len(levs[:,0])))
    return(idx,new_idx,rel_levs)
#%%
def populate_lhs(lhs_empty,idx,rel_levs):
    """
    SJB
    This function takes in the "idx" array containing all of the rates (emission, absorption, stimulated)
    and loaded them into a matrix which will be solved for the populations.
    """
    construct_start = time.time()
    print('Starting Population of LHS Matrix ...')
    for i in range(0,len(idx[:,0])):
        lower_idx = int(idx[i,0]) #lower level indx for line lines[i,0], belonging to level index idx[i,0] in nist arr
        upper_idx = int(idx[i,1]) #similarly, upper level indx
        #troubleshooting flags:
        low_idx_mapped = int(rel_levs[lower_idx,2])
        up_idx_mapped = int(rel_levs[upper_idx,2])
        aval = float(idx[i,4])
        absorp_rate = float(idx[i,5])
        stim_rate = float(idx[i,6])
        """
        Upper level denoted by j, lower by i:
        Spont. emission, A_j->i values populate two cells: positive term to [i,j], negative to [j,j]
        Stim. emission, B_j->i values populate two cells: positive term to [i,j], negative to [j,j]
        Absorption, B_i->j values populate two cells: negative to [i,i], positive to [j,i]
        """
        #A Value; 2 locations:
        lhs_empty[low_idx_mapped,up_idx_mapped] = lhs_empty[low_idx_mapped,up_idx_mapped] + aval #Correct
        lhs_empty[up_idx_mapped,up_idx_mapped] = lhs_empty[up_idx_mapped,up_idx_mapped] - aval #Correct
        #stimulated emission:
        lhs_empty[low_idx_mapped,up_idx_mapped] = lhs_empty[low_idx_mapped,up_idx_mapped] + stim_rate #Correct
        lhs_empty[up_idx_mapped,up_idx_mapped] = lhs_empty[up_idx_mapped,up_idx_mapped] - stim_rate #Correct
        #absorption:
        lhs_empty[low_idx_mapped,low_idx_mapped] = lhs_empty[low_idx_mapped,low_idx_mapped] - absorp_rate #Correct
        lhs_empty[up_idx_mapped,low_idx_mapped] = lhs_empty[up_idx_mapped,low_idx_mapped] + absorp_rate #Correct
    construct_end = time.time()
    print('Rate Matrices Populated in {:} seconds'.format(round(construct_end - construct_start,3)))
    return(lhs_empty)
#%%
def calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm):
    """
    Updated 04-29-2022
    SJB
    This function is called by "calculate_pops_and_gfactors" and calculates the g-factors for all of the lines in the model.
       g = N * hc/lambda *  Aval
       since hc/lambda is energy, we'll just calculate the energy (in Joules) from the vacuum wavelength
       from E = hc/lamba  
    """
    intens_start_time = time.time()
    print('Starting Fluorescence Efficiencies Calculation...')
    pops = np.real(pops) #Bad matrix solutions can return complex #s. This tries to avoid any problems related to that.
    intens = np.zeros((len(idx[:,0]),3))   
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) g-factor. 
    for i in range(0,len(intens[:,0])):
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        upper_pop = pops[mapped_index] #Grab the (normalized) population
        energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
        aval = idx[i,4]
        line_intens = upper_pop * energy * aval #calculate g factor from above
        intens[i,0] = doppler(idx[i,2],geo_vel) #Save (air) wavelength
        intens[i,1] = doppler(idx[i,3],geo_vel) #Save vac wavelength
        intens[i,2] = line_intens #g factor in J/s/mol
    #Normalize to max line = 1 if desired:
    if (renorm == True):
        max_intens = np.nanmax(intens[:,2])
        intens[:,2] = intens[:,2] / max_intens
    intens_end_time = time.time()
    print('Efficiencies calculated in {:} seconds'.format(round(intens_end_time - intens_start_time,3)))
    return(intens)


def calculate_pops_and_gfactors(flor,element_string,orbit_id,renorm=False):  
    """
    SJB 06-14-2021
    This function encompasses the actual execution of the model. It generates the internal indexing scheme, 
    loads the rate matrix & solves it for the populations, and calculates the g-factors.

    if renorm == False, the gfactors are not renormalized. 
    if renorm == True, then the largest g-factor is set to 1. 

    """
    #Define output dict structure:
    flor[element_string][orbit_id]['outputs'] = {}

    #Define arrays for lines/levels after we process them:
    #Generate the internal indexing scheme:    
    lines = flor[element_string][orbit_id]['relevant_lines']
    levs = flor[element_string]['levs_float']
    aval_col = flor[element_string]['column_indices']['aval_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    ritz_col = flor[element_string]['column_indices']['ritz_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    lowerj_col = flor[element_string]['column_indices']['lowerj_col'] #Added on Sep 23, 2022
    upperj_col = flor[element_string]['column_indices']['upperj_col'] #Added on Sep 23, 2022
    geo_vel = flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] * 1000 #m/s required for doppler formula
    fluxes = flor[element_string][orbit_id]['fluxes']
    #
    model_indexing = generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col, lowerj_col, upperj_col)
    idx = model_indexing[0] # The "idx" array contains all of the line information (indices, wavelengths, and 3 atomic rates) needed as input
    new_idx = model_indexing[1] #Number of levels involved in the model
    rel_levs = model_indexing[2] #Relevant levels for the model; this contains the original level indexing, an indicator if a level is in the final model, and the internal indexing scheme
    lhs_empty = np.zeros((new_idx,new_idx)) #Matrix A of matrix equation A * [X] = B; here lhs_empty will become matrix A
    lhs = populate_lhs(lhs_empty,idx,rel_levs) #Pass in idx array (containing atomic data) and populate the LHS matrix (matrix A)
    rhs = np.zeros((len(lhs[:,0]),1)) #Generate the 'B' Matrix
    rhs[0] = 1
    #Normalize to total pop = 1; On RHS, we have:
    #rhs[0] = matrix_rescale_factor
    #and enforce in LHS: set entire first row of lhs to 1.
    lhs_copy = np.copy(lhs)
    flor[element_string][orbit_id]['outputs']['lhs_matrix_raw'] = lhs_copy
    #REMOVED NORM FOR TESTING:
    for i in range(0,len(lhs[:,0])):
        lhs[0,i] = 1
    #TEMPORARY:
    #SOLVE Ax=b via x = A^-1 * B
    pops = matrix_solve(lhs,rhs)   #Attempt to solve matrix equation using, in order, matrix inversion, pseudo-inverse, and Singular value decomp.
    #Using the populations, calculate the g factors:
    line_intensities = calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm)
    #Save outputs to the dictionary structure:
    flor[element_string][orbit_id]['outputs']['rates'] = idx
    flor[element_string][orbit_id]['outputs']['relevant_levels'] = rel_levs
    flor[element_string][orbit_id]['outputs']['pops'] = pops
    flor[element_string][orbit_id]['outputs']['gfactors'] = line_intensities   
    flor[element_string][orbit_id]['outputs']['lhs_matrix'] = lhs
    flor[element_string][orbit_id]['outputs']['rhs_matrix'] = rhs
    return(flor)

#%%
def doppler_dist(wav_line, t_comet, m_species):
    """
    SJB
    This function calculates and returns a doppler line profile at temperature (in Kelvin) "t_comet", for a particle with mass "m_species" (in kg)
    The function assumes the profile covers wav_line +/- 0.005 nm. 
    
    Will be modified in the future to auto-detect the necessary wavelength grid. 
    """
    wav_coverage = 0.005
    c = 2.99792458e8  #SI
    kb = 1.3864852e-23 #SI
    center = wav_line #nm.    
    lin = np.linspace(center-wav_coverage, center + wav_coverage,500)
    xdiff = center - lin
    sigma = center * np.sqrt(2 * kb * t_comet / (m_species * c**2))
    exp_term = np.exp(-(xdiff)**2 / (2 * sigma**2))
    norm = scipy.integrate.simps(exp_term, lin) #Normalization calculated to ensure sum over all space is 1
    profile = exp_term / norm
    output = np.empty((500,2))
    #Save the WAVELENGTH grid and line profile in 2 col format:
    output[:,0] = lin
    output[:,1] = profile    
    return(output)

def generate_rad_grid(rad, wavelengths):
    """
    SJB
    This function takes in the radiation field "rad" (2 column; col0 = wavelength, col1 = flux in W/m^3)
    and interpolates the values onto the line profile in variable wavelengths. 
    The wavelengths array (the profile!) has the same format: col0 = wavelength, col1 = rho(lambda)
    """   
    try:
        output_arr = np.empty((len(wavelengths),2))
        interpd = interpolate.interp1d(rad[:,0],rad[:,1])
        output_arr[:,0] = wavelengths
        output_arr[:,1] = interpd(wavelengths)
    except ValueError:
        print('Problem with line profile / radiation field for wavelength {:}'.format(wavelengths))
        pass
    return(output_arr)
#%%
def error_iterations(flor,element_string,orbit_id, bbsupp=True, bbtemp = 5777, lower_cutoff = 1e30, upper_cutoff=1e30, aval_min = 0, rad_choice = 'default', profiles = 'doppler', err_id = 'default_setting', num_iterations=100,process=True):
    """
    08-10-2021
    SJB
    This functions takes in the 'flor' dict and iterates the model to estimate uncertainties in the output from uncertainties in transition rates.
    
    Required Inputs:
        1. 'flor' dict: see elsewhere.
        2. 'element_string': see elsewhere.
        2. 'orbit_id'': string indicating which dict of orbital conditions, radiation choice etc. to use in iterations.
        3. 'rad_choice': See documentation elsewhere; options are 'Default' (high-res compiled spectrum provided), or user-provided radiation fields.
        4. 'profiles': line profiles; options are delta-profile (option 'delta'), or doppler profiles 'doppler'. Temperature for the doppler profile
            is taken from the orbit data in orbit_id dict.
        5. To be continued later.
    
    Optional Inputs:
        1. num_iterations: number of times to iterate model. Recommended size is ~10^3 - 10^4 or greater. Default is 100 for testing purposes.
        3. 'process': Bool; defaults to True. If True, the raw iteration data are converted to 4 column format in
            flor[orbit_id][error_calcs][err_id]['proccesed_iteration_output']
            with the following columns: Col 0 (Wave Air), Col1 (Wave Vac), Col2 g-factor, Col3 standard deviation of g-factor. 
            Thus, the g-factors would be plottable as Col2 +/- Col3. 
    """
    print('Starting Error Calculation Run #1...')
    #First, we need to run the model once to get the fluxes and line data.
    flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = rad_choice, profiles=profiles, bbtemp=bbtemp, bbsupp=bbsupp,lower_cutoff=lower_cutoff,upper_cutoff=upper_cutoff,aval_min=aval_min)
    flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)
    #
    first_output_dict = flor[element_string][orbit_id]['outputs'].copy()
    
    print('Completed. Defining output array dimensions...')
    if ('error_calcs' not in flor[element_string][orbit_id]):
        flor[element_string][orbit_id]['error_calcs'] = {}
    flor[element_string][orbit_id]['error_calcs'][err_id] = {}
    flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'] = np.zeros((len(first_output_dict['gfactors'][:,0]),3+num_iterations))
    #
    flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][:,0] = flor[element_string][orbit_id]['outputs']['gfactors'][:,0]
    flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][:,1] = flor[element_string][orbit_id]['outputs']['gfactors'][:,1]
    flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][:,2] = flor[element_string][orbit_id]['outputs']['gfactors'][:,2]
    #
    print('Completed. Starting {:} Iteration Run.'.format(num_iterations))
    #Grab the lines array, uncertainty column, and transition rate column:
    original_lines_arr = flor[element_string][orbit_id]['relevant_lines']
    aval_col = flor[element_string]['column_indices']['aval_col']
    uncert_col = flor[element_string]['column_indices']['uncert_col']
    #Make a copy of the lines array to work with; we'll replace the 'relevant_lines' array at the end of the iterations w/ the original so the 
    #original transition rates are preserved.
    #Grab the original A values:
    original_avals = np.copy(original_lines_arr[:,aval_col])
    uncerts = original_lines_arr[:,uncert_col]    
    #Now we can iterate!
    for i in range(0,num_iterations):
        #Reload original A values:
        #Generate new avals:
        new_avals = generate_new_avals_float(original_avals,uncerts)
        #The new_avals has dimension (#lines,1), and we want (#lines). Reshape:
        flor[element_string][orbit_id]['relevant_lines'][:,aval_col] = np.reshape(new_avals,len(new_avals))
        #Re-calc pops & g-factors:
        flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)
        flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][:,3+i] = flor[element_string][orbit_id]['outputs']['gfactors'][:,2]        
    #Now, we re-save over the original output array to make things easy:
    flor[element_string][orbit_id]['relevant_lines'][:,aval_col] = np.reshape(original_avals,len(original_avals))
    flor[element_string][orbit_id]['outputs'] = first_output_dict    
    #Now, we just process the data if we want that option:
    print('Iterations Completed and Output Saved.')
    if (process == True):
        print('Generating min, max, and stdev for each line...')
        #Define 'processed' array:
        flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'] = np.zeros((len(flor[element_string][orbit_id]['outputs']['gfactors'][:,0]),6))
        flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][:,0] = flor[element_string][orbit_id]['outputs']['gfactors'][:,0]
        flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][:,1] = flor[element_string][orbit_id]['outputs']['gfactors'][:,1]
        flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][:,2] = flor[element_string][orbit_id]['outputs']['gfactors'][:,2]
        for i in range(0,len(flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'])):
            #Calculate the MINIMUM value for a given line:
            flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][i,3] = np.amin(flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][i,2:])
            #Calculate the MAXIMUM value for a given line:
            flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][i,4] = np.amax(flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][i,2:])
            #Calculate the STANDARD DEVIATION for a given line:
            flor[element_string][orbit_id]['error_calcs'][err_id]['processed_iteration_output'][i,5] = np.std(flor[element_string][orbit_id]['error_calcs'][err_id]['raw_iteration_output'][i,2:]) 
    return(flor)



#%%

def calculate_pops_and_gfactors_time_dep(flor,element_string,orbit_id,td_id,time_grid,init_conds=[],renorm=False,time_integrated=False):  
    """
    SJB
    Time-dependent version of calculate_pops_and_gfactors.    
    
    td_id == unique dict key that defines the run.
    
    time_grid == 1d grid of times to solve at. a np.geomspace spanning several orders of magnitude, with the maximum value large enough to be in
    fluorescence equilibrium is recommended
    
    init_conds == array of initial population. Default is all population in ground state.    
    Example of user defined initial conditions for a model with say 30 levels. Note that the population in the initial conditions MUST sum to 1.
            initial_conditions = np.zeros((30,1))
            initial_conditions[0] = 0.5
            initial_conditions[1] = 0.25
            initial_conditions[2] = 0.25
    
    """
    #Define output dicts:
    if (td_id not in flor[element_string][orbit_id]):
        flor[element_string][orbit_id][td_id] = {}
    #Define arrays for lines/levels after we process them:
    #Generate the internal indexing scheme:
    lines = flor[element_string][orbit_id]['relevant_lines']
    levs = flor[element_string]['levs_float']
    aval_col = flor[element_string]['column_indices']['aval_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    ritz_col = flor[element_string]['column_indices']['ritz_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    lowerj_col = flor[element_string]['column_indices']['lowerj_col'] #Added on Sep 23, 2022
    upperj_col = flor[element_string]['column_indices']['upperj_col'] #Added on Sep 23, 2022
    #    geo_vel = flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] * 1000 #m/s required for doppler formula
    fluxes = flor[element_string][orbit_id]['fluxes']
    #
    model_indexing = generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col,lowerj_col,upperj_col)
    idx = model_indexing[0] # The "idx" array contains all of the line information (indices, wavelengths, and 3 atomic rates) needed as input
    new_idx = model_indexing[1] #Number of levels involved in the model
    rel_levs = model_indexing[2] #Relevant levels for the model; this contains the original level indexing, an indicator if a level is in the final model, and the internal indexing scheme
    lhs_empty = np.zeros((new_idx,new_idx)) #Matrix A of matrix equation A * [X] = B; here lhs_empty will become matrix A
    lhs = populate_lhs(lhs_empty,idx,rel_levs) #Pass in idx array (containing atomic data) and populate the LHS matrix (matrix A)
    rhs = np.zeros((len(lhs[:,0]),1)) #Generate the 'B' Matrix
    #Normalize to total pop = 1; On RHS, we have:
    if (len(init_conds) ==0):
        rhs[0] = 1
    else: #If true, the user has provided their own initial conditions.
        rhs[:] = init_conds[:]
        #Save outputs:    
    flor[element_string][orbit_id][td_id]['rates'] = idx
    flor[element_string][orbit_id][td_id]['relevant_levels'] = rel_levs
    flor[element_string][orbit_id][td_id]['lhs_matrix'] = lhs
    flor[element_string][orbit_id][td_id]['rhs_matrix'] = rhs #initial conditions saved.
    flor = time_dep_solve(flor,element_string,orbit_id,lhs,rhs, td_id=td_id,time_grid=time_grid)   
    #Calculate time-dependent g-factors:
    flor = calc_line_intensities_time_dep(flor,element_string,orbit_id,td_id,time_int=time_integrated)
    #Done.
    return(flor)
#%
def time_dep_solve(flor,element_string,orbit_id,lhs,rhs,td_id,time_grid):
    """
    09-22-2021
    Time-dependent solution using matrix exponentiation method.
    Recommended to use np.geomspace to generate the time grid that is ultimately passed here.
    
    Our equation takes the form dx/dt = Ax(t=0)   where x(t=0) are the initial conditions.
    If the user doesnt pass in initial conditions, it is assumed that all initial population is in the ground state.
   
    Eq:
    x(t) = exp(t * A) * x(t=0)
    
    If A is diagonalizable, it can be written as A = PDP^-1
    
    Thus,
    x(t) = P*[diagonal matrix with exp(eigenvalue * t) as diagonals] * P^-1 * x(t=0)
    where P = | eigvec1   eigvec2   ... |
              | eigvec1   eigvec2   ... |
              |   ...       ...     ... |  
    """
    td_start = time.time()
    #Calculate eigenvalues / eigenvectors:
    eigenvals, eigenvectors = np.linalg.eig(lhs)
    #Calculate inverse of eigenvector matrix:
    eigenvectors_inverse = np.linalg.inv(eigenvectors)
    #Output array
    pops_t = np.zeros((len(rhs),len(time_grid)))
    pops_t[:,0] = rhs.reshape(-1)
    for i in range(0,len(time_grid[:])):
        diag_eig_vals = np.diag(np.exp(eigenvals * time_grid[i]))
        mult1 = np.matmul(eigenvectors,diag_eig_vals)
        mult2 = np.matmul(mult1,eigenvectors_inverse)
        pops = np.dot(mult2,rhs)
        pops_t[:,i] = pops.reshape(-1)
    flor[element_string][orbit_id][td_id]['time_grid'] = time_grid
    flor[element_string][orbit_id][td_id]['time_dep_pops'] = pops_t
    td_end = time.time()
    print('Time-dependent solution completed. Elapsed: {:} seconds'.format( round( abs(td_end - td_start),2)))
    return(flor)    
#%%
def calc_line_intensities_time_dep(flor,element_string,orbit_id,td_id,time_int=False):
    """
    SJB
    09-22-2021
    Time-dependent calculation of line intensities.
    Modified on Feb 7, 2023
    The user-provided optional parameter "time_int" was added.
    If set to True (default), the time-dependent g-factors are calculated in two ways:
        
        1. The first way is to first produce the instantaneous time-dependent gfactors as simply n(t) * hc/lambda. This occurs for all 
        values (either True or False) for time_int
        2. If time_int=True, the time-integrated g-factors from t=0 to time t = t_j will be computed. As expected, this step takes a very long
        time as (# transitions x # timesteps) integrations need to be carried out.

    In the interest of saving needless computation time, the default behavior is false even though the integrated quantity is more likely
    the quantity of interest.
    
    If only a few time-integrated quantities are desired, a separate function "gen_gfacs_time_integrated" is provided (see function header for use)
    """
    intens_start_time = time.time()
    print('Starting Time-Dependent Fluorescence Efficiencies Calculation...')
    #Grab the rate information array:
    idx = flor[element_string][orbit_id][td_id]['rates']
    rel_levs = flor[element_string][orbit_id][td_id]['relevant_levels']
    td_pops = flor[element_string][orbit_id][td_id]['time_dep_pops']
    geo_vel = flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] * 1000 #m/s required for doppler formula
    #Grab the time grid:
    time_grid = flor[element_string][orbit_id][td_id]['time_grid']
    td_ints = np.zeros((len(idx[:,0]),len(time_grid)))    
    waves = np.empty((len(idx[:,0]),2))
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) intensity (arb. units.)
    for i in range(0,len(td_ints[:,0])):
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        waves[i,0] = doppler(idx[i,2],geo_vel) #Save (air) wavelength
        waves[i,1] = doppler(idx[i,3],geo_vel) #Save vac wavelength
        energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
        aval = idx[i,4]
        #To compute the instantaneous time-dependent g-factors, we don't need to loop over the time grid:
        td_ints[i,:] = td_pops[mapped_index,:] * energy * aval
    #If the user is interested in the time-integrated g-factors, this is a bit more involved. We loop over the time array and integrate accordingly
    #Note: Be prepared to wait some time. Large structures with many transitions can take up to several hours for large time grids (~10,000 steps)
    #I'm sure there is  way to make this faster, but a faster approach was not needed for my work at the time this was written.
    if (time_int == True):
        #Define the output array:
        time_integrated_gfacs = np.zeros((len(idx[:,0]),len(time_grid)))    
        for i in range(0,len(td_ints[:,0])):
            upper_pop_index_unmapped = int(idx[i,1])
            mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
            #index that we stored in the rel_levs array.
            waves[i,0] = doppler(idx[i,2],geo_vel) #Save (air) wavelength
            waves[i,1] = doppler(idx[i,3],geo_vel) #Save vac wavelength
            energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
            aval = idx[i,4]   
            #Presaving the first element so we don't have to perform an extra if check for each transition:
            dt = time_grid[0]
            pop = dt * td_pops[mapped_index,0]
            time_integrated_gfacs[i,0] = pop * energy * aval
            for j in range(1,len(time_grid[:])):
                #Modified to time-integral on 04-20-2022
                #Fixed index issue on Feb 7, 2023
                time_integrated_gfacs[i,j] = scipy.integrate.simps(td_pops[mapped_index,:j+1], time_grid[:j+1]) * energy * aval / time_grid[j] #The joys of python indexing! This bug took awhile to spot...
    
        flor[element_string][orbit_id][td_id]['time_integrated_gfactors'] = time_integrated_gfacs 
    intens_end_time = time.time()
    flor[element_string][orbit_id][td_id]['time_dep_gfactors'] = td_ints 
    flor[element_string][orbit_id][td_id]['waves'] = waves    
    print('Efficiencies calculated in {:} seconds'.format(round(intens_end_time - intens_start_time,3)))
    return(flor)

def calculate_pops_gfactors_lte(flor,element_string,orbit_id,orbital_conds, renorm=False, temperature=None):
    """
    SJB
    08-15-2022
    This function uses the apparatus' of FlorPy to generate populations and gfactors for LTE level populations. User-provided temperature.
    If a temperature is not provided, the temperature defaults to the blackbody temperature at the heliocentric distance given by:
        
        T = 297 * 1/sqrt(r) | r in AU
    
    flor: input dict after loading atomic data and radiation fields
    element_string: specifies which element is being calculated.
    orbital_conds: list of orbital parameters; distance in AU, velocity in km/s
    
    In LTE, the populations are determined by the degeneracies and the boltzman factor ONLY. 
    
    Output is stored in flor[element_string][orbit_id]['lte']
    
    Some existing code has been modified or re-used to take advantage of existing functionality. For example, the gfactor code, index-mapping, and others
    were used with "dummy" arrays in some cases so that we don't have to reinvent the wheel
    
    To use the LTE code, the user must run, in order:
        build_model
        load_nist_data
        load_orbit_data
        this function
    """
    k_b = 8.617333262e-5 #eV / K
    cm_inv_to_ev = 1 / 8065.544 #1 cm inverse = 1.23981e-4 eV    Used for computing boltzmann factors
    #
    print('Starting calculation of LTE populations')
    solar_dist = flor[element_string][orbit_id]['orbit_params']['helio_dist_au']
    if (temperature == None):
        print('No LTE temperature specified. Using default blackbody temperature')
        temperature = 279 * (1/np.sqrt(solar_dist)) #Kelvin; comet temperature for doppler profiles
    #
    flor[element_string][orbit_id]['lte'] = {}
    flor[element_string][orbit_id]['lte']['temperature'] = str(temperature) + str(' Kelvin')
    flor[element_string][orbit_id]['lte']['relevant_lines'] = np.copy(flor[element_string]['lines_data_float'])
    flor[element_string][orbit_id]['lte']['outputs'] = {}  
    lines = flor[element_string][orbit_id]['lte']['relevant_lines']
    levs = flor[element_string]['levs_float']
    aval_col = flor[element_string]['column_indices']['aval_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    ritz_col = flor[element_string]['column_indices']['ritz_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    lowerj_col = flor[element_string]['column_indices']['lowerj_col']
    upperj_col = flor[element_string]['column_indices']['upperj_col']
    geo_vel = flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] * 1000 #m/s required for doppler formula
    #We make a dummy "fluxes" array so we can use the same generate_index_scheme code we have from before. This will give us a mapping from transitions -> levels we can use
    #for calculations the gfactors, since the gfactor function will read that in. 
    fluxes = np.zeros((len(lines[:,0]),1))
    model_indexing = generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col, lowerj_col, upperj_col)
    idx = model_indexing[0] # The "idx" array contains all of the line information (indices, wavelengths, and 3 atomic rates) needed as input
    new_idx = model_indexing[1] #Number of levels involved in the model
    rel_levs = model_indexing[2] #Relevant levels for the model; this contains the original level indexing, an indicator if a level is in the final model, and the internal indexing scheme
    #Dummy array for populations:
    pops = np.zeros((len(levs[:,0]),1))
    pops[0] = 1 #Set ground to 1
    ground_degen = 2 * levs[0,2] + 1
    kt = k_b * temperature
    lte_pop_start = time.time()
    print('Computing LTE Populations...')
    #We compute the populations with respect to the ground state:
    for i in range(0,len(levs[:,0])):
        degen = 2 * levs[i,2] + 1
        lev_energy_ev = levs[i,3] * cm_inv_to_ev
        pops[i] = (degen / ground_degen) * np.exp(-1 * (lev_energy_ev) / kt)
    #Re-normalize such that total population is 1:
    pops[:] = pops[:] / np.sum(pops[:])
    lte_pop_end = time.time()
    print('LTE Populations Calculated in {:} seconds'.format(round((lte_pop_end - lte_pop_start)/60,3)))
    print('Calculating LTE G-factors at temperature T = {:} Kelvin'.format(temperature))
    line_intensities = calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm)
    #Save outputs to the dictionary structure:
    flor[element_string][orbit_id]['lte']['outputs']['rates'] = idx
    flor[element_string][orbit_id]['lte']['outputs']['relevant_levels'] = rel_levs
    flor[element_string][orbit_id]['lte']['outputs']['pops'] = pops
    flor[element_string][orbit_id]['lte']['outputs']['gfactors'] = line_intensities   
    return(flor)

#%
def gen_gfacs_time_integrated(flor,element_string,orbit_id,td_id,input_inds):
    """
    SJB
    Feb 7 2023
    This is a custom function intended to allow the user to only calculate time-integrated fluorescence efficiencies for lines they want. This is necessary in some cases
    as calculating the time-integrated gfactors for ALL lines is very computationally demanding.
    This function will will compute the time-integrated gfactors for lines with indices in the input array "input_inds". These indices correspond to the same index scheme in the 
    "relevant_lines_str" dictionary structure and the output equilibrium gfactors
    """
    #Need h and c physical constants in SI units:
    h = 6.62607004e-34
    c = 299792458
    idx = flor[element_string][orbit_id][td_id]['rates']
    rel_levs = flor[element_string][orbit_id][td_id]['relevant_levels']
    td_pops = flor[element_string][orbit_id][td_id]['time_dep_pops']
    #Grab the time grid:
    time_grid = flor[element_string][orbit_id][td_id]['time_grid']
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) intensity (arb. units.)
    time_integrated_gfacs = np.zeros((len(input_inds[:]),len(time_grid)))    
    for i in range(0,len(input_inds[:])):
        upper_pop_index_unmapped = int(idx[input_inds[i],1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        energy = h * c / (idx[input_inds[i],3] * 1e-9) #Calculate transition energy from VAC wavelength.
        aval = idx[input_inds[i],4]   
        #Presaving the first element so we don't have to perform an extra if check for each transition:
        dt = time_grid[0]
        pop = dt * td_pops[mapped_index,0]
        time_integrated_gfacs[i,0] = pop * energy * aval
        for j in range(1,len(time_grid[:])):
            #Modified to time-integral on 04-20-2022
            #Fixed index issue on Feb 7, 2023
            time_integrated_gfacs[i,j] = scipy.integrate.simps(td_pops[mapped_index,:j+1], time_grid[:j+1]) * energy * aval / time_grid[j] #The joys of python indexing! This bug took awhile to spot...
    return(time_integrated_gfacs)


