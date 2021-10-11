# -*- coding: utf-8 -*-
"""
Source code for "FlorPy" fluorescence model.
by Steven J. Bromley at Auburn University
Contact: sjb0068@auburn.edu
The dictionary version of this code was released on 10-11-2021 on GitHub
    

    This program is distributed under the GNU General Public License.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, version 3 of the GNU General Public License is available at
    https://www.gnu.org/licenses/gpl-3.0.en.html
    
    
You, the downloader and potential users of this code, are free to use any and all parts of this code in other programs,
or as stand-alone code. Given the complexity of the below codes, I would appreciate if those who intend to use or 
modify this code reach out to me for help if needed. I am happy to help or instruct others 
that wish to use, modify, and build on this work.

Best wishes,
Steven J. Bromley

------------------------
"""
##############################
#          Imports           #
##############################

import numpy as np #The OG
import matplotlib.pyplot as plt
import time
import math
import sys
import scipy
import os
from random import random 
from scipy.integrate import quad
from scipy import interpolate

##############################
#          Constants         #
##############################

# Some useful constants All required units are defined
#inside the functions where they are used
Rsun = 6.957e8 # meters
h = 6.62607004e-34
hbar = 1.054571817e-34
k = 1.38064852e-23
c = 299792458
cm_to_J = 1.986491855e-23


##############################
#         Functions          #
##############################

def add_element_to_model(flor,element_string):
    """
    SJB
    This function takes in the flor dict and adds to it a dict for a new element with key "element_string"
    i.e. flor['fe0'] for neutral iron. The key of the dict is added in as a string;
    
    Note there is no requirement on the naming convention; if I add in silicon data but name it iron,
    nothing breaks. The element_string is a labeling convention ONLY and it is best if YOU understand what each entry into flor contains.
    
    Only one entry with a given key can exist. 
    """
    if (element_string not in flor):
        flor[element_string] = {}
    else:
        print('Element key "{:}" already exists in dictionary! Ignoring addition request'.format(element_string))
    return(flor)

#%%
def blackbody_wave(wavelength, temp=5777):
    #Adjusted on 05-25-2021 to accomodate changes made to flux calc. 
    #temp in Kelvin, all in SI
    #This function is not a typical blackbody; it reads in a wavelength (in nm), and outputs
    #Output intensity is W/m^3 
    #Assuming defined in nm, we convert:
    wavelength = wavelength * 1e-9
    h = 6.62607004e-34
    c = 299792458
    k = 1.38064852e-23
    exp_factor = h*c / (wavelength * k * temp) #compute item inside exponential factor
    exp_term = 1/(np.exp(exp_factor) - 1) #calculate exponential. Then, calculate BB intensity (W/m^2 per wavelength(m) per ster.):
    intensity = (2 * h * c**2 / (wavelength**5) )* exp_term
    intensity = np.pi * intensity
    #The factor of pi corrects for solid angle; 
    #Above function integrates to 1358 W/m^2 when integrated and projected to 1AU by (Rsun (m) / 1 AU (m))^2    
    return intensity

#%%
def build_model():
    #06-11-2021 SJB Initiates Dict structure for flor model
    flor = {}
    return(flor)

#%%
def calculate_pops_and_gfactors(flor,element_string,orbit_id,renorm=False,minor_corr=False):  
    """
    This function calculates the equuilibrium populations of levels assuming transition, absorption, and stimulated emission rates,
    and then calculates fluorescence efficiencies, i.e. "g-factors"
    Inputs:
        flor: dictionary structure. It is assumed that the following functions were called before this function:
            build_model, add_element_to_model, load_nist_data, define_rad_fields, load_orbit_data, generate_fluxes_general        
        element string, orbit id: identifiers (keys) for the element and orbital parameters.
    Optional parameters:
        renorm: if "True", the g-factors are normalized such that the strongest line is 1. Beneficial for quick plotting and scaling etc.
        minor_corr: old optional parameter to test sensitivity of matrix solvers to perturbations. Currently does nothing and may be
            re-implemented in a future release.

    """
    #Define arrays for lines/levels after we process them:
    #Generate the internal indexing scheme:
    lines = flor[element_string][orbit_id]['relevant_lines']
    levs = flor[element_string]['levs_float']
    aval_col = flor[element_string]['column_indices']['aval_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    ritz_col = flor[element_string]['column_indices']['ritz_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    geo_vel = flor[element_string][orbit_id]['orbit_params']['geo_vel_kms'] * 1000 #m/s required for doppler formula
    fluxes = flor[element_string][orbit_id]['fluxes']
    #Generate the internal index mapping scheme.
    model_indexing = generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col)
    idx = model_indexing[0] # The "idx" array contains all of the line information (indices, wavelengths, and 3 atomic rates) needed as input
    new_idx = model_indexing[1] #Number of levels involved in the model
    rel_levs = model_indexing[2] #Relevant levels for the model; this contains the original level indexing, an indicator if a level is in the final model, and the internal indexing scheme
    lhs_empty = np.zeros((new_idx,new_idx)) #Matrix A of matrix equation A * [X] = B; here lhs_empty will become matrix A
    lhs = populate_lhs(lhs_empty,idx,rel_levs) #Pass in idx array (containing atomic data) and populate the LHS matrix (matrix A)
    rhs = np.zeros((len(lhs[:,0]),1)) #Generate the 'B' Matrix
    #Normalize to total pop = 1; On RHS, we have:
    rhs[0] = 1
    #and enforce in LHS: set entire first row of lhs to 1.
    for i in range(0,len(lhs[0,:])):
        lhs[0,i] = 1
    #SOLVE Ax=b via x = A^-1 * B
    #THE FOLLOWING IS TEMPORARY FOR TESTING:
    #We try to first solve using standard methods. If necessary, we try the pseudo-inverse and Singular Value Decomposition:
    model_solution = matrix_solve(lhs,rhs)   
    solution_flag = model_solution[0] 
    #If the solution flag is 0, we have a valid solution. If = 1, we have a bad solution and need to revisit the atomic data. 
    if (solution_flag == 0):
        print('Solution is valid and non-negative.')
        pops = model_solution[1]
        #ground_norm = pops[0]
        #pops = pops / np.sum(pops)
    else:
        print('Problem with solution (levels with 0 population). Be wary of fluorescence efficiencies.')
        pops = model_solution[1]
    #All said and done, we now calculate line intensities:
    line_intensities = calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm)
    #Done.
    #Save outputs:    
    flor[element_string][orbit_id]['outputs'] = {}
    flor[element_string][orbit_id]['outputs']['rates'] = idx
    flor[element_string][orbit_id]['outputs']['relevant_levels'] = rel_levs
    flor[element_string][orbit_id]['outputs']['pops'] = pops
    flor[element_string][orbit_id]['outputs']['gfactors'] = line_intensities   
    flor[element_string][orbit_id]['outputs']['lhs_matrix'] = lhs
    flor[element_string][orbit_id]['outputs']['rhs_matrix'] = rhs
    return(flor)

def calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm):
    """
    Called by "calculate_pops_and_gfactors"
    
    We compute the g-factor as (see Bromley et al 2021 PSJ, arxiv):
       g = N * hc/lambda *  Aval
       since hc/lambda is energy, we'll just calculate the energy (in Joules) from the vacuum wavelength
       from E = hc/lamba  
    """
    intens_start_time = time.time()
    print('Starting Fluorescence Efficiencies Calculation...')
    pops = np.real(pops) #Just in case; some populations may be complex depending on which method is used to compute them. This is an insurance policy for stability
    intens = np.zeros((len(idx[:,0]),3))   
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) g-factor.
    for i in range(0,len(intens[:,0])):
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        upper_pop = pops[mapped_index] #Grab the (normalized) population
        energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
        aval = idx[i,4]
        line_intens = upper_pop * energy * aval #calculate line intensity from above
        intens[i,0] = doppler(idx[i,2],geo_vel) #Save (air) wavelength
        intens[i,1] = doppler(idx[i,3],geo_vel) #Save vac wavelength
        intens[i,2] = line_intens #save intensity
    #Normalize to max line = 1 if desired:
    if (renorm == True):
        max_intens = np.nanmax(intens[:,2])
        intens[:,2] = intens[:,2] / max_intens
    intens_end_time = time.time()
    print('Efficiencies calculated in {:} seconds'.format(round(intens_end_time - intens_start_time,3)))
    return(intens)


def convert_to_float(frac_str):
    """
    Updated 06-11-2021 by SJB
    This converts a string to float; it is intended to do 2 things:
     1. Remove flag characters from "questionable" levels in NIST level files, and 
     2. Convert half-integer J's, e.g. 3/2, to decimals (1.5)
    Known questionable flags (to date) are: 'a', brackets i.e. '[ ]', '*', and '?'
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
def define_rad_fields(flor,element_string,rad_source,bbtemp = 5777, radfiles=[], vac_or_air = []):
    """
    SJB
    06-11-2021
    This function defines the radiation field(s) used. 
    
    IF rad_source = 'default', the following spectra is used with all wavelengths in VACUUM NM:
          0 - 168 nm: We use the "13x1" corona / chromosphere model of Fontenla et al 2014, which has good agreement from 168 down
          with data from the SORCE mission.
          168 - 200.06 nm: SOLSTICE spectral measurements with 0.1 nm resolution
          200.06 - 202: Hall & Anderson (1991) UV measured spectra
          202 - 2730 nm: Composite high-resolution benchmarked solar spectrum from Coddington et al (2021). Their dataset is compiled
          from various space- and ground-based measurements of the solar spectrum, including the famous measurements of Kurucz.
      The provided "composite" spectrum (filename "0-2730_vac.txt") contains the above spectra in 2 column format.
            
    A blackbody choice is available. Default temp is  5777 (~solar). User is free to set any temp they wish.

    IF rad_source = 'user', the user must provide two lists to this function call:
        1. "radfiles" variable: list of files containing the spectra to be used. The files must be in 2 column format, with
            Col 0 = Wavelength (air or vacuum), and Col 1 = Watts / m^2 / s
        2. "medium" variable: a list of media for the above radfiles. 
        
     
        For example, if the user provided 2 files called spec1 and spec2:
            spec1: 0 - 300 nm (vac), and 
            spec2: 500 - 2000 nm (air), 
        The inputs 'radfiles' and vac_or_air would look like:
                radfiles = ['spec1.txt', 'spec2.txt']
                vac_or_air = ['vac', 'air']
           
        For points outside these ranges, the user can either set rates to 0 (not recommended; leads to potentially bad solutions) by setting bbsupp = False later,
        or use a blackbody to estimate fluxes at wavelengths outside those bounds; for the above case, this would be for lines between 300 - vac_to_air(500),
        and above 2000. Using a blackbody to estimate fluxes would be setting the variable bbsupp = True later.
    """    
    if ('rad_choices' not in flor[element_string]):
        flor[element_string]['rad_choices'] = {}
    if (rad_source == 'default'):
        print('Default Radiation fields selected. Loading data...')
        flor[element_string]['rad_choices']['default'] = {}
        flor[element_string]['rad_choices']['default']['all'] = {}
        flor[element_string]['rad_choices']['default']['all']['medium'] = 'vac'
        flor[element_string]['rad_choices']['default']['all']['lims'] = np.zeros((1,2))
        flor[element_string]['rad_choices']['default']['all']['fluxes'] = np.genfromtxt('0-2730_vac.txt', delimiter='\t')
        flor[element_string]['rad_choices']['default']['all']['lims'][0,0] = flor[element_string]['rad_choices']['default']['all']['fluxes'][0,0]
        flor[element_string]['rad_choices']['default']['all']['lims'][0,1] = flor[element_string]['rad_choices']['default']['all']['fluxes'][-1,0]
        flor[element_string]['rad_choices']['blackbody'] = {}
        flor[element_string]['rad_choices']['blackbody']['temperature_K'] = bbtemp
    print ('Default solar spectrum loaded. See documentation for spectral details.')
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
    return(flor)

#%%
def doppler(wavelength,velocity):
    #positive velocity defined as away from sun
    #Updated 06-21-2021 to relativist form for completeness
    #Wavelength in nm; velocity in m/s. All calls in the model convert from km/s to m/s before calling this function.
    c = 299792458
    vel = velocity #assuming m/s
    rest_wavelength = wavelength #nm
    frac = np.sqrt((1 + velocity/c)/(1 - velocity/c))
    shifted_wavelength = rest_wavelength * frac
    return(shifted_wavelength)

def doppler_dist(wav_line, t_comet, m_species):
    """
    Normalized Gaussian Line profile function; typically called internally. Can be called with the inputs:
        
        wav_line == center wavelength
        t_comet == temperature of ABSORBING SPECIES for thermal width
        m-species == species mass in kg
    
    Default wavelength coverage is +/- 0.005 from center wavelength on a 500 pt grid.
    """
    wav_coverage = 0.005
    c = 2.99792458e8  #SI
    c_nm = c * 1e9
    kb = 1.3864852e-23 #SI
    center = wav_line #nm.    
    lin = np.linspace(center-wav_coverage, center + wav_coverage,500)
    xdiff = center - lin
    #OLD sigma:
    #sigma = center * np.sqrt(kb * t_comet / (m_species * c**2))
    #NEW sigma:
    sigma = center * np.sqrt(2 * kb * t_comet / (m_species * c**2))
    exp_term = np.exp(-(xdiff)**2 / (2 * sigma**2))
    norm = numerically_integrate_2col(lin,exp_term)
    profile = exp_term / norm
    output = np.empty((500,2))
    #Save the WAVELENGTH grid and line profile in 2 col format:
    output[:,0] = lin
    output[:,1] = profile    
    #Normalize profile so intregal of profile over all wavelengths = 1 !
    return(output)

def error_iterations(flor,element_string,orbit_id, bbsupp=True, bbtemp = 5777, lower_cutoff = 1e30, upper_cutoff=1e30, aval_min = 0, rad_choice = 'default', profiles = 'doppler', err_id = 'default_setting', num_iterations=100,process=True):
    """
    Monte-Carlo iteration function to generate approximate uncertainties of g-factors calculated with these codes.
    
    This functions takes in the 'flor' dict and iterates the model to estimate uncertainties in the output from uncertainties in transition rates.
    The regular inputs are much the same as the "normal" functions, "generate_fluxes_general" and "calculate_pops_and_gfactors" 
    
    See Bromley et al 2021 in PSJ or Arxiv for further details.
   
    For these Monte-Carlo iterations, the additional parameters are:
        1. num_iterations: number of times to iterate model. Recommended size is ~10^3 - 10^4 or greater. Default is 100.
        3. 'process': Bool; defaults to True. If True, the raw iteration data are converted to 4 column format in
            flor[orbit_id][error_calcs][err_id]['proccesed_iteration_output']
            with the following columns: Col 0 (Wave Air), Col1 (Wave Vac), Col2 g-factor, Col3 standard deviation of g-factor. 
            Thus, the g-factors would be plottable as Col2 +/- Col3. 
    """
    print('Starting Error Calculation: Generating Fluxes and "normal" g-factors...')
    #First, we need to run the model once to get the fluxes and line data.
    flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = rad_choice, profiles=profiles, bbtemp=bbtemp, bbsupp=bbsupp,lower_cutoff=lower_cutoff,upper_cutoff=upper_cutoff,aval_min=aval_min)
    flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)
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
    print('Completed Dimensioning. Starting {:} Iterations...'.format(num_iterations))
    #Grab the lines array, uncertainty column, and transition rate column:
    original_lines_arr = flor[element_string][orbit_id]['relevant_lines']
    aval_col = flor[element_string]['column_indices']['aval_col']
    uncert_col = flor[element_string]['column_indices']['uncert_col']
    #Make a copy of the lines array to work with; we'll replace the 'relevant_lines' array at the end of the iterations w/ the original so the 
    #original transition rates are preserved.
    modified_lines_arr = np.copy(original_lines_arr)
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
def find_nearest_index(array,value):
    diff_temp = abs(value-array[0])
    index = 0 
    for i in range(0,len(array[:])):
        diff = abs(array[i] - value)
        if (diff < diff_temp):
            index = i
            diff_temp = diff
    return index

def fniiv(array,value,start_idx):   
    """
    " Find Nearest Index Iterative Version"
    This is essentially a "find nearest index" code in a sorted array, with the modification that 
    the found index is used to inform the next search. 
    
    It is assumed that the input values being looked for are sorted.
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
def generate_fluxes_general(flor,element_string,orbit_id,rad_choice='default', bbsupp=True, bbtemp = 5777, profiles='delta', lower_cutoff = 1e30, upper_cutoff=1e30, aval_min = 0):
    """
    In the detailed balance of transition, absorption ,and stimulated emission rates (See Bromley et al 2021 in PSJ)
    
    the absorption and stimulated emission rates have two components: the B coefficient, and the integral of the radiation field over a line profile
    This function computes the radiation-field-dependent portions of the rates. This calculation is handled separately as it
    increases the speed of the code when running Monte-Carlo iterations to estimate model uncertanties.
    
    Variables:
        flor: dictionary structure. The following functions must have been called before this for normal execution:
            flor = build_model()
            flor = add_element_to_model(flor, element_string)
            flor = load_nist_data(flor,element_string,lines_file, levs_file)
            flor = define_rad_fields(flor,element_string,rad_source='default')
        bbsup: True or False. If True, a blackbody is used to estimates data for transition outside the scope of the chosen radiation field(s).
        bbtemp: if bbsup is True, this sets the temperature of the blackbody used in bbsup.
        profiles: choice of line profile. Two options are availabe: "delta" for a Dirac Delta, and "doppler" for a doppler-broadened (thermal)
            line profile. Note that the temperature in the line profile is a different temperature than the blackbody temperature of the radiation field

    Users may place restrictions on atomic data as follows:
        lower_cutoff (optional): sets the MAXIMUM allowed energy (in cm^-1) of the LOWER levels of the transitions.
        upper_cutoff (optional): sets the MAXIMUM allowed energy (in cm^-1) of the UPPER levels of the transitions.
        aval_min (optional): Sets the LOWER LIMIT of the transition rates (Note: NOT the B coefficients) included in the model
    Default values for the atomic data restrictions are set such that all atomic data is used by default.
    
    """
    print('Starting Calculation of Integrated Fluxes...')
    start = time.time()
    cm_to_J = 1.986491855e-23
    Rsun = 6.957e8 # meters
    h = 6.62607004e-34
    c = 299792458
    #Grab necessary data from dict structure:
    raw_levs = flor[element_string]['levs_float']
    flor = prep_indexing(flor,element_string,orbit_id,rad_choice, bbsupp, lower_cutoff, upper_cutoff, aval_min)
    raw_lines = flor[element_string][orbit_id]['relevant_lines']
    #Column indexes
    aval_col = flor[element_string]['column_indices']['aval_col']
    lower_col = flor[element_string]['column_indices']['lower_col']
    ritz_col = flor[element_string]['column_indices']['ritz_col']
    uncert_col = flor[element_string]['column_indices']['uncert_col']
    upper_col = flor[element_string]['column_indices']['upper_col']
    integrated_fluxes = np.zeros((len(raw_lines[:,0]),1),dtype=float)
    #Grab orbital data and convert as necessary to SI:
    solar_dist = 1.496e+11 * flor[element_string][orbit_id]['orbit_params']['helio_dist_au'] #meters
    obj_vel = 1e3 * flor[element_string][orbit_id]['orbit_params']['helio_vel_kms'] #m/s
    #We want to loop over the choice of radiation fields and generate fluxes for each line in that range.
    upper_energies = (raw_lines[:,upper_col])
    lower_energies = (raw_lines[:,lower_col])
    vac_waves = 1e9 * h * c / (cm_to_J * (upper_energies - lower_energies)) #Ritz wavelength in VACUUM
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
                    #Added 0.1 here to prevent problems at edge; TEMPORARY
                    if (flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,0] + 1 < doppler(all_waves[i,wave_col],obj_vel) < -1 + flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,1]):
                        dist_scale_1au = (1.496e+11/solar_dist)**2
                        searched_rad_indices = fniiv(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0],doppler(all_waves[i,wave_col],obj_vel),search_idx_iter)
                        #start_idx is passed to next iteration.
                        search_idx_iter = searched_rad_indices[1]
                        #Index for the lines radiation field:
                        rad_idx = searched_rad_indices[0]
                        rad_at_wav = dist_scale_1au * flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][rad_idx,1] * 1e9 #converts from W/m^2 per nm to W /m^2 per meter (wavelength)
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
                    #Added 0.1 here to prevent problems at edge; TEMPORARY
                    if (flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,0] + 1 < doppler(all_waves[i,wave_col],obj_vel) < -1 + flor[element_string]['rad_choices'][rad_choice][key]['lims'][0,1]):
                        profile = doppler_dist(doppler(all_waves[i,wave_col],obj_vel),comet_temp,m_species)
                        dist_scale_1au = (1.496e+11/solar_dist)**2
                        searched_rad_indices = fniiv(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][:,0],doppler(all_waves[i,wave_col],obj_vel),search_idx_iter)
                        rad_idx = searched_rad_indices[0]
                        search_idx_iter = searched_rad_indices[1]
                        #Indices for matching radiation field to line profile:
                        #We take +/- 20 bins of the radiation field:
                        rad_start_idx = searched_rad_indices[0] - 70 #Save the index of the smallest wavelength from profile array
                        rad_end_idx = searched_rad_indices[0] + 70 
                        #print(rad_start_idx)
                        #print(rad_end_idx)
                        new_rad_grid = generate_rad_grid(flor[element_string]['rad_choices'][rad_choice][key]['fluxes'][rad_start_idx:rad_end_idx,:],profile[:,0])
                        rad_at_wav = dist_scale_1au * 1e9 * numerically_integrate_2col(profile[:,0], new_rad_grid[:,1] * profile[:,1])
                        integrated_fluxes[i] = rad_at_wav
            if (bbsupp == True):
                dist_scale  = (Rsun/solar_dist)**2
                for i in range(0,len(integrated_fluxes[:])):
                  if (integrated_fluxes[i] == 0):
                      profile = doppler_dist(doppler(all_waves[i,wave_col],obj_vel),comet_temp,m_species)
                      bb = blackbody_wave(profile[:,0],temp=bbtemp) #Only need the profile wavelengths at temp=5777K (default); assumes nm in
                      rad_at_wav = dist_scale * numerically_integrate_2col(profile[:,0], bb * profile[:,1])
                      integrated_fluxes[i] = rad_at_wav
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
                rad_at_wav = dist_scale * numerically_integrate_2col(profile[:,0], bb * profile[:,1])
                integrated_fluxes[i] = rad_at_wav
            
    end = time.time()
    print('Calculation of Integrated Fluxes finished. Elapsed Time {:} s'.format(round(end - start,2)))
    flor[element_string][orbit_id]['fluxes'] = integrated_fluxes
    return(flor)   

#%%
def generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col):
    """
    SJB
    Used internally by fluorescene_model to generate matrix "idx" array used fto encapsulate all atomic data that goes into the LHS
    matrix of our matrix equation.
    """
    #First, find the LEVELS for the LINES
    # Array 'idx' is the same length as the lines array and will hold the indexes corresponding to the index of the 
    # lower level (col0) and upper level (col1)
    rel_levs = np.zeros((len(levs[:,0]),3)) #removed -1 from len
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
        lower_lev = lines[i,lower_col]
        #find level indices:
        for j in range(0,len(levs[:,0])):
            #For each line, we have  the upper/lower level indices.
            if (upper_lev == levs[j,3]):
                idx[i,1] = j 
                rel_levs[j,1] = 1
                upper_energy = levs[j,3]
                #print('Upper Level found for line {:}, index {:}'.format(lines[i,0],i))
            if (lower_lev == levs[j,3]):
                idx[i,0] = j
                rel_levs[j,1] = 1
                lower_energy = levs[j,3]
                #print('Lower Level found for line {:}, index {:}'.format(lines[i,0],i))
        #Degeneracies and Avalue:
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))
        wavelength_air = vac_to_air(wavelength_vac)
        g_i = 2 *float(levs[int(idx[i,0]),2]) + 1 #lower level degeneracy
        g_j = 2 *float(levs[int(idx[i,1]),1]) + 1 #upper level degeneracy
        aval = float(lines[i,aval_col]) #spont. emission rate
        #
        stim_coeff = (wavelength_vac * 1e-9)**5 * aval * 1/(8 * np.pi * h * c**2)
        absorp_coeff = stim_coeff * g_j / g_i     #calc'd absorption coeff. from stim. coeff. This is regardless of form of rad. field
        #Calculate rates from integrated fluxes that were pre-calculated:
        rad_at_wav = fluxes[i]
        absorp_rate = absorp_coeff * rad_at_wav
        stim_rate = stim_coeff * rad_at_wav #
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
def generate_new_avals_float(avals,uncerts):
    """
    Generates new a values for array of input array values and uncertainties. 
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

#%%
def generate_rad_grid(rad, wavelengths):
    """
    SJB
    Fast scipy interpolator for projecting solar spectrum onto line profile grid. 
    """   
    output_arr = np.empty((len(wavelengths),2))
    interpd = interpolate.interp1d(rad[:,0],rad[:,1])
    output_arr[:,0] = wavelengths
    output_arr[:,1] = interpd(wavelengths)
    return(output_arr)

#%%
def grab_uncertainty(lines,col):
    """
    Added 01-08-2021. This function is used for converting A value uncertainty ratings in NIST to numbers
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
    return(a_uncert)
#%

#%%
def load_nist_data(flor,element_string,lines_file, levels_file):
    """
    SJB
    06-11-2021
    This function reads in NIST lines and levels files from ASD and adds the data to the model.
    
    This function assumes the data is saved in tab-delimited format. Please only use lines with transition rates.
    """
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
                                              'uncert_col' : uncert_col}
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
    #
    flor[element_string]['levs_str'] = raw_levs
    #Added sort to following line on 06-22-2021
    flor[element_string]['levs_float'] = levs[np.argsort(levs[:,3])]
    flor[element_string]['bad_levels'] = bad_levels_out #Saved for later
    print('Lines and levels loaded for model "{:}"'.format(element_string))
    return(flor)



#%%
def load_orbit_data(flor,element_string,orbit_id,orbital_conds, t_comet = 280, m_species = 9.74620109e-26):
    """
    This function loads the orbital data into dictionary structure.
    
    If using a doppler line profile, the default temperature is 280K (~1 AU). The default mass is set to that of nickel in kg.
    
    Each set of distance/velocity is assigned a unique id number, and all calculated data will be saved under that ID number.
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
    return(flor)

#%%
def numerically_integrate(array):
    #Simple function for numerically integrating values using rectangles:
    #Ex: array has columns: (0) wavelength, (1) intensity
    integrated = 0
    #Note that the limits are not included; the limits must be taken care of by ONLY passing in values to be integrated
    for i in range(0,len(array[:,0])-1):
        lower_x = array[i,0]
        upper_x = array[i+1,0]
        diff = upper_x - lower_x
        avg = (array[i,1] + array[i,1])/2
        integrated = integrated + (avg*diff)
    return(integrated)    

#%%
def populate_lhs(lhs_empty,idx,rel_levs):
    """
    SJB
    Added on 03-02-2021 to functionalize construction/generate etc of matrices.
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
        #A Value; 2 locations. lower level:
        lhs_empty[low_idx_mapped,up_idx_mapped] = lhs_empty[low_idx_mapped,up_idx_mapped] + aval #Correct
        lhs_empty[up_idx_mapped,up_idx_mapped] = lhs_empty[up_idx_mapped,up_idx_mapped] - aval #Correct
        #stimulated emission:
        lhs_empty[low_idx_mapped,up_idx_mapped] = lhs_empty[low_idx_mapped,up_idx_mapped] + stim_rate #Correct
        lhs_empty[up_idx_mapped,up_idx_mapped] = lhs_empty[up_idx_mapped,up_idx_mapped] - stim_rate #Correct
        #absorption:
        lhs_empty[low_idx_mapped,low_idx_mapped] = lhs_empty[low_idx_mapped,low_idx_mapped] - absorp_rate #Correct
        lhs_empty[up_idx_mapped,low_idx_mapped] = lhs_empty[up_idx_mapped,low_idx_mapped] + absorp_rate #Correct
    construct_end = time.time()
    print('Rate Matrices Populated in {:} seconds'.format(round(construct_end - construct_start),3))
    return(lhs_empty)

def prep_indexing(flor,element_string,orbit_id,rad_choice, bbsupp, lower_en_cutoff, upper_en_cutoff, aval_min):
    """

    The purpose of this function is to prepare the lines/levels list for ONLY those that are within the limits of the radiation fields.
    If we include the lines without radiation field components and the blackbody suppelementation is applied, we'll have
    a singular matrix and the solution methods will struggle. 
    """
    print('Checking raw line list for relevant Lines.')
    #First, we make a copy of the TOTAL lines file.
    checklines = np.copy(flor[element_string]['lines_data_float'])
    ritz_col = flor[element_string]['column_indices']['ritz_col']
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
            wave_vac = doppler(1e9 * h * c / (cm_to_J * (checklines[i,upper_col] - checklines[i,lower_col])) ,obj_vel) #
            wave_air = doppler(vac_to_air(wave_vac),obj_vel)
            for j in range(0,len(lims[:,0])):
                if ('v' in media[j,0]):
                    wave_to_use = wave_vac
                else:
                    wave_to_use = wave_air
                if ((lims[j,0] < wave_to_use < lims[j,1]) and (checklines[i,upper_col] < upper_en_cutoff) and (checklines[i,lower_col] < lower_en_cutoff) and (checklines[i,aval_col] > aval_min)):  #Check wavelength limits and upper energy cutoff
                    idx_to_remove[i,1] = 1
    #The only lines w/o a '1' to indicate keeping the line are outside the radiation field. We flag those for keeping if bbsupp == True:
    for i in range(0,len(checklines[:,0])):
        if ((bbsupp == True) and (checklines[i,upper_col] < upper_en_cutoff) and (checklines[i,lower_col] < lower_en_cutoff) and (checklines[i,aval_col] > aval_min)):
            idx_to_remove[i,1] = 1
    #After this procedure, all lines w/ a valid radfield will have a '1' in their corresponding entry in idx_to_remove
    #We set by descending order of col0:
    idx_to_remove = idx_to_remove[idx_to_remove[:,0].argsort()[::-1]]    
    for i in range(0,len(idx_to_remove[:,0])):
        if (idx_to_remove[i,1] == 0):
            checklines = np.delete(checklines,idx_to_remove[i,0],0)
    #The checklines array with deleted rows is now what we want as our "input" line data for the model.
    flor[element_string][orbit_id]['relevant_lines'] = checklines
    return(flor)

def vac_to_air(wavelength_vac):
    #IAU Standard Conversion (Morton 2000)
    #NOTE: inputs assumed in nm, which we convert to angstroms required for the conversion:
    #NOTE: Previous versions of this function in my codes had 's' instead of 's**2'; please verify you are using correct version
    #if copy+pasting from my older codes.
    wavelength_vac = wavelength_vac * 10 #convert to ang
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength = wavelength_vac / n
    wavelength_conv = wavelength / 10 #convert back to nm
    return(wavelength_conv)

#%%
def error_calc(fluxes,raw_lines,raw_levs,ritz_col,uncert_col,lower_col,upper_col,aval_col,num_samples,model_lines,renorm=False):
    """
    This is a custom function for iterating the fluorescence model and deriving approximate uncertainties on line intensities.
    There are a large number of user inputs; these inputs include the original model inputs, and additional ones that the iteration procedure requires.
    
        USER INPUTS:
    uncert_col: column # (note 0 indexing in python) of the NIST A value uncertainty ratings.
    num_samples = number of model iterations to run
    model_lines: lines array (tuple element 1 from "fluorescence_spectra" output); will be changed
    in future release.
    """
    
    scatter = np.zeros((len(model_lines[:,0]),3+num_samples))
    scatter[:,0] = model_lines[:,0]
    scatter[:,1] = model_lines[:,1]
    scatter[:,2] = model_lines[:,2]
    original_avals = raw_lines[:,aval_col]
    uncerts = grab_uncertainty(raw_lines,uncert_col)
    for i in range(0,num_samples):
        new_avals = generate_new_avals(original_avals,uncerts)    
        new_lines = np.copy(raw_lines)
        #We save the new avals over the old ones in new_lines arr.
        #The new_avals has dimension (#lines,1), and we want (#lines). Reshape:
        new_lines[:,aval_col] = np.reshape(new_avals,len(new_avals))
        spectra = fluorescence_spectra(fluxes,new_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,renorm=False)  
        new_spec = spectra[1]
        new_intens = new_spec[:,2]
        scatter[:,3+i] = new_intens
    return(scatter) 

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
def numerically_integrate_2col(x_arr,y_arr):
    """
    Simple function for numerically integrating values using rectangles
    The inputs are x_arr, the x values, and y_arr, the y valuues.
    Note that the limits are not included; the limits must be taken care of by ONLY passing in values to be integrated.
    This code will integrate over ALL of the data input here.
    """
    integrated = 0
    for i in range(0,len(x_arr[:])-1):
        lower_x = x_arr[i]
        upper_x = x_arr[i+1]
        diff = upper_x - lower_x
        avg = (y_arr[i] + y_arr[i])/2
        integrated = integrated + (avg*diff)
    return(integrated)      

#%%
def matrix_solve(lhs,rhs,off_diag=False,diag=False):
    """
    Custom function for solving Ax=B in this context. 
    In order, this function tries:
            1. Regular numpy linear algebra inverse calculation.
            2. Psuedo-inverse
            3. Explicit SVD
        If these attempts fail, it will return a "solve_flag" = 1, which will trigger the code to look for and trim 'bad levels 1 by 1
    """
    solve_flag = 1
    solve_start_time = time.time()
    #Added on 03-02-2021 to better handle the fluorescence model.
    pops = np.zeros((len(rhs),1))
    #We need to do this iteratively. For most systems, the solution easily follows from np.linalg.inv.    
    try:
        a_inverse = np.linalg.inv(lhs)
        pops = np.dot(a_inverse,rhs)    
        if (zeros_check(pops) == False):
            solve_end = time.time()
            print('Matrix solution Completed in {:} seconds'.format(round(solve_end - solve_start_time,3)))
            solve_flag = 0
            neg_check = 0
            for i in range(0,len(pops[:])):
                if (pops[i] < 0):
                    if (neg_check == 0):
                        print('!!! WARNING !!! Negative population for level {:}'.format(i))
                        print('Multiplying ALL pops by -1..')
                    pops = pops * -1
                    neg_check = neg_check + 1
            if (neg_check > 1):  
                print('Negative populations observed again; solution likely invalid. Check Atomic Data inputs before re-running.')
        else:
            print('Zero populations detected in populations! Likely problem in missing or corrupt atomic data. Exiting code.')
            sys.exit()
    #In *some* unique cases, the matrix can be singular due to floating point rounding, or bad atomic data / large fluxes etc. 
    #We will try a couple different approaches:
    except:
        print('!!! WARNING !!! Singular or Near-Singular Matrix unsuitable for np.inv. Attempting alternative methods...')
        lhs_pert = np.copy(lhs)
        solution_check = False #Used to check for adequate solution:
        #We will try, in order, PINV, SVD.
        a_inverse = np.linalg.pinv(lhs_pert)
        pops = np.dot(a_inverse,rhs)
        #Check for negative populations:
        for i in range(0,len(pops[:])):
            if (pops[i] < 0):
                pops = pops * -1
                neg_check = neg_check + 1
        if (neg_check > 1):  
            print('Negative populations observed in Pseudo-Inverse. Attempting SVD Solution...')
        else:
            if (zeros_check(pops) == False):
                solve_flag = 0
                print('Valid solution found with Numpy Pseudo-Inverse')
        if (solve_flag > 0):
            U,S,V=np.linalg.svd(lhs_pert)
            a_inverse = np.dot(V.transpose(), np.dot(np.diag(S**-1),U.transpose()))
            pops = np.dot(a_inverse,rhs)
            neg_check = 0
            for i in range(0,len(pops[:])):
                if (pops[i] < 0):
                    pops = pops * -1
                    neg_check = neg_check + 1
            if (neg_check > 1):  
                print('Negative populations observed in SVD Solution. Solution is likely invalid.')
            else:
                if (zeros_check(pops) == False):
                    solve_flag = 0
                    print('Valid solution found with Singular Value Decomposition')
    return(solve_flag,pops)    

def zeros_check(arr):
    #Added on 03-02-2021; searches quickly for 0s in the input array. If present, break and return True in the bool.
    check_bool = False
    for p_check in range(0,len(arr[:])):
        if (arr[p_check] == 0):
           check_bool = True
           break
    return(check_bool)










