# -*- coding: utf-8 -*-
"""
Fluorescence Model using NIST data
SJB
Released tentatively on 01-04-2021
------------------------
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
# Some useful constants:
Rsun = 6.957e8 # meters
h = 6.62607004e-34
hbar = 1.054571817e-34
k = 1.38064852e-23
c = 299792458
cm_to_J = 1.986491855e-23

def doppler(wavelength,velocity):
    #positive velocity defined as away from sun
    #Updated 1-12-2021
    c = 299792458
    vel = velocity #assuming m/s
    rest_wavelength = wavelength #nm
    frac = 1 + (velocity/c)
    shifted_wavelength = rest_wavelength * frac
    return(shifted_wavelength)

def blackbody_wave(wavelength, temp=5777):
    #temp in Kelvin, all in SI
    #Adjusted on 01-13-2021 to clarify units
    #This function is not a typical blackbody; it reads in a wavelength (in nm), and outputs
    #an intensity in W/m^2 per nm; NOTE: We have taken care of the factor of ster in the 4pi/c multiplication
    #Assuming defined in nm, we convert:
    wavelength = wavelength * 1e-9
    h = 6.62607004e-34
    c = 299792458
    k = 1.38064852e-23
    exp_factor = h*c / (wavelength * k * temp) #compute item inside exponential factor
    exp_term = 1/(np.exp(exp_factor) - 1) #calculate exponential. Then, calculate BB intensity (W/m^2 per wavelength(m) per ster.):
    intensity = (2 * h * c**2 / (wavelength**5) )* exp_term
    intensity = intensity / 1e9
    return intensity

def find_nearest_index(array,value):
    diff_temp = abs(value-array[0])
    index = 0 
    for i in range(0,len(array[:])):
        diff = abs(array[i] - value)
        if (diff < diff_temp):
            index = i
            diff_temp = diff
    return index

def find_nearest_index_mod(array,value):
    #This is a modified version of find_nearest_index that cuts off searching through entire array once value is found
    #Requires that input array is sorted in increasing order.
    #In theory, this could be sped up by saving the index and passing that to the iteration for the next line
    #So that early values aren't re-searched. May implement in future if speed is concern:
    for i in range(0,len(array[:])):
        diff = array[i] - value
        if (diff > 0):
            lower_val = array[i-1]
            if (abs(lower_val) < abs(diff)):
                index = i-1
                break
            else:
                index = i
                break
    return index

def vac_to_air(wavelength_vac):
    #IAU Standard Conversion (Morton 2000)
    #NOTE: inputs assumed in nm, which we convert to angstroms required for the conversion:
    #NOTE: Previous versions of this function in my codes had 's' instead of 's**2'; please very you are using correct version!
    wavelength_vac = wavelength_vac * 10 #convert to ang
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength = wavelength_vac / n
    wavelength_conv = wavelength / 10 #convert back to nm
    return(wavelength_conv)
#%%
def convert_to_float(frac_str):
    #Updated on 03-03-21 (SJB) to remove question marks from NIST level lists after they were encountered and found to be problematic for Fe I
    #This converts a string to float; it is intended to do 2 things:
    # 1. Remove '?' characters from "questionable" levels in NIST level files, and 
    # 2. Convert half-integer J's, e.g. 3/2, to decimals (1.5)
    number_string = (frac_str.split('?'))[0]
    try:
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
def numerically_integrate(array):
    #Simple function for numerically integrating values; used rectangles:
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
    err_out[:,1] = errdat[:,1] #Vac Wave 
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
def fluorescence_spectra(fluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,renorm=True,geo_vel=0,minor_corr=True):  
    """
    Updated 03-02-2021 to trim 'bad' atomic data for small heliocentric distances
    by large, complex matrices with similar columns 
    This function calculates fluorescence spectra from NIST data. The function returns a tuple; tuple elements are:
    Element 0: level populations, indexed the same as "raw_levs"
    Element 1: line wavelengths and intensities normalized to maximum intensity.
                                                USER INPUTS:                           
    fluxes: Integrated flux in W/m^2; calculated with "grab_integrated_fluxes" function
    raw_lines: string array of NIST line data. Wavelengths must be in nm
    raw_levs: string array of NIST level data. Ensure both lines and levels have no headers. Energies must be in cm^-1
    ritz_col: column for ritz wavelengths (in VACUUM)
    lower_col: column of lower energy in lines array. Note 0 indexing, and energies must be in cm^-1
    upper_col: column of upper energy. Same limitations as 'lower_col'
    aval_col: column for A values in lines array.
    renorm: optional; assumed True. Normalizes line intensities such that the maximum intensity = 1.  
        -If scaling to REAL spectra to derive ABSOLUTE quantities e.g. column densities, set renorm=False
    geo_vel: geocentric velocity of the emitter in m/s. Default is 0.
    minor_corr: parameter for enabling minor perturbations added to the matrix solution in cases of singularity.
        -Used in upcoming version; default to True
    """
    #Define arrays for lines/levels after we process them:
    lines = np.zeros((len(raw_lines[:,0]),len(raw_lines[0,:])),dtype=float)
    levs = np.zeros((len(raw_levs[:,0]),len(raw_levs[0,:])),dtype=float)
    bad_levels_out = np.empty((0)) #Returns nothing for a good run; for a bad run, will contain the probable 'bad' levels found through model iterations.
    #Save data in usable form and convert half-integer J's to floats where necessary in the lines and levels arrays:
    print('Loading Atomic Data...')
    for i in range(0,len(lines[:,0])):
        for j in range(0,len(lines[0,:])):
            try:
                lines[i,j] = convert_to_float(raw_lines[i,j])
            except ValueError:
                pass
    for i in range(0,len(levs[:,0])):
        for j in range(0,len(levs[0,:])):
            try:
                levs[i,j] = convert_to_float(raw_levs[i,j])
            except ValueError:
                pass
    #Generate the internal indexing scheme:
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
    print('Starting Matrix Solution...')
    model_solution = matrix_solve(lhs,rhs)   #Attempt to solve matrix equation using, in order, matrix inversion, pseudo-inverse, and Singular value decomp.
    solution_flag = model_solution[0] 
    #If the solution flag is 0, we have a valid solution. If = 1, we have a bad solution and need to trim atomic data.
    if (solution_flag == 0):
        print('Solution is VALID and non-negative!')
        pops = model_solution[1]
    #If the model has run correctly, none of the below code will execute. However, if it hasn't run correctly, we start to look for bad atomic data and trim it:     
    else:
        #If we have gotten to this point, I am sorry. Good luck!
        pops = np.zeros(len(lhs[:,0])) #so the code can actually return info to user, we define an empty 'pops' variable
        while (solution_flag != 0):
            """
                    INTEGRITY CHECKS:
            """
            data_check = identify_bad_data(lhs)
            #Data check looks for problematic columns, i.e. identical, and we  then make them NOT identical.
            if (len(data_check) == 0):
                print('Bad Level detection failed. Revisit atomic data. Exiting code.')
                break
            bad_levels = np.sort(data_check)
            bad_level = np.amax(data_check)
            if ((len(bad_levels) == 1) and (bad_level==0)):
                print('Ground state found as only "bad" level. Please check atomic data for validity. Exiting code.')
                break
            else:
                if (minor_corr == True): #This is the default; could also do diag / non-diag perturbs.
                #If we have made it here, we have identified a column(s) that is not the 0th column that are identical.
                #We will slightly modify these columns w/ a near-zero perturbation to allow a solution to proceed.
                    print('Option for minor perturbative corrections not yet available. Exiting code.')
                    break
                
            """    
                if (off_diag == True):
                    print('Off-diagonal perturbative method not available yet. Exiting code.')
                    sys.exit(0)
                if (diag == True):
                    print('Diagonal perturbative method not available yet. Exiting code.')
                    sys.exit(0)
            if ((bad_level == 0) and (len(bad_levels) != 1)): #If true, the ground is incorrectly selected as the 'bad' and only 'bad' level. This was caused by the normalization condition.
                bad_level = bad_levels[-1] #Take highest level as the 'bad' level    
            print('(Approximate) Bad Level at original level index {:}. Removing from model.'.format(bad_level))
            #Save the 'bad' levels ORIGINAL index!
            #We look in "rel_levs" before we modify it to grab the original level index:
            for x in range(0,len(rel_levs[:,0])):
                if (rel_levs[x,2] == bad_level):
                    original_index = int(rel_levs[x,0]) #Must be int; float messes everything up.
                    break
            #Save the original index of the bad level for troubleshooting purposes after code end
            bad_levels_out = np.append(bad_levels_out, np.array([original_index]),axis=0)
            new_idx_and_rel_levs = trim_bad_atomic(idx,rel_levs,original_index)
            idx = new_idx_and_rel_levs[0]
            rel_levs = new_idx_and_rel_levs[1]
            #With the new line array and level mapping, we generate + populate the LHS.
            #Grab the "new_idx" from maximum of rel_levs array col2:
            #max_lev_idx = int(np.amax(rel_levs[:,2])) + 1 #To account for zero-indexing in python:
            new_idx = new_idx_and_rel_levs[2]
            lhs_empty = np.zeros((new_idx,new_idx)) #A * [X] = B; here lhs_empty will become matrix A
            lhs = populate_lhs(lhs_empty,idx,rel_levs) #Pass in lines, idx array (containing atomic data), and empty array to populate
            print('Starting Matrix Solution on modified matrix after removing level {:}.'.format(bad_level))
            rhs = np.zeros((len(lhs[:,0]),1))
            #Normalize to total pop = 1; enforce in first row:
            rhs[0] = 1
            for i in range(0,len(lhs[0,:])):
                lhs[0,i] = 1
            model_solution = matrix_solve(lhs,rhs,off_diag)   #Attempt to solve matrix equation
            solution_flag = model_solution[0]     
            pops = model_solution[1]
            if (solution_flag == 0):
                print('Solution Found!')
                break
            """
    #All said and done, we now calculate line intensities:
    line_intensities = calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm)
    #Done.
    return(pops,line_intensities,idx,lhs,rhs,rel_levs,bad_levels_out)

#%%
def fluxes_with_profiles(raw_rad,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist,t_comet, m_species,wav_coverage=0.005,supp=True,temp=5777, num_profile_points=500,profiles='Doppler'):
    """
    This function is used to integrate a radiation field convolved with a line profile for each line provided.
    The output is intended to be fed into the fluorescence_model function
    The inputs are:
        raw_rad: The 2-col solar spectrum used. Format: col0 wavelength (nm), col1 flux (W / m^2 / nm)
        raw_lines: lines array imported from a file downloaded for NIST; see documentation.
        raw_levs: levels array imported from a file downloaded for NIST; see documentation.
        lower_col: column index for lower level energy (cm^-1) in lines array
        upper_col: column index for upper level energy (cm^-1) in lines array
        aval_col: column index for einstein A values in lines array
        obj_vel: heliocentric velocity in m/s of object; positive is AWAY from sun.
        solar_dist: heliocentric distance of emitting object
        t_comet: temperature of emmitting species in Kelvin; 100K is approx. for comets.
        m_species: mass of emitting species in kg

    OPTIONAL inputs:
        wav_coverage: Defaulted to 0.005; this sets the range over which the radiation field is integrated.
            It is advised that this is slightly larger than the bulk of the line profile. Also be aware that this
            is expected to be span many pixels in the flux wavelength grid
        supp: defaulted to True; if the wavelength is outside the range of the radiation field "raw_rad", a blackbody
            at 5777K ('temp' default value) is used.
        num_profile_points: number of grid points to use for the line profiles. If changed, be sure that the # points
            is sufficient to resolve the structure of the profile. If set too low, the peak of the profile will be too broad
        profiles: choice of line profile. As of 01-28-2021, only 'Doppler' (i.e. a Doppler-broadened gaussian) is available.
            Future choices will be 'Lorentzian' or 'Voigt'
    """
    print('Starting Calculation of Integrated Fluxes...')
    start = time.time()
    cm_to_J = 1.986491855e-23
    Rsun = 6.957e8 # meters
    integrated_fluxes = np.zeros((len(raw_lines[:,0]),1),dtype=float)
    start_idx = 0 #This variable will be used as a 'tracker' of sorts when we search the raw_rad variable for the right x variable.
    #This is necessary as "raw_raw" is ~2 million lines (if using the providede Kurucz spectrum), and using this method DRAMATICALLY improves the speed!
    for i in range(0,len(raw_lines[:,0])):     
        upper_energy = float(raw_lines[i,upper_col])
        lower_energy = float(raw_lines[i,lower_col])
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))
        doppler_wave = doppler(wavelength_vac,obj_vel)
        #We need to define the lineshape for integration purposes. We will use only a gaussian at temperature "comet_temp" in K
        if (profiles == 'Doppler'):
            profile = doppler_dist(doppler_wave, t_comet=t_comet, m_species=m_species, num_points=num_profile_points, wav_coverage=wav_coverage)
        if (profiles == 'Lorentzian'):
            print('Option `Lorentzian` profile (natural linewidth) not yet available. Exiting code execution.')
            print('Use Option `Doppler`.')
            sys.exit(0)
        if (profiles == 'Voigt'):
            print('Option `Voigt` profile not yet available. Exiting code execution.')
            print('Use Option `Doppler`.')
            sys.exit(0)
        if (profiles == 'Custom'):
            print('Option `Custom` not yet available. Exiting code execution.')
            print('Use Option `Doppler`.')
            sys.exit(0)
        #Changed from rest wavelength to doppler wavelength on 01-27-2021
        if ((supp == True) and ((doppler_wave < raw_rad[0,0]) or (doppler_wave > raw_rad[-1,0]))):
            #Updated on 01-27-2021 to account for line profile
            #Generate the blackbody intensities on the same x grid as the doppler profile:
            dist_scale  = (Rsun/solar_dist)**2
            bb = blackbody_wave(profile[:,0])
            #Calculate the integrated product of the radiation field and profile:
            rad_at_wav = dist_scale * numerically_integrate_2col(profile[:,0], bb * profile[:,1])
            integrated_fluxes[i] = rad_at_wav
            print('Supplemental Blackbody invoked for line at {:} nm'.format(wavelength_vac))
        if ((supp == False) and ((doppler_wave < raw_rad[0,0]) or (doppler_wave > raw_rad[-1,0]))):
              #Only matters when 'supp=False' AND radiation field is not available for all relevant wavelengths
               print('!!! WARNING !!! No Radiation field for line {:} nm; all rates for line set to 0! Consider blackbody supplementation'.format(wavelength_vac))
               integrated_fluxes[i] = 0
        if (raw_rad[0,0] < doppler_wave < raw_rad[-1,0]):
                dist_scale_1au = (1.496e+11/solar_dist)**2
                #We want the radiation field on the same grid as the line profile, i.e. a fine grid.
                #Therefore, we first search the radiation field for the START of the line profile, i.e. the min x in 'profile' variable.
                #We use the improved function 'fniiv' (Find Nearest Index Iteration Version), and pass in the "start_idx" :)
                searched_rad_indices = fniiv(raw_rad[:,0],profile[0,0],start_idx)
                start_idx = searched_rad_indices[0]
                rad_start_idx = searched_rad_indices[0] - 10 #Save the index of the smallest wavelength from profile array
                start_idx = rad_start_idx #Save the suggested index for the next iteration, i.e. the search for the next wavelength.
                rad_center_idx = fniiv(raw_rad[:,0],doppler_wave,start_idx)[0] 
                rad_end_idx = fniiv(raw_rad[:,0], profile[-1,0],start_idx)[0] + 10 #find end of profile wavelength range in radiation array w/ 1 pixel buffer
                #Pass this information on and use interpolation to find the radiation field values over the wavelength grid of the line profile:        
                new_rad_grid = dist_scale_1au * generate_rad_grid(raw_rad[rad_start_idx:rad_end_idx,:],profile[:,0])
                #We now integrate over the line profile using rectangles:
                rad_at_wav = numerically_integrate_2col(new_rad_grid[:,0], new_rad_grid[:,1] * profile[:,1])
                #Save the fluxes belonging to line 'i' :
                integrated_fluxes[i] = rad_at_wav
    end = time.time()
    print('Calculation of Integrated Fluxes finished. Elapsed Time {:} s'.format(round(end - start,2)))
    return(integrated_fluxes)   

#%%
def generate_rad_grid(rad, wavelengths):
    """
    Updated documentation on 03-03-2021
    Introduced on 01-27-2021 to transform the radiation values in rad array to those of the fine-grid line profile. 
    This function takes in a radiation field in 2 col format; col0 wavelength, col1 flux., and a wavelength grid "wavelengths"
    It is presumed that the grid is MUCH finer for wavelengths compared to rad.
    This function takes the values in 'grid' and places them on the same grid spacing as "wavelengths" using interpolation.     
    """   
    output_arr = np.empty((len(wavelengths),2))
    output_arr[:,0] = wavelengths
    lower_radx = 0
    upper_radx = 0
    for i in range(0,len(wavelengths[:])):
        #The # of points in rad is very small, so we'll just brute force the search for the appropriate indices when they need updating:
        if ((rad[lower_radx,0] > wavelengths[i]) or (rad[upper_radx,0] < wavelengths[i])):
            for j in range(0,len(rad[:,0])):
                check1 = find_nearest_index(rad[:,0], wavelengths[i])
                #really we just need the lower wavelength found, then the upper index is just lower + 1
                if (rad[check1,0] < wavelengths[i]):
                    lower_radx = check1
                    upper_radx = check1 + 1
                if (rad[check1,0] > wavelengths[i]):
                    upper_radx = check1
                    lower_radx = check1 - 1
        #Interpolate and save the new y value:
        x_val = wavelengths[i]
        lower_x = rad[lower_radx,0]
        upper_x = rad[upper_radx,0]
        lower_y = rad[lower_radx,1]
        upper_y = rad[upper_radx,1]
        slope = (upper_y - lower_y) / (upper_x - lower_x)
        interp_y = slope*(wavelengths[i] - lower_x) + lower_y
        output_arr[i,1] = interp_y
    return(output_arr)
#%%
def fniiv(array,value,start_idx):   
    """
    "Find Nearest Index Iterative Version"
    Added on 01-27-2021
    This is essentially a "find nearest index" code in a sorted array,
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
def numerically_integrate_2col(x_arr,y_arr):
    """
    SJB Documentation updated 03-03-2021
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
def doppler_dist(wav_line, t_comet, m_species, num_points, wav_coverage=0.005):
    """
    02-05-2021: Updated to properly normalize line profile. 
    This function returns a doppler broadened line profile at temperature "t_comet" (K), mass "m_species" (kg) and wavelength "wav_line" (nm). 
    This function is called by the "fluxes with profiles" function; it takes in a (doppler-shifted) wavelength and returns a line profile.
    "wav_coverage" is the range of the profile, e.g +/- wav_coverage in nm units.
    num_points is the number of points in the wavelength grid returned. 
    The output is in 2 column format:
    col0: wavelength, col1: profile (set to 1 at the center wavelength)
    """
    c = 2.99792458e8  #SI
    c_nm = c * 1e9
    kb = 1.3864852e-23 #SI
    wavelength_nm = wav_line
    wavelength_m = wav_line * 1e-9 #convert o m
    freq_rest = c / wavelength_m #calcute frequency
    gauss_std = np.sqrt(kb * t_comet / (m_species)) * (freq_rest / c) #standard deviation of Gaussian
    #The wavelength grid in units of nm
    x_range = np.linspace(wav_line - wav_coverage, wav_line + wav_coverage, num_points)
    #convert to FREQUENCY for ease of calculation:
    gauss_x = c_nm * (wavelength_nm - x_range) / (wavelength_nm * x_range)
    #Calculate HWHM from standard deviation:
    alpha = gauss_std * np.sqrt(2 * np.log(2))  
    #Profile:
    """
    gauss_std = np.sqrt(kb * temp / (mass)) * (freq_rest / c)
    x_range = np.linspace(wavelength_nm - 0.005, wavelength_nm + 0.005, 100)
    gauss_x = c / (x_range - wavelength_nm) * 1e-9      #FREQUENCY SPACE
    lor_x = c_nm * (wavelength_nm - x_range) / (wavelength_nm * x_range)
    test_gauss = G(lor_x,gauss_std)       
    """
    profile = np.sqrt(np.log(2) / np.pi)/alpha * np.exp(- (gauss_x / alpha)**2 * np.log(2))
    output = np.empty((num_points,2))
    #Save the WAVELENGTH grid and line profile in 2 col format:
    output[:,0] = x_range
    output[:,1] = profile    
    #Normalize profile so intregal of profile over all wavelengths = 1 !
    gauss_norm = numerically_integrate_2col(x_range,profile[:]) #Normalization factor check!
    output[:,1] = output[:,1] / (gauss_norm)
    return(output)

def lorenztian_dist():
    #Code option available in future release.
    return()

def voigt_dist():
    #Code option available in future release    
    return()

#%%
#Some experimental perturbative functions.    

def diag_perturb(arr,perturb=1e-10):
    #Added on 03-01-2021
    #Assumes a square matrix. Adds a perturbation (default 1e-10) to diagonal elements.
    for r in range(0,len(arr[:,0])):
        arr[r,r-1] = arr[r,r-1] + perturb
    return(arr)  

def off_diag_perturb(arr,perturb=1e-10):
    #Added on 03-02-2021
    #Assumes a square matrix. Adds a perturbation to OFF diagonal elements!
    for r in range(1,len(arr[:,0])):
        arr[r,r] = arr[r,r] + perturb
    return(arr)  
      
#%%
def zeros_check(arr):
    #Added on 03-02-2021; searches quickly for 0s in the input array. If present, break and return True in the bool.
    check_bool = False
    for p_check in range(0,len(arr[:])):
        if (arr[p_check] == 0):
           check_bool = True
           break
    return(check_bool)
#%%
def matrix_solve(lhs,rhs,off_diag=False,diag=False):
    """
    Added on 03-02-2021 to functionalize matrix solving and return an indicator that the solution didn't work
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
                        print('Multiplying ALL pops by -1 (scalar multiple of vector is still a vector).')
                    pops = pops * -1
                    neg_check = neg_check + 1
            if (neg_check > 1):  
                print('Negative populations observed again. Rate Matrix returned in tuple element 3. Check Atomic Data inputs before re-running.')
        else:
            print('Zero populations detected in populations! Likely problem in missing or corrupt atomic data. Exiting code.')
    #In *some* unique cases, the matrix can be singular due to floating point rounding, or bad atomic data / large fluxes etc. 
    #We will try a couple different approaches:
    except:
        print('!!! WARNING !!! Singular or Near-Singular Matrix Detected. Using alternative methods...')
        lhs_pert = np.copy(lhs)
        solution_check = False #Used to check for adequate solution:
        #We will try, in order, PINV, SVD.
        a_inverse = np.linalg.pinv(lhs_pert)
        pops = np.dot(a_inverse,rhs)
        if (zeros_check(pops) == False):
            sol_string = 'Numpy Psuedo-inverse'
            solve_flag = 0
        else:
            U,S,V=np.linalg.svd(lhs_pert)
            a_inverse = np.dot(V.transpose(), np.dot(np.diag(S**-1),U.transpose()))
            pops = np.dot(a_inverse,rhs)
            if (zeros_check(pops) == False):
                sol_string = 'Singular Value Decomposition'
                solve_flag = 0 
        """
        perturb_iteration = 0 #Iteration counter for perturbation. 
        if (diag==True):
            diagonal_perturb = 1e-6 #starting value of diagonal perturb
        else:
            diagonal_perturb = 0
        #For each iteration, we try a diagonal perturbation and attempt to solve using pseudo-inverse
        lhs_pert = np.copy(lhs)
        if (diag == True):
            max_num_iter = 50
        else:
            max_num_iter = 0
        while (perturb_iteration < max_num_iter):            
            lhs_pert = diag_perturb(lhs_pert,diagonal_perturb)
            a_inverse = np.linalg.pinv(lhs_pert)
            pops = np.dot(a_inverse,rhs)
            perturb_iteration = perturb_iteration + 1 #Add to iteration variable.
            #Now we check that the solution was valid.
            #zeros_check is a function that returns true if ANY value of passed in array is = to 0.
            if (zeros_check(pops) == False):
                sol_string = 'Numpy Psuedo-inverse'
                solution_check = True
            else: #If zero check is True, pinv solution is no good.
                #Try an SVD Solution:
                U,S,V=np.linalg.svd(lhs_pert)
                a_inverse = np.dot(V.transpose(), np.dot(np.diag(S**-1),U.transpose()))
                pops = np.dot(a_inverse,rhs)
                if (zeros_check(pops) == False):
                    sol_string = 'Singular Value Decomposition'
                    solution_check = True
                else:
                    a_inverse = ip.invert(lhs_pert,rhs,100,1)
                    #ip.invert returns a tuple; we want the element [1] 
                    pops = np.dot(a_inverse[1],rhs)
                    if (zeros_check(pops) == False):
                        solution_check = True
                        sol_string = 'Tikhonov Regularization'
                    else:
                        diagonal_perturb = diagonal_perturb *1.1          
            if (solution_check == True):                   
                iteration_end = time.time()
                solve_flag = 0
                print('Matrix solution Completed in {:} seconds after {:} iterations using {:} with a perturbation of \
                      +{:}'.format(round(iteration_end - solve_start_time),perturb_iteration,sol_string,diagonal_perturb))
                break
            if (perturb_iteration == max_num_iter + 1):
                print('No solution found using iterative *diagonal* perturbation or Tikhonov Regularization. All pops set to zero.')
        if ((solution_check == False) and (off_diag == True)):
            print('Exhausted Diagonal Perturbation Options. Attempting minor *off-diagonal* correction.')
            perturb_iteration = 0 #Iteration counter for perturbation. 
            diagonal_perturb = 1e-3 #starting value of diagonal perturb
            #For each iteration, we try an OFF diagonal perturbation and attempt to solve using pseudo-inverse
            lhs_pert = np.copy(lhs)
            max_num_iter = 50
            while (perturb_iteration < max_num_iter):            
                lhs_pert = off_diag_perturb(lhs_pert,diagonal_perturb)
                a_inverse = np.linalg.pinv(lhs_pert)
                pops = np.dot(a_inverse,rhs)
                perturb_iteration = perturb_iteration + 1 #Add to iteration variable.
                #Now we check that the solution was valid.
                #zeros_check is a function that returns true if ANY value of passed in array is = to 0.
                if (zeros_check(pops) == False):
                    sol_string = 'Numpy Psuedo-inverse'
                    solution_check = True
                else: #If zero check is True, pinv solution is no good.
                    #Try an SVD Solution:
                    U,S,V=np.linalg.svd(lhs_pert)
                    a_inverse = np.dot(V.transpose(), np.dot(np.diag(S**-1),U.transpose()))
                    pops = np.dot(a_inverse,rhs)
                    if (zeros_check(pops) == False):
                        sol_string = 'Singular Value Decomposition'
                        solution_check = True
                    else:
                        a_inverse = ip.invert(lhs_pert,rhs,100,1)
                        #ip.invert returns a tuple; we want the element [1] 
                        pops = np.dot(a_inverse[1],rhs)
                        if (zeros_check(pops) == False):
                            solution_check = True
                            sol_string = 'Tikhonov Regularization'
                        else:
                            diagonal_perturb = diagonal_perturb *1.1          
                if (solution_check == True):                   
                    iteration_end = time.time()
                    solve_flag = 0
                    print('Matrix solution Completed in {:} seconds after {:} iterations using {:} with a perturbation of \
                          +{:}'.format(round(iteration_end - solve_start_time),perturb_iteration,sol_string,diagonal_perturb))
                    break
                if (perturb_iteration == max_num_iter + 1):
                    print('No solution found using iterative off-diagonal perturbation or Tikhonov Regularization. All pops set to zero.')            
        #EXPERIMENTAL
        #EXPERIMENTAL
                """
    return(solve_flag,pops)    
#%%
def generate_index_scheme(lines,levs,fluxes,upper_col,lower_col,aval_col,ritz_col):
    """
    Added 03-02-2021
    SJB
    Used internally by fluorescene_model to generate matrix "idx" used  to encapsulate all atomic data that goes into the LHS
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
        wavelength_air = lines[i,ritz_col]
        #freq = (upper_energy - lower_energy)*cm_to_J/h
        #Degeneracies and Avalue:
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))
        g_i = 2 *float(levs[int(idx[i,0]),2]) + 1 #lower level degeneracy
        g_j = 2 *float(levs[int(idx[i,1]),1]) + 1 #upper level degeneracy
        aval = float(lines[i,aval_col]) #spont. emission rate
        #Old stim_coeff = aval * np.pi**2 * c**3 / (hbar * omega**3)
        #New, added 01-29-2021 to switch to wavelength units:
        stim_coeff = (wavelength_vac * 1e-9)**3 * aval * 1/(2 * h * c)
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
#%%
def identify_bad_data(lhs_arr):    
    """
    Added on 03-02-2021 by SJB
    Codes on 03-02-2021 Brought to you by 8 hours of King Crimson's "21st Century Schizoid Man" on repeat
    This function is intended to identify potentially 'bad' atomic data based on the number of similarities in the LHS matrix
    If this code is called, you are in trouble; there is a big problem, and that requires big resources to solve/figure out.
    
    This code does not work well.
    
    """
    #Default for output; returns empty if no 'bad' data found.
    bad_levels = np.empty((0,1))
    bad_cols = []
    bad_rows = []
    # If the previous approaches failed, we need to identify the problem columns/rows. Now things get fun.
    lhs_rank = np.linalg.matrix_rank(lhs_arr)
    print('No solution possible. LHS dimension is {:}, rank is {:}. Searching for bad data (Slow)!'.format(len(lhs_arr[:,0]),lhs_rank))
    #First, let's check if any columns or rows are all similar. Let's do this iteratively, and remove the largest level index first
    #Since those are most likely to have the least # of lines; we'd prefer not to remove low (e.g. metastable) levels from the 
    #rate matrix. 
    #New method added on 03-02-2021:
    check0 = False
    #First, we check if any columns are identical using a quick numpy routine:
    cols_equal = np.empty((0,2))
    for i in range(0,len(lhs_arr[0,:])):
        for j in range(i+1,len(lhs_arr[0,:])):
            if (np.array_equal(lhs_arr[:,i],lhs_arr[:,j]) == True):
                temp = np.empty((1,2))
                temp[0,0] = i
                temp[0,1] = j
                print('Potential Problem: Columns {:} and {:} identical!'.format(i,j))
                cols_equal = np.append(cols_equal,temp,axis=0)
    if (len(cols_equal) > 0): #If True, we have identical columns
        check0 = True #flag to use the later method using histograms; if True, we ignore that for now.
        #Flatten the list and grab unique level #'s. These are the 'bad' columns:
        cols_flat = np.append(cols_equal[:,0], cols_equal[:,1], axis=0)
        bad_cols = np.unique(cols_flat)
        bad_levels = np.copy(bad_cols)
    """
    rows_equal = np.empty((0,2))
    for i in range(0,len(lhs_arr[:,0])):
        for j in range(i+1,len(lhs_arr[:,0])):
            if (np.array_equal(lhs_arr[i,:],lhs_arr[j,:]) == True):
                temp = np.empty((1,2))
                temp[0,0] = i
                temp[0,1] = j
                print('Potential Problem: Rows {:} and {:} identical!'.format(i,j))
                rows_equal = np.append(rows_equal,temp,axis=0)
    if (len(rows_equal) > 0): #If True, we have identical columns
        check0 = True #flag to use the later method using histograms; if True, we ignore that for now.
        #Flatten the list and grab unique level #'s. These are the 'bad' columns:
        rows_flat = np.append(rows_equal[:,0], rows_equal[:,1], axis=0)
        bad_rows = np.unique(rows_flat)    
    
    bad_levels = np.append(bad_cols,bad_rows)
    """    
    if (check0 == False):
        print('No identical columns found. Attempting to identical *similar* columns...')
        compare_thresh = 1e2
        comp_out = np.empty((0,1))
        #First, we compare  column by column:
        for i in range(0,len(lhs_arr[0,:])):
            init_col = lhs_arr[:,i]
            for j in range(i+1,len(lhs_arr[0,:])):
                temp = np.empty((1,1))
                temp2 = np.empty((1,1))
                comp_col = lhs_arr[:,j]
                col_diff = init_col - comp_col
                col_check = 0
                for k in range(1,len(col_diff[:])):
                    if (0 < abs(col_diff[k]) < compare_thresh):
                        col_check = col_check + 1
                if (col_check > 0):
                    temp[0,0] = i
                    temp2[0,0] = j
                    comp_out = np.append(comp_out,temp,axis=0)
                    comp_out = np.append(comp_out,temp2,axis=0)
        #We compare rows:   
        comp_out = np.empty((0,1))
        for i in range(0,len(lhs_arr[:,0])):
            init_row = lhs_arr[i,:]
            for j in range(i+1,len(lhs_arr[:,0])):
                temp = np.empty((1,1))
                temp2 = np.empty((1,1))
                comp_row = lhs_arr[j,:]
                row_diff = init_row - comp_row
                row_check = 0
                for k in range(1,len(row_diff[:])):
                    if (0 < abs(row_diff[k]) < compare_thresh):
                        row_check = row_check + 1
                if (row_check > 0):
                    temp[0,0] = i
                    temp2[0,0] = j
                    comp_out = np.append(comp_out,temp,axis=0)
                    comp_out = np.append(comp_out,temp2,axis=0)               
                    
        #The col_comp_out array contains the indices of levels that give us trouble.
        #We'll search for the largest numbers of occurences, and return those for later processing.
        hist_occur = (np.histogram(comp_out,bins=len(lhs_arr[:,0])))[0]
        #This calculates a histogram of the occurences of the level indices that had columns within the "compare_thresh"
        hist_occur_copy = np.copy(hist_occur)
        #bad_levels = np.sort(hist_occur_copy)
        
        max_of_hist = np.amax(hist_occur)
        #We search within +/- 10% of the maximum to be safe:
        hist_search_min = max_of_hist*.5
        hist_search_max = max_of_hist*1.1
        for i in range(0,len(hist_occur[:])):
             if (hist_search_min < hist_occur[i] < hist_search_max):
                 temp = np.empty((1,1))
                 temp[0,0] = i
                 bad_levels = np.append(bad_levels,temp,axis=0)
    return(bad_levels)
#%%
def trim_bad_atomic(idx_arr,rel_levs,bad_level_original_index):
    """
    03-02-2021
    Added by SJB
    This function is specialized, and used internally by the model function to remove 'bad' atomic data. This function does NOT 'find' bad atomic data,
    it only contains the functionality to remove select levels from the model. In the future, I may code this up so that the user can selectively remove levels
    if they wish.
    
    This code would handle the 'removal' from the model.
    """
    new_idx_arr = np.copy(idx_arr)
    new_rel_levs = np.copy(rel_levs)
    #Since we are going from len(arr) to >= 0, we subtract 1 to correct for 0 indexing:
    i = len(new_idx_arr[:,0]) - 1
    while (i >= 0):
        if ((new_idx_arr[i,0] == bad_level_original_index) or (new_idx_arr[i,1] == bad_level_original_index)):
            #If true, the bad level is involved in emission line 'i'. We remove it:
            new_idx_arr = np.delete(new_idx_arr, i, axis=0)
        i = i - 1
    #With that done, we turn our attention to the 'rel_levs' array. 
    #First, we identify where the 'bad level' is located. Thankfully, the rel_levs array spans the entirety of the original levels file.
    #In rel_levs, the columns are: col0 (original index), col1 (=1 if present in model), col2 = new index:
    #We essentially just set the col1 of the bad level to 0, and re-generate the new index scheme:
    #Note, we are assuming the BAD level is indexed according to the ORIGINAL indexing scheme, i.e. row # = bad level # in original index scheme
    new_rel_levs[bad_level_original_index,1] = 0
    new_idx = 0
    for i in range(0,len(new_rel_levs[:,0])):
        if (new_rel_levs[i,1] == 1):
            new_rel_levs[i,2] = new_idx
            new_idx = new_idx + 1
    return(new_idx_arr,new_rel_levs,new_idx)
#%%
def calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm):
    """
    Added on 03-02-2021
    SJB
    Functionalized version of line intensity calculation.
    Also removed references to "lines" and "levels" arrray; everything is now based on the internal 'idx' array containing
    all the atomic rate data.
    ------------------------------------------------------
    Now, we need to calculate intensities. The intensity of a line is defined as:
       I = N hbar omega Aval
       since hbar*omega is energy, we'll just calculate the energy (in Joules) from the vacuum wavelength
       from E = hc/lamba  
    """
    intens_start_time = time.time()
    print('Starting Intensity Calculation...')
    pops = np.real(pops) #Just in case; I've seen some cases early on that yielding complex numbers; no idea if this is still relevant.
    intens = np.zeros((len(idx[:,0]),3))   
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) intensity (arb. units.)
    for i in range(0,len(intens[:,0])):
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        upper_pop = pops[mapped_index] #Grab the (normalized) population
        energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
        aval = idx[i,4]
        line_intens = upper_pop * energy * aval #calculate line intensity from above
        intens[i,0] = doppler(idx[i,2],geo_vel) #Save (air) wavelength
        intens[i,1] = idx[i,3] #Save vac wavelength
        intens[i,2] = line_intens #save intensity
    #Normalize to max line = 1
    if (renorm == True):
        max_intens = np.nanmax(intens[:,2])
        intens[:,2] = intens[:,2] / max_intens
    intens_end_time = time.time()
    print('Intensity Calculations finished in {:} seconds'.format(round(intens_end_time - intens_start_time,3)))
    return(intens)
