# -*- coding: utf-8 -*-
"""
Fluorescence Model using NIST data
SJB
Released tentatively on 06-07-2021
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
    #Updated 05-24-2021
    #Wavelength in nm; velocity in m/s
    #Positive defined as away from object i.e. redshift
    c = 299792458
    vel = velocity #assuming m/s
    rest_wavelength = wavelength #nm
    frac = 1/(1 - velocity/c)
    shifted_wavelength = rest_wavelength * frac
    return(shifted_wavelength)

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

def find_nearest_index(array,value):
    diff_temp = abs(value-array[0])
    index = 0 
    for i in range(0,len(array[:])):
        diff = abs(array[i] - value)
        if (diff < diff_temp):
            index = i
            diff_temp = diff
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
    else:
        print('Problem with solution (levels with 0 population). Be wary of fluorescence efficiencies.')
        pops = model_solution[1]
    #If the model has run correctly, none of the below code will execute. However, if it hasn't run correctly, we start to look for bad atomic data and trim it:     
    
    #All said and done, we now calculate line intensities:
    line_intensities = calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm)
    #Done.
    return(pops,line_intensities,idx,lhs,rhs,rel_levs,bad_levels_out)

#%%
def fluxes_with_profiles(raw_rad_vac,raw_rad_air,raw_rad_ir,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist,supp=True,bbtemp=5777):
    """
    This function is used to integrate a radiation field convolved with a line profile for each line provided.
    The output is intended to be fed into the fluorescence_model function
    The inputs are:
        rad_uv: 2col solar fluxes; Col0 in VACUUM nm, col1 flux in W / m^2 / nm in vacuum wavelengths.
        rad_vis: 2col solar fluxes; Col0 in AIR nm, col1 flux in W / m^2 / nm in vacuum wavelengths.
        rad_ir: 2col solar fluxes; Col0 in VACUUM nm, col1 flux in W / m^2 / nm in vacuum wavelengths.
        
        For radiation fields, conversion to W/m^3 is handled internally. 
        
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

        supp: defaulted to True; if the wavelength is outside the range of the radiation field "raw_rad", a blackbody
            at 5777K ('temp' default value) is used.
    """
    print('Starting Calculation of Integrated Fluxes...')
    start = time.time()
    cm_to_J = 1.986491855e-23
    Rsun = 6.957e8 # meters
    integrated_fluxes = np.zeros((len(raw_lines[:,0]),1),dtype=float)
    start_idx_air = 0 #This variable will be used as a 'tracker' of sorts when we search the raw_rad variable for the right x variable.
    start_idx_vac = 0# Variable used as a similar 'tracker' for vacuum
    start_idx_ir = 0# Variable as a tracker for vacuum infrared
    for i in range(0,len(raw_lines[:,0])):     
        #Grab upper/lower energies and vac wavelength
        upper_energy = convert_to_float(raw_lines[i,upper_col])
        lower_energy = convert_to_float(raw_lines[i,lower_col])
        #Vac wavelength in meters:
        wavelength_vac_m = h * c / (cm_to_J * (upper_energy - lower_energy))
        wavelength_vac_nm = 1e9 * wavelength_vac_m
        wavelength_air_nm = vac_to_air(wavelength_vac_nm)        
        #Doppler-shifted wavelengths:
        doppler_wave_air = doppler(wavelength_air_nm,obj_vel)
        doppler_wave_vac = doppler(wavelength_vac_nm,obj_vel)
        #Three cases:
        #1: in air:
        if ((raw_rad_air[0,0] < doppler_wave_air < raw_rad_air[-1,0])):
            dist_scale_1au = (1.496e+11/solar_dist)**2
            searched_rad_indices = fniiv(raw_rad_air[:,0],doppler_wave_air,start_idx_air)
            #start_idx is passed to next iteration.
            start_idx_air = searched_rad_indices[1]
            #Index for the lines radiation field:
            rad_idx = searched_rad_indices[0]
            rad_at_wav = raw_rad_air[rad_idx,1] * 1e9 #converts from W/m^2 per nm to W /m^2 per meter (wavelength)
            integrated_fluxes[i] = rad_at_wav
        #2: in Vac:
        if ((raw_rad_vac[0,0] < doppler_wave_vac < raw_rad_vac[-1,0])):
            dist_scale_1au = (1.496e+11/solar_dist)**2
            searched_rad_indices = fniiv(raw_rad_vac[:,0],doppler_wave_vac,start_idx_vac)
            #start_idx is passed to next iteration.
            start_idx_vac = searched_rad_indices[1]
            #Index for the lines radiation field:
            rad_idx = searched_rad_indices[0]
            rad_at_wav = raw_rad_vac[rad_idx,1] * 1e9 #converts from W/m^2 per nm to W /m^2 per meter (wavelength)
            integrated_fluxes[i] = rad_at_wav
        #3: Vacuum infrared:
        if ((raw_rad_ir[0,0] < doppler_wave_vac < raw_rad_ir[-1,0])):
            dist_scale_1au = (1.496e+11/solar_dist)**2
            searched_rad_indices = fniiv(raw_rad_ir[:,0],doppler_wave_vac,start_idx_ir)
            #start_idx is passed to next iteration.
            start_idx_ir = searched_rad_indices[1]
            #Index for the lines radiation field:
            rad_idx = searched_rad_indices[0]
            rad_at_wav = raw_rad_ir[rad_idx,1] * 1e9 #converts from W/m^2 per nm to W /m^2 per meter (wavelength)
            integrated_fluxes[i] = rad_at_wav    
        #4: Outside all bounds; if supp=True, use blackbody:
        if ((supp == True) and ((doppler_wave_vac < raw_rad_vac[0,0])  or (doppler_wave_vac > raw_rad_ir[-1,0]))):
            #Only matters when 'supp=False' AND radiation field is not available for all relevant wavelengths
            print('!!! WARNING !!! No Radiation field for line {:} nm; all rates for line set to 0! Consider blackbody supplementation'.format(wavelength_vac_nm))
            integrated_fluxes[i] = 0
            dist_scale  = (Rsun/solar_dist)**2
            bb = blackbody_wave(doppler_wave_vac,temp=bbtemp) #Only need that specific wavelength at temp=5777K (default); assumes nm in
            rad_at_wav = dist_scale * bb
            integrated_fluxes[i] = rad_at_wav
            print('Supplemental Blackbody invoked for line at {:} nm'.format(wavelength_vac_nm))
        if ((supp == False) and ((doppler_wave_vac < raw_rad_vac[0,0])  or (doppler_wave_vac > raw_rad_ir[-1,0]))):
            integrated_fluxes[i] = 0
    end = time.time()
    print('Calculation of Integrated Fluxes finished. Elapsed Time {:} s'.format(round(end - start,2)))
    return(integrated_fluxes)   
#%%
def fniiv(array,value,start_idx):   
    """
    " Nearest Index Iterative Version"
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
        #
        stim_coeff = (wavelength_vac * 1e-9)**5 * aval * 1/(8 * np.pi * h * c**2)
        #
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
def calc_line_intensities(pops,idx,rel_levs,geo_vel,renorm):
    """
    Added on 03-02-2021
    SJB
    Functionalized version of line intensity calculation.
    Also removed references to "lines" and "levels" arrray; everything is now based on the internal 'idx' array containing
    all the atomic rate data.
    ------------------------------------------------------
    Now, we need to calculate intensities. The intensity of a line is defined as:
       I = N * hc/lambda *  Aval
       since hc/lambda is energy, we'll just calculate the energy (in Joules) from the vacuum wavelength
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
