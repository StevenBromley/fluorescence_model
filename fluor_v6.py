# -*- coding: utf-8 -*-
"""
Fluorescence Model using NIST data
SJB
Released tentatively on 01-04-2021
------------------------
"""
import numpy as np #The OG
import matplotlib.pyplot as plt
import math
from random import random #Used for generating new A values according to NIST uncertainty rating
import time
from scipy.integrate import quad #for integrating radiation field
import scipy
import os

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
    print(wavelength_vac)
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength = wavelength_vac / n
    wavelength_conv = wavelength / 10 #convert back to nm
    return(wavelength_conv)

def convert_to_float(frac_str):
    #This attempts to convert a string to float; made to handle J values that are half-integers:
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac

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
    #Added 01-08-2021
    #Used for converting A value uncertainty ratings in NIST to numbers
    #pass in lines array, and 'col' is the col# of the accuracy ratings
    a_uncert = np.zeros((len(lines[:,0]),1),dtype=float)
    for i in range(0,len(lines[:,0])):
        temp = np.zeros((1,1))
        letter = lines[i,4]
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
    #Added on 01-08-2021
    #Pass in str array of a values, and a 1d array of uncertainties (floats)
    #Then, generate a new A value within the range Aval +/- uncert.
    adjusted_avals = np.empty((len(avals[:]),1),dtype='U256')
    for i in range(0,len(avals[:])):
        uncert = uncerts[i]
        aval = float(avals[i])
        low_aval = aval - aval*uncert*0.01 #change from percent to decimal
        high_aval = aval + aval*uncert*0.01 #change from percent to decimal
        adjusted_avals[i] = (str(np.random.uniform(low_aval,high_aval)).replace('[','')).replace(']','')
    return(adjusted_avals)

#%%
def grab_integrated_fluxes(raw_rad,raw_lines,lower_col,upper_col,obj_vel,solar_dist,binsize=2,binsize_nm=0.002,supp=True,temp=5777):
    #This functions read in a radiation field, velocity, and distance, and integrates the rad field around each line
    #The default binsize (pixels) for input radiation field is 2 pixels, corresponding to default binsize_nm = 0.002nm
    #Default is to supplement the provided radiation field w/ a blackbody IF 'supp=True' (default) and radiation field
    #data is not provided at a wavelength; 
    print('Starting Calculation of Integrated Fluxes...')
    start = time.time()
    cm_to_J = 1.986491855e-23
    Rsun = 6.957e8 # meters
    integrated_fluxes = np.zeros((len(raw_lines[:,0]),1),dtype=float)
    for i in range(0,len(raw_lines[:,0])):
        upper_energy = float(raw_lines[i,upper_col])
        lower_energy = float(raw_lines[i,lower_col])
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))
        if ((supp == True) and ((wavelength_vac < raw_rad[0,0]) or (wavelength_vac > raw_rad[-1,0]))):
          rad_at_wav = ((Rsun/solar_dist)**2)*(quad(blackbody_wave,doppler((wavelength_vac - binsize_nm),obj_vel),doppler((wavelength_vac + binsize_nm),obj_vel))[0])
          integrated_fluxes[i] = rad_at_wav
        else:
            if (wavelength_vac > raw_rad[0,0]):
                rad_idx = find_nearest_index(raw_rad[:,0],doppler(wavelength_vac,obj_vel))
                #New; integrate from wav_vac - 0.02nm to wav_vac + 0.02nm
                rad_at_wav = numerically_integrate(raw_rad[rad_idx-binsize:rad_idx+binsize,:])
                integrated_fluxes[i] = rad_at_wav
            if (wavelength_vac < raw_rad[0,0]):
                #Only matters when 'supp=False' AND radiation field is not available for all relevant wavelengths
                integrated_fluxes[i] = 0
    end = time.time()
    print('Calculation of Integrated Fluxes finished. Elapsed Time {:} s'.format(end - start))
    return(integrated_fluxes)
#%%
def fluorescence_spectra(fluxes,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp=5777,supp=True,renorm=True):  
    """
    This function calculates fluorescence spectra from NIST data. The function returns a tuple; tuple elements are:
    Element 0: level populations, indexed the same as "raw_levs"
    Element 1: line wavelengths and intensities normalized to maximum intensity.
            
                                                USER INPUTS:            
    fluxes: Integrated flux in W/m^2; calculated with "grab_integrated_fluxes" function
    obj_vel: object heliocentric velocity in m/s; positive is defined AWAY from sun
    solar_dist: distance in meters from sun
    raw_lines: string array of NIST line data. Wavelengths must be in nm
    raw_levs: string array of NIST level data. Ensure both lines and levels have no headers. Energies must be in cm^-1
    ritz_col: column for ritz wavelengths (in VACUUM)
    lower_col: column of lower energy in lines array. Note 0 indexing, and energies must be in cm^-1
    upper_col: column of upper energy. Same limitations as 'lower_col'
    aval_col: column for A values in lines array.
    temp: optional; temperature of blackbody; assumed 5777K. Only used when radfield does not have flux at required wavelengths
    supp: optional; assumed True. Supplements flux array w/ blackbody for missing wavelengths. If unused, missing
    transition wavelengths in radfield are assigned rates of 0 for all radiative processes    
    renorm: optional; assumed True. Normalizes line intensities such that the maximum intensity = 1.    
    """
    Rsun = 6.957e8 # meters
    h = 6.62607004e-34
    hbar = 1.054571817e-34
    k = 1.38064852e-23
    c = 299792458
    cm_to_J = 1.986491855e-23
    #Define arrays for lines/levels after we process them:
    lines = np.zeros((len(raw_lines[:,0]),len(raw_lines[0,:])),dtype=float)
    levs = np.zeros((len(raw_levs[:,0]),len(raw_levs[0,:])),dtype=float)
    #Save data in usable form and convert half-integer J's to floats where necessary in the lines and levels arrays:
    print('Cleaning J Values...')
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
    #First, find the LEVELS for the LINES
    # Array 'idx' is the same length as the lines array and will hold the indexes corresponding to the index of the 
    # lower level (col0) and upper level (col1)
    rel_levs = np.zeros((len(levs[:,0]-1,),3))
    for i in range(0,len(rel_levs[:,0])):
        rel_levs[i,0] = i
    #We go through and match transition energies to those in NIST
    #At the same time, we keep track of which levels have transitions w/ A vals and generate the rates (product of einstein coeffs. and rad. field)
    #If a level has a transition, we store a value '1' in the second column of "rel_levs" array
    idx = np.zeros((len(lines[:,0]),7))
    print('Matching Transitions to Levels and generating rates...')
    print(' !!! WARNING !!! May take several minutes if using large spectral atlas! ')
    datagen_start = time.time()
    for i in range(0,len(lines[:,0])):
        #Check we are not on a header line; otherwise, proceed!
        upper_lev = lines[i,upper_col]
        lower_lev = lines[i,lower_col]
        #find level indices:
        for j in range(0,len(levs[:,0])):
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
        freq = (upper_energy - lower_energy)*cm_to_J/h
        #Degeneracies and Avalue:
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))
        g_i = 2 *float(levs[int(idx[i,0]),2]) + 1 #lower level degeneracy
        g_j = 2 *float(levs[int(idx[i,1]),1]) + 1 #upper level degeneracy
        aval = float(lines[i,aval_col]) #spont. emission rate
        omega = freq * 2 * np.pi #for stim. coeff.
        stim_coeff = aval * np.pi**2 * c**3 / (hbar * omega**3)
        absorp_coeff = stim_coeff * g_j / g_i     #calc'd absorption coeff. from stim. coeff. This is regardless of form of rad. field
        #Calculate rates from integrated fluxes that were pre-calculated:
        rad_at_wav = fluxes[i]
        absorp_rate = absorp_coeff * rad_at_wav
        stim_rate = stim_coeff * rad_at_wav #
        #
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
    print('Total time generating rates: {:} s'.format(datagen_end - datagen_start))
    print('Total Number of Relevent Levels: {:} out of {:}'.format(new_idx,len(levs[:,0])))
    """
    The final array will take the form [X][n] = 0, which we solve via linear algebra
    #It has dimensions new_idx x new_idx (i.e. spans the space of 'relevant' levels):
    """
    lhs = np.zeros((new_idx,new_idx)) #[X] = lhs
    pop_start = time.time()
    print('Starting Matrix Construction...')
    for i in range(0,len(lines[:,0])):
        lower_idx = int(idx[i,0]) #lower level indx for line lines[i,0], belonging to level index idx[i,0] in nist arr
        upper_idx = int(idx[i,1]) #similarly, upper level indx
        low_idx_mapped = int(rel_levs[lower_idx,2])
        up_idx_mapped = int(rel_levs[upper_idx,2])
        aval = idx[i,4]
        absorp_rate = idx[i,5]
        stim_rate = idx[i,6]
        """
        Upper level denoted by j, lower by i:
        Spont. emission, A_j->i values populate two cells: positive term to [i,j], negative to [j,j]
        Stim. emission, B_j->i values populate two cells: positive term to [i,j], negative to [j,j]
        Absorption, B_i->j values populate two cells: negative to [i,i], positive to [j,i]
        
        """
        #A Value; 2 locations. lower level:
        lhs[low_idx_mapped,up_idx_mapped] = lhs[low_idx_mapped,up_idx_mapped] + aval #Correct
        lhs[up_idx_mapped,up_idx_mapped] = lhs[up_idx_mapped,up_idx_mapped] - aval #Correct
        #stimulated emission:
        lhs[low_idx_mapped,up_idx_mapped] = lhs[low_idx_mapped,up_idx_mapped] + stim_rate #Correct
        lhs[up_idx_mapped,up_idx_mapped] = lhs[up_idx_mapped,up_idx_mapped] - stim_rate #Correct
        #absorption:
        lhs[low_idx_mapped,low_idx_mapped] = lhs[low_idx_mapped,low_idx_mapped] - absorp_rate #Correct
        lhs[up_idx_mapped,low_idx_mapped] = lhs[up_idx_mapped,low_idx_mapped] + absorp_rate #Correct
    pop_end = time.time()
    print('Rate Matrices Populated in {:} seconds'.format(pop_end - pop_start))
    print('Starting Matrix Solution...')
    
    #Calculate Eigenvectors of A^T * A, and use the eigenvector corresponding to the smallest eigenvalue:
    #This is equivalent to using scipy or numpy singular value decomposition (SVD) methods:
    eigenvalues, eigenvectors = np.linalg.eig(np.dot(np.transpose(lhs),lhs))
    #
    # Following code added on 01-13-2021
    # When using the 
    #
    pops = eigenvectors[:, np.argmin(eigenvalues)]
    if (pops[0] < 0):
        #Since any scalar multiple of an eigenvector is also an eigenvector, we force pops positive by default.
        pops = pops * (-1)
    #Confirm that populations are all positive. If not, you have a big problem:
    for i in range(0,len(pops[:])):
        if (pops[i] < 0):
            print('!!! WARNING !!! Negative population for level {:}'.format(i))
            print('Possible cause: bad atomic data, or bad blackbody approximation. Consider setting supp=False in function call!')
    solution_end = time.time()
    print('Matrix solution Completed in {:} seconds'.format(solution_end - pop_end))   
    """
    Now, we need to calculate intensities. The intensity of a line is defined as:
       I = N hbar omega Aval
       since hbar*omega is energy, we'll just calculate the energy (in Joules) from the vacuum wavelength
       E = hc/lamba
    """
    print('Starting Intensity Calculation...')
    intens = np.zeros((len(lines[:,0]),3))   
    #We'll save it as: (0) wavelength in air, (1) wavelength in vac, (2) intensity (arb. units.)
    for i in range(0,len(intens[:,0])):
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        upper_pop = pops[mapped_index] #Grab the (normalized) population
        energy = h * c / (idx[i,3] * 1e-9) #Calculate transition energy from VAC wavelength.
        line_intens = upper_pop * energy * lines[i,aval_col] #calculate line intensity from above
        intens[i,0] = idx[i,2] #Save (air) wavelength
        intens[i,1] = idx[i,3] #Save vac wavelength
        intens[i,2] = line_intens #save intensity
    #Normalize to max line = 1
    if (renorm == True):
        max_intens = np.nanmax(intens[:,2])
        intens[:,2] = intens[:,2] / max_intens
    intens_end = time.time()
    print('Intensity Calculations finished in {:} seconds'.format(intens_end - solution_end))
    return(pops,intens)

