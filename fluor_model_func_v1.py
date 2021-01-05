# -*- coding: utf-8 -*-
"""
Fluorescence Model using NIST data
SJB
Released tentatively on 01-04-2021
------------------------
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from random import random
import time
from scipy.optimize import nnls
import os

def doppler(wavelength,velocity):
    #positive velocity defined as away from sun
    c = 299792458
    vel = velocity #assuming m/s
    rest_wavelength = wavelength
    top = 1 + (velocity/c)
    bottom = 1 - (velocity/c)
    frac_term = math.sqrt(top/bottom)
    shifted_wavelength = rest_wavelength * frac_term
    return(shifted_wavelength)

def blackbody_wave(wavelength_range, temp):
    #temp in Kelvin, all in SI
    #wavelength_range must be in METERS. Assuming defined in nm, we convert:
    wavelength_range = wavelength_range/(1e9)
    h = 6.62607004e-34
    c = 299792458
    k = 1.38064852e-23
    exp_term = 1/(np.exp(h*c/(wavelength_range * k * temp)) - 1)
    intensity = 2*h*c**2 / (wavelength_range**5) * exp_term
    #convert to ENERGY DENSITY: MULTIPLY BY 4pi / c:
    intensity = 4* np.pi * intensity / c    
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
    #NOTE: inputs assumed in nm, which we convert to angstroms for the IAU standard conversion:
    wavelength_vac = wavelength_vac * 10 #convert to ang
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength_air = wavelength_vac / n
    wavelength_air = wavelength_air / 10 #convert back to nm
    return(wavelength_air)

def convert_to_float(frac_str):
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

h = 6.62607004e-34
hbar = 1.054571817e-34
k = 1.38064852e-23
c = 299792458
cm_to_J = 1.98630e-23

def fluorescence_spectra(radfield,obj_vel,solar_dist,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col,temp,lower_cutoff,upper_cutoff,supp):
    """
                USER INPUTS
    ritz_col: column (note: python is 0 indexed) of ritz wavelengths in NIST lines file
    lower_col: column for lower level ENERGY (in inverse cm)
    upper_col: column for upper level ENERGY (in inverse cm)
    aval_col: column for Einstein A values
    obj_vel: comet/object heliocentric velocity in meters for person (-36.7 x 1e3 [m/s] for Hyakutake)
    raw_lines: loaded in NIST lines data as strings
    raw_levels: loaded in NIST levels data as strings. 
    solar_dist: distance from sun in meters; assumed 
    radfield: high resolution Kurucz solar spectrum
    supp: if 'True', supplements the imported solar spectra by using a blackbody intensity value if a given lines wavelength
    is outside the range of values in the imported solar spectrum
    
    This function calculates synthetic spectra and returns a tuple. The tuple elements are:
    Element 0: level populations, indexed the same as "raw_levs"
    Element 1: line wavelengths and intensities normalized to maximum intensity.
    """
    Rsun = 6.957e8 # meters
    h = 6.62607004e-34
    hbar = 1.054571817e-34
    k = 1.38064852e-23
    c = 299792458
    cm_to_J = 1.98630e-23
    #The end result is W/m^3 (W/m**2 per wavelength in meters)
    #Define arrays for lines/levels after we adjust them:
    lines = np.zeros((len(raw_lines[:,0]),len(raw_lines[0,:])),dtype=float)
    levs = np.zeros((len(raw_levs[:,0]),len(raw_levs[0,:])),dtype=float)
    #Clean up J values. We convert half-integer J's to floats where necessary in the lines and levels arrays:
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
        #Swap in Ritz wavelength (air) where necessary:
        if (math.isnan(lines[i,0]) == True):
            lines[i,0] = lines[i,ritz_col] #set to Ritz value is no observed wavelength, e.g. lines[i,0] is nan
        if (int(lines[i,0]) == 0):
            lines[i,0] = raw_lines[i,ritz_col]
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
        wavelength_vac = 1e9 * h * c / (cm_to_J * (upper_energy - lower_energy))#E = hc/lambda -> lambda = hc / E
        if (wavelength_vac > 200):
            wavelength_air = vac_to_air(wavelength_vac) #nm
        else:
            wavelength_air = wavelength_vac
        freq = (upper_energy - lower_energy)*cm_to_J/h
        #Degeneracies and Avalue:
        g_i = 2 *float(levs[int(idx[i,0]),2]) + 1 #lower level degeneracy
        g_j = 2 *float(levs[int(idx[i,1]),1]) + 1 #upper level degeneracy
        aval = float(lines[i,aval_col]) #spont. emission rate
        omega = freq * 2 * np.pi #for stim. coeff.
        stim_coeff = aval * np.pi**2 * c**3 / (hbar * omega**3)
        absorp_coeff = stim_coeff * g_j / g_i     #calc'd absorption coeff. from stim. coeff. This is regardless of form of rad. field
        #mapped index from all levs -> relevant levs in lhs and rhs arrays:        
        #The user can 'supplement' the imported solar spectra w/ a blackbody if Supp=True and wavelength is outside range
        if ((supp == True) and ((wavelength_vac < radfield[0,0]) or (wavelength_vac > radfield[-1,0]))):
            #If true, we calculate the blackbody (wavelength form) flux at the doppler-shifted wavelength.
            rad_at_wav = blackbody_wave(doppler(wavelength_vac,obj_vel),temp)*4*np.pi/c * Rsun**2 / solar_dist**2      
            absorp_rate = absorp_coeff * rad_at_wav
            stim_rate = stim_coeff * rad_at_wav #
        else:
            if (wavelength_vac > radfield[0,0]):
                #If not true, we find and use the flux from Kurucz at the doppler-shifted wavelength
                rad_at_wav = radfield[find_nearest_index_mod(radfield[:,0],doppler(wavelength_vac,obj_vel)),1]
                absorp_rate = absorp_coeff * rad_at_wav
                stim_rate = stim_coeff * rad_at_wav #
            if (wavelength_vac < radfield[0,0]):
                aval=0
                absorp_rate = 0
                stim_rate = 0
        #Save relevant data to idx array:
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
    print('Total Number of Relevent Levels: {:} out of {:}'.format(new_idx,len(levs[:,0])-1))
    """   
        Now we set up the array that we will invert to find the equilibrium level populations.    
        We have [X][n] = [dn/dt]
        Where we set the RHS to 0 by requiring all dn/dt = 0 in equilibrium. Rearranging, we have
        the final array will take the form [X][n] = Y, where X contains all 'populating' temrs and Y contains all 'depopulating' terms.
    #It has dimensions #new_idx by new_idx (i.e. spans the space of 'relevant' levels):
    """
    lhs = np.zeros((new_idx,new_idx)) #[X] = lhs
    #The second array is the RHS above, e.g. [Y]; it has 1 col and length #relevant levels:
    rhs = np.zeros((new_idx,1))
    """
    Let's do an example. Say we have a line that goes from level 35 to level 12. The transition rate is A_(35 -> 12)
    The COLUMN here would be col 35, as the effective rate is (A(35->12) * n(35))
    The A value can be thought of as belonging to (A_ji) where j = col, i = row, and col/row given by upper(col) and lower(row) indices
    The Avalue above would be POSITIVE for populating lower level 12; similarly, the Aval would also be placed
    on the RHS in the 35th ROW for depopulating level 35.
    A similar thing occurs for the absorption. Absorption LOWERS the population of lower state, and increases pop
    of upper state, while stimulated emission produces the opposite.
    We now populate the matrices by placing our A vals and previously calculated stim/absorp rates in the correct locations:
    """
    pop_start = time.time()
    print('Starting Matrix Construction...')
    # Updated version to only include relevant levels:
    for i in range(0,len(lines[:,0])):
        #Only use air wavelengths (for now):
        if (lower_cutoff < lines[i,0] < upper_cutoff):
            lower_idx = int(idx[i,0]) #lower level indx for line lines[i,0], belonging to level index idx[i,0] in nist arr
            upper_idx = int(idx[i,1]) #similarly, upper level indx
            low_idx_mapped = int(rel_levs[lower_idx,2])
            up_idx_mapped = int(rel_levs[upper_idx,2])
            aval = idx[i,4]
            absorp_rate = idx[i,5]
            stim_rate = idx[i,6]
            #A_ji --> 
            lhs[low_idx_mapped,up_idx_mapped] = lhs[low_idx_mapped,up_idx_mapped] + aval #Correct
            rhs[up_idx_mapped] = rhs[up_idx_mapped] + aval #Correct
            #stimulated emission:
            lhs[low_idx_mapped,up_idx_mapped] = lhs[low_idx_mapped,up_idx_mapped] + stim_rate #correct
            rhs[up_idx_mapped] = rhs[up_idx_mapped] + stim_rate #correct
            #absorption:    
            rhs[low_idx_mapped] = rhs[low_idx_mapped] + absorp_rate #correct
            lhs[up_idx_mapped,low_idx_mapped] = lhs[up_idx_mapped,low_idx_mapped] + absorp_rate #correct now (4PM 12-15-2020)
    pop_end = time.time()
    print('Rate Matrices Populated in {:} seconds'.format(pop_end - pop_start))
    print('Starting Matrix Solution...')
    pops = nnls(np.array(lhs),np.concatenate(rhs))[0]
    #Renormalize w.r.t. total population
    total_pop = np.sum(pops[:])
    pops = pops / total_pop
    solution_end = time.time()
    print('Matrix solution Completed in {:} seconds'.format(solution_end - pop_end))   
    """
    Now, we need to calculate intensities. The intensity of a line is defined as:
       I = N hbar omega Aval
       since hbar*omega is energy, we'll just calculate the energy (in Joules) from the wavelength!
       E = hc/lamba
    """
    print('Starting Intensity Calculation...')
    intens = np.zeros((len(lines[:,0]),2))   
    #We'll save it as: (0) wavelength, (1) intensity (arb. units.)
    for i in range(0,len(intens[:,0])):
        wavelength = idx[i,2]
        upper_pop_index_unmapped = int(idx[i,1])
        mapped_index = int(rel_levs[upper_pop_index_unmapped,2]) #the line upper level index is in the original 'NIST' mapping. We need the 'relevant'
        #index that we stored in the rel_levs array.
        upper_pop = pops[mapped_index] #Grab the (normalized) population
        energy = h * c / (wavelength * 1e-9) #Calculate transition energy from air wavelength.
        line_intens = upper_pop * energy * lines[i,aval_col] #calculate line intensity from above
        intens[i,0] = wavelength #Save (air) wavelength
        intens[i,1] = line_intens #Save intensity.
    #rescale to maximum wavelength
    #Normalize to max line = 1
    max_intens = np.nanmax(intens[:,1])
    intens[:,1] = intens[:,1] / max_intens
    intens_end = time.time()
    print('Intensity Calculations finished in {:} seconds'.format(intens_end - solution_end))
    return(pops,intens)
