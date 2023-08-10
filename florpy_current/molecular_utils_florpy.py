#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Utility codes related to FlorPy

Updated Aug 9 to include "separate_vib_bands" and some time-dependent output processing.

"""
import numpy as np
import time

def separate_vib_bands(flor, element_string, orbit_id):
    """
    04-27-2022
    SJB
    Function for splitting up an imported transitions file for MOLECULES into the relevant vibrational bands.
    NOTE: This function now utilizes only the lines included in the model, and thus indices have a 1:1 mapping
    with the indices in the final gfactor array.
    """
    start = time.time()
    #Define an output dict for all of the bands:
    #Grab the transition array from florpy dict structure:
    #temporary dict for saving indices for sorting:
    #OLD:
    #lines = flor[element_string]['lines_data_str']
    #Correct:
    lines = flor[element_string][orbit_id]['relevant_lines_str']
    #Unique band labels. Will be used to define dict entries for each array of transitions belonging to a band:
    #First, we need to find all the possible unique combinations of upper and lower electronic states:
    bands = {}
    #Loop over lines array:
    for i in range(0,len(lines[:,0])):    
        #Grab lower and upper electronic states from the lines array:
        lower_state = lines[i,5]
        upper_state = lines[i,10]
        lower_vib = int(lines[i,7])
        upper_vib = int(lines[i,12])
        #generate a state string, i.e. "upper_lower"
        labstring = '{:}_{:}'.format(upper_state,lower_state)
        #Same for vibrational:
        vibstring = '{:}_{:}'.format(upper_vib,lower_vib)
        #if existing, add to the indices array:
        if (labstring in bands):
            bands[labstring]['indices'] = np.append(bands[labstring]['indices'], [i], axis = 0)
            if vibstring in bands[labstring]['vib']:
                bands[labstring]['vib'][vibstring]['indices'] = np.append(bands[labstring]['vib'][vibstring]['indices'], [i], axis=0)
                bands[labstring]['vib'][vibstring]['trs'] = np.vstack((bands[labstring]['vib'][vibstring]['trs'], lines[i,:]))            
            else:
                bands[labstring]['vib'][vibstring] ={}
                bands[labstring]['vib'][vibstring]['indices'] = [i]
                bands[labstring]['vib'][vibstring]['trs'] = np.empty((0,len(lines[0,:])),dtype=object)
                bands[labstring]['vib'][vibstring]['trs'] = np.vstack((bands[labstring]['vib'][vibstring]['trs'], lines[i,:]))            
        else:
            bands[labstring] = {}
            bands[labstring]['indices'] = [i] #replacing:
            bands[labstring]['vib'] = {}
    end = time.time()
    print('Band post-processing completed in {:} seconds'.format((end-start)/60))
    flor[element_string]['bands'] = bands
    return(flor)
#
def grab_trs_in_range(flor,element_string, xlow, xhigh):
    """
    04-27-2022
    SJB
    This is a post-processor utility used to grab the transition information for a given wavelength range for plotting purposes.
    xlow and xhigh are the wavelength limits you want data between. Units: nm
    
    Required user has ran "separate_vib_bands" before running. 
    """
    lines = flor[element_string]['lines_data_str']
    bands = flor[element_string]['bands']
    relevant_trs = np.empty((0,len(lines[0,:])),dtype=object)
    temp_dict = {}

    for key in bands:
        for key2 in bands[key]['vib']:
            band_trs = bands[key]['vib'][key2]['trs']
            temp_trs =  np.empty((0,len(lines[0,:])),dtype=object)
            for i in range(0,len(band_trs[:,0])):
                if (xlow < vac_to_air(float(band_trs[i,0])) < xhigh):
                    relevant_trs = np.vstack((relevant_trs,band_trs[i,:]))
                    temp_trs = np.vstack((temp_trs,band_trs[i,:]))
                    #Following for testing purposes only:
                    #print('{:}-{:}: {:} - {:} nm N = {:} to N = {:}'.format(band_trs[i,10], band_trs[i,5], key2,round(float(band_trs[i,0]),2),band_trs[i,14], band_trs[i,9]))
                    #
                    if key not in temp_dict:
                        temp_dict[key] = {}
                        temp_dict[key][key2] = {}
                    if (key2 not in temp_dict[key]):
                        temp_dict[key][key2] ={}
            if (len(temp_trs[:,0]) > 0):
                temp_dict[key][key2] = temp_trs
    return(temp_dict)

def vac_to_air(wavelength_vac):
    #IAU Standard Conversion (Morton 2000)
    #NOTE: inputs assumed in nm, which we convert to angstroms required for the conversion:
    wavelength_vac = wavelength_vac * 10 #convert to ang
    s = 1e4 / wavelength_vac
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s**2) + (0.00015998)/(38.9 - s**2)
    wavelength = wavelength_vac / n
    wavelength_conv = wavelength / 10 #convert back to nm
    return(wavelength_conv)

def split_branches_rotational(rel_waves):
    """
    SJB
    05-02-2022
    Function for splitting the output of grab_trs_in_range function into individual branches.
    Requires numpy
    
    Ex for CO+:
    For Pi - X:
    12 branches:
        P11, Q11, R11, P22, Q22, R22, P12, Q12, R12, P21, Q21, R21
    delta N = delta J for: P11 Q11 R11 P22 Q22 R22
    P12 and R21: deltaN = -2 and +2
    Q12 and Q21 = delta N = -1 and +1
    R12 and P21 -> delta N = 0
    
    For B-X:
        P11, R11, P22, R22.
        Satellite branches: Q12, Q21 -> delta N = -1 and +1
    """
    out_dict = {}
    for key in rel_waves:
        out_dict[key] = {}
        for key2 in rel_waves[key]:
            lins = rel_waves[key][key2]
            out_dict[key][key2] = {}
            out_dict[key][key2]['lines'] = lins
            out_dict[key][key2]['branches'] = {}
            out_dict[key][key2]['branches']['Q'] = {} 
            out_dict[key][key2]['branches']['Q']['all_lines'] = np.empty((0,len(lins[0,:])))
            out_dict[key][key2]['branches']['Q']['split'] = {}
            out_dict[key][key2]['branches']['P'] = {} 
            out_dict[key][key2]['branches']['P']['all_lines'] = np.empty((0,len(lins[0,:])))
            out_dict[key][key2]['branches']['P']['split'] = {}
            out_dict[key][key2]['branches']['R'] = {} 
            out_dict[key][key2]['branches']['R']['all_lines'] = np.empty((0,len(lins[0,:])))
            out_dict[key][key2]['branches']['R']['split'] = {}
            #
            for i in range(0,len(lins[:,0])):
                if ('P' in lins[i,15]):
                    out_dict[key][key2]['branches']['P']['all_lines'] = np.vstack((out_dict[key][key2]['branches']['P']['all_lines'],lins[i,:]))
                if ('Q' in lins[i,15]):
                    out_dict[key][key2]['branches']['Q']['all_lines'] = np.vstack((out_dict[key][key2]['branches']['Q']['all_lines'],lins[i,:]))
                if ('R' in lins[i,15]):
                    out_dict[key][key2]['branches']['R']['all_lines'] = np.vstack((out_dict[key][key2]['branches']['R']['all_lines'],lins[i,:]))
    #%
    for key in out_dict:
        for key2 in out_dict[key]:
            #Q:
            lins = out_dict[key][key2]['branches']['Q']['all_lines']        
            for i in range(0,len(lins[:,0])):
                label = lins[i,15]
                split_lab = (label.split('('))[0]
                split_lab2 = split_lab[1:]
                if (split_lab2 not in out_dict[key][key2]['branches']['Q']['split']):
                    out_dict[key][key2]['branches']['Q']['split'][split_lab2] = np.empty((0,len(lins[0,:])))
                else:
                    out_dict[key][key2]['branches']['Q']['split'][split_lab2] = np.vstack((out_dict[key][key2]['branches']['Q']['split'][split_lab2],lins[i,:]))
            #P:
            lins = out_dict[key][key2]['branches']['P']['all_lines']        
            for i in range(0,len(lins[:,0])):
                label = lins[i,15]
                split_lab = (label.split('('))[0]
                split_lab2 = split_lab[1:]
                if (split_lab2 not in out_dict[key][key2]['branches']['P']['split']):
                    out_dict[key][key2]['branches']['P']['split'][split_lab2] = np.empty((0,len(lins[0,:])))
                else:
                    out_dict[key][key2]['branches']['P']['split'][split_lab2] = np.vstack((out_dict[key][key2]['branches']['P']['split'][split_lab2],lins[i,:]))
            #R:
            lins = out_dict[key][key2]['branches']['R']['all_lines']        
            for i in range(0,len(lins[:,0])):
                label = lins[i,15]
                split_lab = (label.split('('))[0]
                split_lab2 = split_lab[1:]
                if (split_lab2 not in out_dict[key][key2]['branches']['R']['split']):
                    out_dict[key][key2]['branches']['R']['split'][split_lab2] = np.empty((0,len(lins[0,:])))
                else:
                    out_dict[key][key2]['branches']['R']['split'][split_lab2] = np.vstack((out_dict[key][key2]['branches']['R']['split'][split_lab2],lins[i,:]))
    #Get the # of branches to plot:
    num_branches = 0
    #This dict structure is... hideous.
    for key in out_dict:
        for key2 in out_dict[key]:
            for key3 in out_dict[key][key2]['branches']:
                for key4 in out_dict[key][key2]['branches'][key3]['split']:
                    num_branches = num_branches + 1
    print('{:} Branches to plot'.format(num_branches))
    return(out_dict)

#%%
def gen_band_lum(flor,element_string,orbit_id,up_state,low_state,vib_up,vib_low):
    """
    Sep 12 2022
    SJB
    
    This function takes in the florpy dict. Requires the code "separate_vib_bands" be run prior to this call.
    The user provides and upper and lower electronic state, and upper and low vibrational levels, ex:
        up_state = 'Pi'
        low_state = 'X'
        vib_up = 2
        vib_low = 0
    The function returns the summed band luminosity for all transitions in this band in J / s / mol
    Note: The argument labels MUST match the nomenclature in the input files.
    
    """
    vibkey = '' #Starting values used to trigger useful feedback following a bad prompt to the function
    state_key = '' #Starting values used to trigger useful feedback following a bad prompt to the function
    #
    for key in flor[element_string]['bands']:
        band_id = key.split('_')
        upper_state = band_id[0]
        lower_state = band_id[1]
        if ( (upper_state == up_state) and (lower_state == low_state) ):
            state_key = key
            break
    #Now that we know the right electronic band key, we search for the correct vibrational dict:
    if (state_key == ''):
        print('States in grab_band_lum() call, Upper = {:} and Lower = {:}, could not be identified in dictionary. Please check state labels.'.format(up_state,low_state))
    else:    
        for key in flor[element_string]['bands'][state_key]['vib']:
            band_vibs = key.split('_')
            upper_vib = int(band_vibs[0])
            lower_vib = int(band_vibs[1])
            if ( (lower_vib == vib_low) and (upper_vib == vib_up)   ):
                vibkey = key
                break
        #Now, we loop over and generate the ratio:
        #If we cant find the vibrational levels, give useful error message. Otherwise return the summed band luminosity.
        if (vibkey == ''):
            print('Vibrational band for {:} (v = {:}) -> {:} (v = {:}) could not be found. Please check grab_band_lum() input'.format(up_state,vib_up,low_state,vib_low))
        else:
            sum_lum = 0
            for i in range(0,len(flor[element_string]['bands'][state_key]['vib'][vibkey]['indices'])):
                #Loop over the indices array and use those to pick out the gfactors to be added:
                sum_lum = sum_lum + flor[element_string][orbit_id]['outputs']['gfactors'][int(flor[element_string]['bands'][state_key]['vib'][vibkey]['indices'][i]),2]
    #Return the summed luminosity
    return(sum_lum)
#%%
def gen_band_lum_photonunits(flor,element_string,orbit_id,up_state,low_state,vib_up,vib_low):
    """
    Sep 12 2022
    SJB
    
    This function takes in the florpy dict. Requires the code "separate_vib_bands" be run prior to this call.
    This returns the same output as "gen_band_lum", but converts to photon units instead of energy units.
    
    The user provides and upper and lower electronic state, and upper and low vibrational levels, ex:
        up_state = 'Pi'
        low_state = 'X'
        vib_up = 2
        vib_low = 0
    The function returns the summed band luminosity for all transitions in this band in J / s / mol
  
    Note: The input state labels (Pi, X, ...) are pulled from the input file. The labels MUST match the exact nomenclature in the input files.
    """
    vibkey = '' #Starting values used to trigger useful feedback following a bad prompt to the function
    state_key = '' #Starting values used to trigger useful feedback following a bad prompt to the function
    #
    h = 6.62607004e-34 #SI
    c = 299792458 #SI
    for key in flor[element_string]['bands']:
        band_id = key.split('_')
        upper_state = band_id[0]
        lower_state = band_id[1]
        if ( (upper_state == up_state) and (lower_state == low_state) ):
            state_key = key
            break
    #Now that we know the right electronic band key, we search for the correct vibrational dict:
    if (state_key == ''):
        print('States in grab_band_lum() call, Upper = {:} and Lower = {:}, could not be identified in dictionary. Please check state labels.'.format(up_state,low_state))
    else:    
        for key in flor[element_string]['bands'][state_key]['vib']:
            band_vibs = key.split('_')
            upper_vib = int(band_vibs[0])
            lower_vib = int(band_vibs[1])
            if ( (lower_vib == vib_low) and (upper_vib == vib_up)   ):
                vibkey = key
                break
        #Now, we loop over and generate the ratio:
        #If we cant find the vibrational levels, give useful error message. Otherwise return the summed band luminosity.
        if (vibkey == ''):
            print('Vibrational band for {:} (v = {:}) -> {:} (v = {:}) could not be found. Please check grab_band_lum() input'.format(up_state,vib_up,low_state,vib_low))
        else:
            sum_lum = 0
            for i in range(0,len(flor[element_string]['bands'][state_key]['vib'][vibkey]['indices'])):
                #Loop over the indices array and use those to pick out the gfactors to be added:
                en_factor = 1e-9 * (1e7 / (float(flor[element_string][orbit_id]['relevant_lines'][i,3]) - float(flor[element_string][orbit_id]['relevant_lines'][i,2]))) / (h * c)
                sum_lum = sum_lum + en_factor * flor[element_string][orbit_id]['outputs']['gfactors'][int(flor[element_string]['bands'][state_key]['vib'][vibkey]['indices'][i]),2]
    #Return the summed luminosity
    return(sum_lum)

#%%
def auto_gen_band_lum(flor,element_string,orbit_id):
    """
    Feb 1st 2023
    SJB

    This function takes in the florpy dict. Requires the code "separate_vib_bands" be run prior to this call.
    This function will take the processed bands dict from separate_vib_bands and loop over it to calculate the band
    luminosities
    
    Output dict is stored in:
        flor[element_string][orbit_id]['outputs']['band_gfacs']
    with the form
        band1
        band2
        band3
        ...
    where e.g. 'band1' might be (X-X) and contain an array (dtype is object) of the form (| indicating column separators):
    Lower state | Upper State | Upper vib num | Lower vib num | band luminosity (J/s/mol)
    """
    band_lums = {}
    bands_dict = flor[element_string]['bands']
    for key in bands_dict:
        band_id = key.split('_')
        upper_state = band_id[0]
        lower_state = band_id[1]
        #
        band_lums[key] = {}
        key_counter = 0 
        output_lums = np.empty((len(flor[element_string]['bands'][key]['vib']),5),dtype=object)
        output_lums[:,0] = upper_state
        output_lums[:,1] = lower_state
        #
        for key2 in flor[element_string]['bands'][key]['vib']:
            vib_tr_inds = flor[element_string]['bands'][key]['vib'][key2]['indices']
            vib_info = key2.split('_')
            lower_v = vib_info[1]
            upper_v = vib_info[0]
            output_lums[key_counter,2] = upper_v
            output_lums[key_counter,3] = lower_v
            sum_lum = 0
            #
            for i in range(0,len(vib_tr_inds)):
                sum_lum = sum_lum + flor[element_string][orbit_id]['outputs']['gfactors'][int(vib_tr_inds[i]),2]
            output_lums[key_counter,4] = sum_lum
            key_counter = key_counter + 1
        #
        output_lums_header = np.empty((1,5),dtype=object)
        output_lums_header[0,0] = 'Upper State'
        output_lums_header[0,1] = 'Lower State'
        output_lums_header[0,2] = 'Upper Vib'
        output_lums_header[0,3] = 'Lower Vib'
        output_lums_header[0,4] = 'Band Lum. (J/s/mol)'
        output_tot = np.append(output_lums_header,output_lums,axis=0)
        band_lums[key] = output_tot
    #
    flor[element_string][orbit_id]['outputs']['band_gfacs'] = band_lums
    return(flor)

def gen_td_band_lum(flor,element_string,orbit_id,td_id,up_state,low_state,up_vib,low_vib):
    """
    SJB
    custom function for computing instantaneous band luminosities for CO+
    Feb 7 2023
    
    The user provides and upper and lower electronic state, and upper and low vibrational levels, ex:
        up_state = 'Pi'
        low_state = 'X'
        vib_up = 2
        vib_low = 0
    """
    bands_dict = flor[element_string]['bands']
    state_key = up_state + '_' + low_state
    vib_key = str(up_vib) + '_' + str(low_vib)
    band_inds = bands_dict[state_key]['vib'][vib_key]['indices']
    time_grid = flor[element_string][orbit_id][td_id]['time_grid']
    #
    td_bandlum = np.zeros((len(time_grid[:]),2))
    td_bandlum[:,0] = time_grid
    for i in range(0,len(td_bandlum[:])):
        td_bandlum[i,1] = np.sum(flor[element_string][orbit_id][td_id]['time_dep_gfactors'][band_inds[:],i])
    return(td_bandlum)

def auto_gen_td_bandlums(flor,element_string,orbit_id,td_id):
    """
    SJB
    Feb 8, 2023
    This is a function to automatically calculate all of the time-dependent band luminosities using
    the instantaneous time-dependent band luminosities
    
    This function takes in the florpy dict. Requires the code "separate_vib_bands" be run prior to this call.
    This function will take the processed bands dict from separate_vib_bands and loop over it to calculate the band
    luminosities using the time-dependent band luminosities computed in a time-dependent FlorPy run
        
    Output dict is stored in:
        flor[element_string][orbit_id]['td']['outputs']['band_gfacs']
    with the form
        band1
        band2
        band3
        ...
    where e.g. 'band1' might be (X-X) and contain an array (the dtype is object) of the form:
    Lower state | Upper State | Upper vib num | Lower vib num | band luminosity (J/s/mol)
    """
    td_processing_start = time.time()
    td_band_lums = {}
    bands_dict = flor[element_string]['bands']
    time_grid = flor[element_string][orbit_id][td_id]['time_grid']
    #
    for key in bands_dict:
        band_id = key.split('_')
        upper_state = band_id[0]
        lower_state = band_id[1]
        #
        td_band_lums[key] = {}
        key_counter = 0 
        output_lums = np.empty((len(flor[element_string]['bands'][key]['vib']),len(time_grid) + 4),dtype=object)
        output_lums[:,0] = upper_state
        output_lums[:,1] = lower_state
        #
        for key2 in flor[element_string]['bands'][key]['vib']:
            vib_tr_inds = flor[element_string]['bands'][key]['vib'][key2]['indices']
            vib_info = key2.split('_')
            lower_v = vib_info[1]
            upper_v = vib_info[0]
            output_lums[key_counter,2] = upper_v
            output_lums[key_counter,3] = lower_v
            #Compute the time-dependent gfactors for just this band:
            #print(key_counter)
            for i in range(0,len(time_grid[:])):
                output_lums[key_counter,i+4] = np.sum(flor[element_string][orbit_id][td_id]['time_dep_gfactors'][vib_tr_inds[:],i])
            key_counter = key_counter + 1
        #
        output_lums_header = np.empty((1,len(time_grid) + 4),dtype=object)
        output_lums_header[0,0] = 'Upper State'
        output_lums_header[0,1] = 'Lower State'
        output_lums_header[0,2] = 'Upper Vib'
        output_lums_header[0,3] = 'Lower Vib'
        output_lums_header[0,4] = 'Band Lum. (J/s/mol)'
        output_lums_header[0,4:] = time_grid.reshape(-1)
        output_tot = np.append(output_lums_header,output_lums,axis=0)
        td_band_lums[key] = output_tot
        #
    flor[element_string][orbit_id][td_id]['band_gfacs'] = td_band_lums
    td_processing_end = time.time()
    print('Calculation of time-dependent band luminosities finished in {:} sec'.format(round((td_processing_end - td_processing_start)/60,2)))
    return(flor)




