# -*- coding: utf-8 -*-
"""
Running Fluorescence Model for Ni I
SJB
01-28-2021
This uses version 7 of the model
"""
from fluor_v8 import *
"""
fluor_v7 contains a number of dependencies and other useful packages:
    numpy (np), matplotlib.pyplot (plt), os, math
    Numerical integratin is handled internally using a simple rectangle code I have written. 
"""
model = np.genfromtxt('Ni_I_synthetic.txt', delimiter = '\t')
#model columns: col0 air wavelength, col1 vacuum wavelength, col2 intensity (arb. units)
hyak0 = np.genfromtxt('offset_0_arcsec.tab',usecols = (0,1))
#Let's try using the 345.86 line to normalize the model:
norm_line = 345.86
#Find index of that line in the relevant arrays:
norm_model_idx = find_nearest_index(model[:,0],norm_line)
norm_comet_idx = find_nearest_index(hyak0[:,0]/10,norm_line)
#Find INTENSITY at that index in the relevant arrays:
norm_int_model = model[norm_model_idx,2]
norm_int_comet = hyak0[norm_comet_idx,1]
#Scale factor from model -> comet:
norm_factor = norm_int_comet / norm_int_model 


#%%
plt.clf()
plt.plot(hyak0[:,0]/10,hyak0[:,1], color = 'green', label = 'Hyakutake 0 Offset')
plt.stem(model[:,0], model[:,2] * norm_factor, label = 'Synthetic Ni I spectra')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity (arb. units)")
plt.xlim(300,600)
