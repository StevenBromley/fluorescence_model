# -*- coding: utf-8 -*-
"""
SJB
Processing MonteCarlo for Ni0
"""
from fluor_v11 import *

#%%
cwd = os.getcwd()
raw_lines = np.genfromtxt('ni0_lines_processed.txt',delimiter='\t',dtype=str,skip_header=0)
raw_levs = np.genfromtxt('ni0_levs_processed.txt',delimiter ='\t',skip_header=0,usecols=(0,1,2,3),dtype=str)
#%%
ritz_col = 1 #required user input
lower_col = 5 #required
upper_col = 6 #required
uncert_col = 4
aval_col = 3 #required
solar_dist = (1.02 * 1.496e+11) #AU converted to meters
obj_vel =  -36.7 * 1e3 #Ikeya-Seki; ~154.5km/s ; Hyakutake: -36.7 * 1e3 #m/s
#%
m_species= 9.746267511e-26 #mass of Nickel in kg
t_comet = 100 # Kelvin; temperature used for doppler profile.
#Import the Kurucz spectra provided; units W / m^2 per nm:
kurucz_flux = np.genfromtxt('kurucz_150nm-81um.txt',delimiter ='\t',dtype=float,skip_header=1) 
#Calculate the integrated fluxes; this takes into account a doppler-broadedning line profile at temp t_comet by default
fluxes = fluxes_with_profiles(kurucz_flux,raw_lines,raw_levs,lower_col,upper_col,aval_col,obj_vel,solar_dist,t_comet, m_species)
one_run_start = time.time()
model = fluorescence_spectra(fluxes,raw_lines,raw_levs,ritz_col,lower_col,upper_col,aval_col, renorm =False)  
model_lines = model[1]
one_run_end = time.time()


#%%
num_samples = 10
err_test = error_calc(fluxes,raw_lines,raw_levs,ritz_col,uncert_col,lower_col,upper_col,aval_col,num_samples, model_lines, renorm=False)
err_processed =  error_process(err_test)
run_100_end = time.time()

print('Time to finish 1 iteration: {:}'.format(round( (one_run_end - one_run_start) ,2)))
print('Time to finish {:} iteration: {:}'.format(num_samples,round( (run_100_end - one_run_end) ,2)))
np.savetxt('ni0_{:}iterations_stats.txt'.format(num_samples),err_processed,delimiter='\t')
#%%
plt.clf()
plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 22})
plt.errorbar(err_processed[:,0],err_processed[:,2],yerr=err_processed[:,-1],fmt='o',ecolor='red',capsize=4, elinewidth=2)
plt.stem(err_processed[:,0],err_processed[:,2], label = 'Ni I synthetic spectra')
plt.grid(True)
plt.legend()
plt.xlim(150,600)
