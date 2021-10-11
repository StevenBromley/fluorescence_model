# -*- coding: utf-8 -*-
"""
SJB
10-11-2021
Running Fluorescence Code for neutral nickel at 1 AU
"""
from florpy_v1 import *
element_string = 'ni0'
#Two input files: lines file (tab-delimited) from NIST, and tab-delimited levels from NIST.
lines_file = 'ni0_lines_txt.txt'
levs_file = 'ni0_levs_txt.txt'
#Define orbital parameters:
solar_dist = 1.0                                        #Astronomical Units
obj_vel = 0                                             #km/s; internally converted to m/s when needed.
geo_vel = 0                                             #km/s;
geo_dist = 1                                            #geocentric distance; 1AU; not used at present.
temp_comet = 279                                        #Kelvin; comet temperature for doppler profiles.
mass = 9.746267511e-26                                  #Atomic mass of Ni in kg. Used for doppler profile. 
orbital_conds = [solar_dist,obj_vel, geo_dist,geo_vel]
orbit_id = 'test'                                       #unique identifier for the calculations in the dictionary structure
#####################################################
#                                                   #
#                 RUNNING THE MODEL                 #
#                                                   #
#####################################################
#In sequence, we run:
flor = build_model()
flor = add_element_to_model(flor, element_string)
flor = load_nist_data(flor,element_string,lines_file, levs_file)
flor = define_rad_fields(flor,element_string,rad_source='default')
flor = load_orbit_data(flor,element_string,orbit_id,orbital_conds,t_comet = temp_comet, m_species = mass)
flor = generate_fluxes_general(flor, element_string, orbit_id, rad_choice = 'default', profiles='doppler', bbsupp=True, lower_cutoff = 1e30, upper_cutoff=1e30, aval_min=0)
#The 'bbsupp' true indicates that a blackbody should be used to approximate wavelengths outside the radiation field, i.e. above 2730 for the 'default' file
# 'upper_cutoff' sets the cutoff value for the higher allowed energy level. Default values are shown such that *all* atomic data is included for this example.
#For alternative radiation field options, see documentation.
#Two line profiles are available: "delta" (dirac delta), and "doppler" broadened.
flor = calculate_pops_and_gfactors(flor,element_string,orbit_id)


lines = flor[element_string][orbit_id]['outputs']['gfactors']
#The 'lines' array contains 3 columns of data: 
#Col0 is wavelength in air nm, Col1 is wavelength in vacuum nm, col3 is gfactor in J/s/particle at the emitting source


plt.clf()
plt.figure(figsize = (8,8))
plt.stem(lines[:,0], lines[:,2])
plt.ylabel(r'g-factor (J s$^{-1}$ particle$^{-1}$)')
plt.xlabel('Air Wavelength (nm)')
plt.grid(True)
plt.xlim(0,1000)
plt.savefig('ni0_gfactors.pdf')


#####################################################
#                                                   #
#      RUNNING THE MONTE-CARLO ERROR ESTIMATION     #
#                                                   #
#####################################################

flor = error_iterations(flor,element_string,orbit_id,err_id = 'err10', num_iterations=10,profiles='doppler')
flor = error_iterations(flor,element_string,orbit_id,err_id = 'err100', num_iterations=100,profiles='doppler')
flor = error_iterations(flor,element_string,orbit_id,err_id = 'err1k', num_iterations=1000,profiles='doppler')
e1 = flor[element_string][orbit_id]['error_calcs']['err10']['processed_iteration_output']
e2 = flor[element_string][orbit_id]['error_calcs']['err100']['processed_iteration_output']
e3 = flor[element_string][orbit_id]['error_calcs']['err1k']['processed_iteration_output']








