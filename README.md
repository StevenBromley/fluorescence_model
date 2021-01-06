# fluorescence_model

EDIT: Removed code on 01-05-2021 to fix issues w/ invalid matrix operations

The model code is contained in "fluor_model_func_v1.py"
The test script for Ni will generate synthetic spectra at 1.02 AU for a comet traveling at -36.5 km/s w.r.t. the sun

The code requires a Kurucz solar spectra; to acquire this, go to:

http://kurucz.harvard.edu/stars/sun/

and download the file marked: fsunallp.100000resam3   13-Apr-2011 12:03   46M  

The test script requires this file be saved as 'kurucz_solar_spectra.txt'
