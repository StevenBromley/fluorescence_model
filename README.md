# fluorescence_model

README Updated 06-04-2021
-- Newest (and recommended) version is fluor_v13_3.

Updates:
-- Broke handling of radiation field into 3 parts to improve handling of vacuum / air wavelengths. Solar spectra now come in 3 parts: vaccum uv 150 - 200 nm, "visible" 200 - 1000 nm, and vacuum infrared 1 - 81 microns. 

                                      Running the model
The test script for Ni will generate synthetic spectra at 1.02 AU for a comet traveling at -36.7 km/s w.r.t. the sun
To run the model, do the following:

1. Download all files in this github repository,

2. The code requires a solar spectra(um) with the following format: Col0 (wavelength; nm), Col1(flux; W / m^2 / nm). The provided files "kurucz_xxx.txt" provide fluxes in the range 150nm out to 81 um from the following:

       a. Kurucz high-resolution *computed* spectra from 150 - 199.935 nm (in vacuum wavelengths). *Computed* fluxes, converted to air wavelengths, are used for 199.935 - 299 nm.
       b. Kurucz high-resolution measured spectra for used for 299 - 1001.27 nm.
  
       c. Kurucz high-resolution *computed* spectra from 1001.27 nm - 81 um. 
  
       d. For wavelengths outside these bounds, a blackbody is used to estimate the flux. 
  
3. Set the parameters for the model. The important parameters are described below:

    (A) solar_dist: heliocentric distance of the emitter; end result must be in meters
  
    (B) obj_velocity: heliocentric velocity (defined positive away from sun). Must be in m/s.

    (C) Column IDs: 
  
         1. ritz_col: column for ritz wavelengths in the lines array
      
         2. lower_col: column for lower level energy in lines array
      
         3. upper_col: column for upper level energy in lines array
      
         4. aval_col: column for A values in lines array
      
     (D) raw_lines: the processed array of line information. Ensure all data is free of extra "" characters and all headers (there may be up to 3) have been removed.
  
    (E) raw_levs: the processed array of level information. Ensure all data is free of extra "" characters and all headers (there may be multiple) have been removed. 
  

----
                                      To run the model for a different species, do the following:
To run a different system (say Ni II), download the line lists and energy levels from NIST:

1. Query https://physics.nist.gov/PhysRefData/ASD/lines_form.html for lines of e.g. "Ni II", and under "Advanced Options" ask for only lines with transition probabilities. Get the data in tab-delimited form and save both lines and levels as .txt.

2. Open the .txt files in excel. Note that if the J values are half-integers, excel defaults to reading them as dates. You want to set those columns to "text". Double-check that excel has opened them correctly. By opening this file in excel, excess "" characters are removed automatically. Remove all headers (Lines files will have up to 3; levels may have multiple depending on if autoionizing levels are listed). Save with your desired file names, and load it into the array "raw_lines" and "raw_levs" (see example script for proper genfromtxt options).

3. Adjust variables (column variables) accordingly and run.
