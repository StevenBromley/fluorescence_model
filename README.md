# fluorescence_model

The newest version of the model code is contained in "fluor_v8.py"
Several problems in previous versions were located and corrected:

-(02-02-2021): Fixed issue with accidental selection of trivial selection after matrix solving. 0 values are now avoided, and eigenvalues/vectors are also returned from the fluorescence_spectra function. Future releases will change the intensity calculation to a separate function to enable user testing of different eigenvectors and their effect on output spectra.

- (New in v8): Fixed an issue with distance scaling, and corrected expression for stimulated emission coefficient.

-A single solar spectrum is now provided to give the user the 'best of both worlds'; it now uses measured solar spectra, and computed high-resolution spectra elsewhere

-Added in handling of line profiles during flux integration; currently, only the 'Doppler' profile is available and is set as the default. For nearly all cases, the natural linewidth is insignificant, and the gaussian broadening dominates. Voigt and Lorentzian profiles will be added in a future release.

-Implemented new search algorithm when generating integrated fluxes; compute time down by factor of ~40. (now ~5s instead of ~3min). 
----
                                      Running the model
The test script for Ni will generate synthetic spectra at 1.02 AU for a comet traveling at -36.7 km/s w.r.t. the sun
To run the model, do the following:

1. Download all files in this github repository,

2. The code requires a solar spectrum with the following format: Col0 (wavelength; nm), Col1(flux; W / m^2 / nm). The provided file "kurucz_150nm-81um.txt" provides fluxes in the range 150nm out to 81 um from the following:

       a. Kurucz high-resolution *computed* spectra from 150 - 300nm, 
  
       b. Kurucz high-resolution *measured* spectra from a FTS spanning 300 - 1000nm;
  
       c. Kurucz high-resolution *computed* spectra from 1 um - 81 um. 
  
       d. For wavelengths outside these bounds, a blackbody is used to estimate the flux. 
  
3. Set the parameters for the model. The important parameters are described below:

    (A) solar_dist: heliocentric distance of the emitter; end result must be in meters
  
    (B) obj_velocity: heliocentric velocity (defined positive away from sun). Must be in m/s
    
    (C) m_species: mass of the emitting species in kg; used for doppler broadening profiles
  
    (D) t_comet: temperature of the emitting species in Kelvin; used for doppler profiles
  
   (E) Column IDs: 
  
         1. ritz_col: column for ritz wavelengths in the lines array
      
         2. lower_col: column for lower level energy in lines array
      
         3. upper_col: column for upper level energy in lines array
      
         4. aval_col: column for A values in lines array
      
     (F) raw_lines: the processed array of line information. Ensure all data is free of extra "" characters and all headers (there may be up to 3) have been removed.
  
    (G) raw_levs: the processed array of level information. Ensure all data is free of extra "" characters and all headers (there may be multiple) have been removed. 
  

----
                                      To run the model for a different species, do the following:
To run a different system (say Ni II), download the line lists and energy levels from NIST:

1. Query https://physics.nist.gov/PhysRefData/ASD/lines_form.html for lines of e.g. "Ni II", and under "Advanced Options" ask for only lines with transition probabilities. Get the data in tab-delimited form. Save as a .txt.

2. Open the .txt file in excel. Note that if the J values are half-integers, excel defaults to reading them as dates. You want to set those columns to "text". Double-check that excel has opened them correctly. By opening this file in excel, excess "" characters are removed automatically.

3. Remove all headers; depending on the species, there may be multiple: (1) at top of page for vacuum wavelengths, (1) at the beginning of air wavelengths, and (1) around 2000nm to indicate the switch to vacuum wavelengths for IR transitions. Remove these headers if present.

4. Download the levels for that system from https://physics.nist.gov/PhysRefData/ASD/levels_form.html

5. Save as tab-delimited, and repeat the same procedure as the lines. Take care with J values, remove headers, and resave. These files will then be passed into the model code in the "raw_lines" and "raw_levs" variables.

6. Change the example script over to the relevant filenames for the lines/levels file and change the "species_string" variable to the appropriate species. This variable is not critical for the model, but determines the names of the output files when the synthetic spectra plot and model outputs are saved.

7. Adjust other variables according and run the entire script.

