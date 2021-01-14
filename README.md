# fluorescence_model

The newest version of the model code is contained in "fluor_v6.py"
Several problems in previous versions were located and corrected:

-Improper integration over solid angle

-Blackbody function unit issues

-Solar flux handling was moved to a separate function; now, the flux is calculated prior to the main function call and passed in, allowing a user to re-run the model
without requiring several minutes to generate the integrated fluxes used to generate absorption/stimulated emission rates.

-Files to generate error uncertainties described in the .pdf are being prepared.

----
                                                      Downloading Solar Spectra and Running the Test Script
The test script for Ni will generate synthetic spectra at 1.02 AU for a comet traveling at -36.5 km/s w.r.t. the sun
To run the model, do the following:

1. Download all files in this github repository,
2. The code requires Kurucz solar spectra; you can use a computed or measured spectra. If a wavelength is present in the lines files but absent in the imported solar spectrum,
the code defaults to using a blackbody to approximate the flux at those wavelengths. Tests for Ni I using both computed and measured solar spectra show some variation in line height, but first tests show scatter comparable to the uncertainty introduced by the Einstein A value uncertainties (see pdf for error estimate procedure). 

The attached script uses both the measured and computed spectra if the user is interested in comparing. 
The computed Kurucz solar spectra can be found at
http://kurucz.harvard.edu/stars/sun/
and download the file marked: fsunallp.100000resam3   13-Apr-2011 12:03   46M  
The test script requires this file be saved as 'kurucz_solar_spectra.txt'

\textbf{TEST}

The measured Kurucz solar spectra can be found at 
http://kurucz.harvard.edu/sun/irradiance2005/irradrelwl.dat

You want to download this file and process it to 2 column format with col0 = wavelength and col1 = flux; for this, you can use excel, Python, or whichever your preferred tool is. Don't change the units or make any adjustments to the values in these files; the Ni I script shows what units are needed and manipulates as necessary. Both files are too large to store on github.

3. With all of the files in the same directory, run the script "Ni I spectra v6.py" It should take ~6min (3min per solar spectra). 
4. A separate script, "Ni I spectra vs Hyakutake v6.py" is built on the script in step 3 and will plot the synthetic spectra vs spectra from comet Hyakutake if the user has that data.

----
                                      To run the model for a different species, do the following:
To run a different system (say Ni II), download the line lists and energy levels from NIST:
1. Query https://physics.nist.gov/PhysRefData/ASD/lines_form.html for lines of e.g. "Ni II", and under "Advanced Options" ask for only lines with transition probabilities. Get the data in tab-delimited form. Save as a .txt.
2. Open the .txt file in excel. Note that if the J values are half-integers, excel defaults to reading them as dates. You want to set those columns to "text". Double-check that excel has opened them correctly. By opening this file in excel, excess "" characters are removed automatically.
3. Remove all headers; depending on the species, there may be multiple: (1) at top of page for vacuum wavelengths, (1) at the beginning of air wavelengths, and (1) around 2000nm to indicate the switch to vacuum wavelengths for IR transitions. Remove these headers if present.
4. Download the levels for that system from https://physics.nist.gov/PhysRefData/ASD/levels_form.html
5. Save as tab-delimited, and repeat the same procedure as the lines. Take care with J values, remove headers, and resave. These files will then be passed into the model code in the "raw_lines" and "raw_levs" variables.
6. Change the example script over to the relevant filenames for the lines/levels file and change the "species_string" variable to the appropriate species. This variable is not critical for the model, but determines the names of the output files when the synthetic spectra plot and model outputs are saved.

