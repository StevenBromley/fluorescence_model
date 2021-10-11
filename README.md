The code contained in this repo is the same version used in Bromley et al 2021 in PSJ (see also arxiv): "Atomic iron and nickel in the coma of C/1996 B2 (Hyakutake): production rates, emissionmechanisms, and possible parents". A theoretical description is available in Bromley et al 2021. 

The code is distributed under a GNU GPL license. If you intend to modify this code, though not required, it is prefered that you contact me for help or recommendations in updating/modifying the codes.

The code is built to use atomic data from the NIST ASD:
https://physics.nist.gov/PhysRefData/ASD/lines_form.html

The lines and levels files are generated by quering the NIST ASD for lines (only those with transition rates), and energies in units of cm^-1. When the lines/levels pages are loaded as Tab-Delimited, simply right click and "Save As" your preferred naming convention. The names of these files are used in the script and passed into a custom data-loading function (see example script, described below).

A sample file, "ni0 example script.py" shows the syntax and order of function calls. This sample script will general g-factors for neutral nickel at 1 AU (heliocentric velocity = 0 km/s).

The syntax for performing Monte-Carlo iterations to generate approximate model uncertainties with 10 iterations is provided. Note that 10^4 or more iterations
are required for good convergence.

The solar spectrum assumed in the model is provided as a compressed file (uncompressed it is ~200mb). The default "0-2730_vac.txt" file, contains the following:
0 - 168 nm: The UV coronol/chromosphere model of Fontenla et al (2014) [1]. This spectrum shows good agreement with measurements of the quiet sun spectrum from the SORCE mission.
168 - 200.06 nm: SOLSTICE quiet-sun measurements with 0.1 nm resolution
200.06 - 202 nm: Hall and Anderson (1991) UV measurements [2]
 202 - 2730 nm: Composite high-resolution benchmarked solar spectrum from Coddington et al (2021). Their dataset is compiled
          from various space- and ground-based measurements of the solar spectrum, including the famous measurements of Kurucz.
          
The code is currently capable of assuming two line profiles when computing absorption/stimulated emission rates: a delta function, and a doppler (thermal) profile.

A detailed User Guide is in the works, as well as improved functionality to use approximate data for e.g. collisional quenching. 

Be advised that the code DOES NOT check that all lines are connected to ground; that is up to the user. As a short-term workaround, it is recommended to set the "upper_cutoff" variable to the ionization potential value in cm^-1 so that auto-ionizing levels are not included. Auto-ionizing level data is typically lacking, and will not have a direct or indirect radiative pathway to ground, leading to failures in the matrix solver. 
