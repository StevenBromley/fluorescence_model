The code contained in this repo is the same version used in Bromley et al 2021 in PSJ (see also arxiv): "Atomic iron and nickel in the coma of C/1996 B2 (Hyakutake): production rates, emission mechanisms, and possible parents". The theoretical description is available in Bromley et al 2021, but an updated User Guide with theoretical description is in the works. We are also investigating ways to approximate data for collisional quenching involving atoms and water vapour.

The code is distributed under a GNU GPL license. If you intend to modify this code, it is prefered that you contact me for help or recommendations in updating/modifying the codes.

Lastly, I cannot confirm this code will work on all platforms. The codes have been confirmed to work on both Windows 10 and MacOS Big Sur in the Spyder IDE. For most users, I would recommend the code is executed in Spyder to facilitate viewing the various parts of the dictionary structure.

--------

The code is built to use atomic data from the NIST ASD:
https://physics.nist.gov/PhysRefData/ASD/lines_form.html

The user must acquire lines and levels files by querying the NIST ASD for lines (only those with transition rates), and energies in units of cm^-1:

-- On the lines page, enter the element of interest, e.g. "Ni I", and under "Show Advanced Settings" select "Only lines with transition probabilities. Under "Format Output" select Tab-Delimited. Once the page is loaded, right click and "Save As" with your preferred naming convention. Repeat the same for the levels (Tab-Delimited -> Save As). These two files will be required inputs for the python script that runs the model. 

A sample file, "ni0 example script.py" shows the syntax and order of function calls. This sample script will general g-factors for neutral nickel at 1 AU (heliocentric velocity = 0 km/s).

The syntax for performing Monte-Carlo iterations to generate approximate model uncertainties with 10 iterations is also provided. Note that 10^4 or more iterations are required for good convergence. 10^4 iterations will execute in around 15 - 20 minutes on compute time on a modern macbook. Be advised that the RAM required to store more runs (e.g. 10^6) at once may not be possible on a personal computer.

-------- 

The solar spectrum assumed in the model is provided as a compressed file (uncompressed it is ~200mb). The default "0-2730_vac.txt" file, contains the following:

0 - 168 nm: The UV coronal/chromosphere model of Fontenla et al (2014) [1]. This spectrum shows good agreement with measurements of the quiet sun spectrum from the SORCE mission.


168 - 200.06 nm: SOLSTICE quiet-sun measurements with 0.1 nm resolution.


200.06 - 202 nm: Hall and Anderson (1991) UV measurements [2].


202 - 2730 nm: Composite high-resolution benchmarked solar spectrum from Coddington et al (2021) [3]. Their dataset is compiled from various space- and ground-based measurements of the solar spectrum, including the famous measurements of Kurucz.
          
The code is currently capable of assuming two line profiles when computing absorption/stimulated emission rates: a delta function, and a doppler (thermal) profile. For typical comet distances > 1 AU, the two profiles produce very similar g-factors.

Be advised that the code DOES NOT check that all lines are connected to ground; that is up to the user. As a short-term workaround, it is recommended to set the "upper_cutoff" variable to the ionization potential value in cm^-1 so that auto-ionizing levels are not included. Auto-ionizing level data is typically lacking transition rate data, and will not have a direct or indirect radiative pathway to ground, leading to failures in the matrix solver. An automated ground-connection check routine is in development and will be added in a future update alongside more complete documentation.


-------
[1] Fontenla, J. M., Landi, E., Snow, M., & Woods, T. 2014, Solar Physics, 289, 515, doi: 10.1007/s11207-013-0431-4

[2] Hall, L. A., & Anderson, G. P. 1991, Journal of Geophysical Research: Atmospheres, 96, 12927, doi: https://doi.org/10.1029/91JD01111

[3] Coddington, O. M., Richard, E. C., Harber, D., et al. 2021, Geophysical Research Letters, 48, e2020GL091709, doi: https://doi.org/10.1029/2020GL091709


