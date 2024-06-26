
This github repo contains versions of the python code 'FlorPy' as described in

Bromley et al PSJ 2021 (Application: atomic Ni and atomic Fe)
Bromley et al MNRAS 2023 (Application: CO+)
Bromley et al Icarus 2024 (Application: CS)


The older versions used for the 2021 and 2023 manuscripts have been merged into a newer updated version, stored in the folder "current_version".
The new version contains numerous bug fixes, speed ups, and also the addition of new utility codes.

The repo is organized as follows:

The main model codes are stored in florpy_dict_v*.py and molecular_utils_florpy.py under /current_version/

The three provided test cases (atomic Ni, CO+, CS) are in the form of .py files, each in their own directory. Running them will provide multiple outputs in the same directory. These examples show how to run some of the basic functionality, such as:

-interacting with the dictionary structure
-plotting gfactors or band luminosities
-running the time-dependent models
-using the MonteCarlo error estimator
-using the thermalized ground state implementation (see CS examples).

The input files for these models, in the form of transition lists and energy level lists, are in formats similar to those provided for data downloaded from the NIST ASD. They are available in /input_files/, along
with some additional information about how they were generated.

For ATOMIC systems, that

In the meantime, if you have concerns, questions, or suggestions for improvements, changes, or added functionality, please contact me at sjb0068@auburn.edu. I am happy to coordinate any of these.

SJB

-------

The code is built to use atomic data from the NIST ASD:
https://physics.nist.gov/PhysRefData/ASD/lines_form.html

For ATOMIC systems, the user must acquire lines and levels files by querying the NIST ASD for lines (only those with transition rates), and energies in units of cm^-1:

-- On the lines page, enter the element of interest, e.g. "Ni I", and under "Show Advanced Settings" select "Only lines with transition probabilities. Under "Format Output" select Tab-Delimited. Once the page is loaded, right click and "Save As" with your preferred naming convention. Repeat the same for the levels (Tab-Delimited -> Save As). These two files will be required inputs for the python script that runs the model. 

A sample file, "ni0 example script.py" shows the syntax and order of function calls. This sample script will general g-factors for neutral nickel at 1 AU (heliocentric velocity = 0 km/s).

The syntax for performing Monte-Carlo iterations to generate approximate model uncertainties with 10 iterations is also provided. Note that 10^4 or more iterations are required for good convergence. 10^4 iterations will execute in around 15 - 20 minutes on compute time on a modern macbook. Be advised that the RAM required to store more runs (e.g. 10^6) at once may not be possible on a personal computer.

-------- 