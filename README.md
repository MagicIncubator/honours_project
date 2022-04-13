# honours_project
The code from Jordan Pagé's honours project.


The following files were written between September and Decemeber 2021:

  deatailedBalance_V.py: 
This file finds the optimal mu_CV for a spceified Eg and plots the 
efficiency and the JV curve vs mu_CV for the specified eG.

  deatailedBalance_eG&V_plotV.py: 
this file finds the optimal mu_CV and Eg for given starting parameters 
and plots the JV and efficiency curves vs mu_CV.

  deatailedBalance_eG&V_plot_eG.py: 
this file finds the optimal mu_CV and Eg for given starting parameters 
and plots the efficiency curve vs Eg.


The following files were written between January and April 2022:

  detailed_balance_IBSC.py: 
This file prompts the user wether to find the efficiency limit of an 
IBSC device (Luque & Marti 1997), recreate figure 2 from the Luque & 
Marti 1997 paper, or do both.

  find_alpha_eff.py: 
this file imports the solar spectrum data and the data on the index 
of refraction of AlGaAs (all of which is located in folder
'absorption&spectrum_data'. There are other options of materials to use in 
the folder.). This file uses the complex part k of the index of refraction 
to calculate the absorption coefficient vs wavelength of the material. It 
then calculates an effective absorption coefficient for the portion of the 
solar spectrum above the minimum energy needed for a transition in that
material (this energy level is given by the user).

  plot_band_diagram.py: 
this file imports data output from a Simudo simulation and plots the band
digram of the device.

  plot_JV.py: 
this file imports data output from a Simudo simulation and plots the 
JV-curve of the device.

  example_output
This folder contains an example of what is output from Simudo. The data in the folder is called by plot_JV.py and plt_band_diagram.py.
