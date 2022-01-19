LUMERICAL CODES
1) convert9array.lsf will create power monitor boxes, all other monitors, and adjust accordingly. It will also duplicate 1 device into a 9-device array
2) varyBR_gauss - v5.lsf will save out the .mat's accordingly. Any changes here propagate out to the rest

MATLAB CODES and PROCESSING WORKFLOW
0) compileFN.m contains the sweep parameters that will be used for all below files unless you specifically edit them in their respective files.
1) Run compress.m for any .mat's that are saved out from Lumerical in order to get the filesize right
2) Run plot_sort_spec.m to plot sorting efficiency spectra for every single .mat savefile
2.5) Run plot_Enorm.m to plot scattering patterns at the focal plane for each .mat savefile
3) Run plot_ang_range.m to sweep through whichever variable is being swept in this instance, and get how the efficiency varies as this variable is swept, and then compile it into one whole file for side-by-side comparison plotting across multiple sweeps of a 2nd parameter
4) Run plot_angRange_suprimp.m to load the composed file and plot the side-by-side

MAKE A FOR LOOP FOR ALL THE QUADRANTS

v9: No optimizations were done, simply sweeps of the sidewall material and thickness to understand the ideal range of parameters for which to run simulations.