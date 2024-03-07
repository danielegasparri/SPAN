# SPAN
SPAN: SPectral ANalysis software V5.0
Daniele Gasparri, February 2024



****Purpose****

SPAN is a Python 3.X graphical interface program designed to perform operations and analyses on astronomical wavelength-calibrated 1D spectra, optimised for the visible and near-infrared bands and galaxy spectra. 

The program accepts ASCII and fits binary tables 1D spectra.
SPAN deals with linear sampled spectra, with wavelength in physical units (A, nm and mu). If you don't have linear sampled spectra and/or with log/ln wavelength scale, SPAN will try to read the spectra, will convert them automatically to linear sampling and will assign a physical wavelength scale. If these operations fails, your spectra will show strange wavelength scale when clicking 'Plot'. If that is the case, you will need to adjust them with other software before load to SPAN.

The program has been tested with IRAF-reduced spectra, SDSS spectra, IRTF (also extended version) spectra, SAURON spectra, X-Shooter library spectra, (E)MILES, GALAXEV and FSPS stellar libraries, and, in general, with the ESO standard for 1D spectra.
The software DOES NOT accept ASCII spectra files with fortran scientific notation, like the PHOENIX synthetic stellar spectra. In this case, you will need to open the file and substitute the scientific notation of flux and wavelength "D" with "E".

Currently, SPAN reads only the wavelength and the flux, discarding the (potential) column with uncertainties.



****What do you need to run SPAN****

- In order to run the source code, you need Python 3.X and the following libraries (some of them are installed by default, others will require: pip install <library>):
    1) Pysimplegui (version >= 4.55. Please, check for updates)
    2) Astropy
    3) Pandas
    4) Numpy
    6) Matplotlib
    7) Scipy
    8) ppxf
    9) scikit-image
    10) PyWavelets
    11) Os (already in your python)
    12) Time (already in your python)
    13) Math (already in your python)
    14) Io (already in your python)

 - A screen resolution of at least 1600X900 pix is required, otherwise the graphical panel will be truncated. If your screen vertical resolution is lower than 900 pix, try the span_5.0_low_res_screen.py version, which works with vertical resolutions as low as 768 pix.
 

    
****How SPAN works****

SPAN can work with just one 1D spectrum (option "I have just one spectrum" in the upper left frame of SPAN), either in fits or ASCII file, with the first column to be wavelength and the second flux.
SPAN can load and process also a list of n 1D spectra, where n must be greater than 1. In order to do this, you need to create and load an ASCII file containing the relative path of the spectra (with respect to the location of the main SPAN program) or the absolute path, with the complete spectra names. The first row of this list file must be commented with # and usually contains something like that: #Spectrum. You can put any type of 1D spectra in this file list, but I strongly suggest to insert spectra with at least the SAME wavelength unit scale.
It seems difficult, but don't worry: the button "Create spectra list" in the main panel of SPAN will help you to create a ready to use spectra list file by just selecting a folder containing the spectra you want to process.

In this repository you can find example file lists in the example_files directory. They are:
1) xshooter_vis_sample_list_spectra.dat, already preloaded in the main application (you just need to click "Load Spectra"), contains 5 spectra of the central regions of nearby galaxies observed with the VIS arm of ESO XShooter spectrograph at resolution of R = 5000. Wavelength units are in nm. Sampling is linear;
2) NGC5320_bins.dat, ngc5806_bins.dat and ic3392_bins.dat contain the spatial bins of three spiral galaxies observed with the TNG telescope at resolution FWHM = 3.5 A, from 470 to 670 nm. Sampling is logarithmic and wavelengths are in log(A). SPAN will take care of everything: you just need to set "A" in the "Open Spectra" frame before clicking "Load spectra";
3) irtf_K-M_giants_list_spectra.dat contains a sample of 31 giants K and M stellar spectra of the IRTF library from 800 to 5000 nm. Sampling is logarithmic. Wavelength units are in 'mu', so you need to check the option "mu" in the "Open Spectra" frame befoce clicking "Load spectra".



****Quick start****

If you want to compile the source code, install all the required modules and then type in the terminal: python3 span_5.0.py, then press the "Load spectra" button to load the example files. 
The span_5.0_win.py file is the same program optimized (I hope) for Windows operating systems, with just some minor graphical adjustments. Feel free to use the version that better adapts to your screen and machine settings.

The spectra loaded will appear in the upper central frame (the white window). 
Just select one spectrum with the mouse, then click "Plot" to see it. Close the plot to activate again the main panel. 
You can select one of the many operations possible in the other frames, for example the "Add noise", then press "Preview spec." to see the result. If you like it, you can click "Process selected" to save the work (but first you need to close the plot window!). If you press "Process all", you will apply the operation selected to all the loaded spectra. The results will be stored in the folder "results/". The output window is your friend: it will tell you all the things that the program is doing.

Now, let's try something in the "Spectral analysis" frame. We activate the "Equivalent width (EW)" task, we click on "Equivalent width parameters" and we select the option "Single index". We confirm the parameters by clicking the button "Confirm" and we preview the results by clicking the "Preview result" button. 
The spectrum does look strange to you? Did you unactivate the "Add noise" task above? If not, the "Equivalent width (EW)" task is analysing this noisy spectrum and not the original one. Congratulations, you just discovered why the "Spectral Analysis" frame is separated from the others: it will consider the spectrum (or the spectra) processed by the other upper frames. If you activate 10 tasks, the spectrum processed in this frame will be the combination of all the activated tasks above. If you don't activate any task, the spectrum processed will be the original loaded. 

If you activated so many tasks that the entropy of the program tends to infinite, don't worry. Just click on the menu: "Edit --> Clear all tasks" to start from fresh. If you want to restore also the default parameters, you can do it with: "File --> Restore default parameters". If you want to save your personal parameters, you can do with: "File --> Save parameters" and load them again whenever you want.

I worked to SPAN during my PhD and my first post-doc fellowship, starting from zero with Python coding and with some Fortran previous experience, so it is very likely that my coding is not pythonic at all: sorry for that. 

Feel free to test and contribute to SPAN!
