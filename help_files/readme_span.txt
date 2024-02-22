SPAN: SPectral ANalysis software V5.0
Daniele Gasparri, February 2024

****Purpose****
SPAN is a Python 3.X graphical interface program designed to perform operations and analyses on astronomical wavelength calibrated 1D spectra.
SPAN has been developed and optimized to analyse galaxy and stellar spectra in the optical and near infra-red (NIR) bands.
The program accepts ASCII and fits binary tables. SPAN deals with linear sampled spectra, with wavelength in physical units (A, nm and mu). If you don't have linear sampled spectra and/or with log/ln wavelength scale, SPAN will try to read the spectra, will convert them automatically to linear sampling and will assign a physical wavelength scale. If these operations fails, your spectra will show strange wavelength scale when clicking 'Plot'. If that is the case, you will need to adjust them with other software before load to SPAN.

The program has been tested with IRAF-reduced spectra, SDSS spectra, IRTF (also extended version) spectra, SAURON spectra, X-Shooter library spectra, (E)MILES, GALAXEV and FSPS stellar libraries, and, in general, with the ESO standard for spectra. The software DOES NOT accept ASCII spectra file with fortran scientific notation, like the PHOENIX synthetic stellar spectra. In this case, you will need to open the file and substitute the scientific notation of flux and wavelength "D" with "E" (you can do this operation even with the embed text editor of SPAN).

Currently, SPAN reads only the wavelength and the flux, discarding the (potential) column with uncertainties.

****What do you need to run SPAN****
- In order to run the source code, you need Python 3.X and the following libraries (some of them are installed by default, others will require: pip install <library>):

    1) Pysimplegui (version >= 4.55. Please, check for updates.)
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
    
    To automatically check and install the missing packages, you can run the script file "install_modules_span.sh" in the terminal:
    ./install_modules_span.sh

    The stand-alone executables applications are ready to use.

 - A screen resolution of at least 1600X900, otherwise the panel will be truncated.
 
    
****How SPAN works****
SPAN can work with just one 1D spectrum, either in Fits or ASCII file, with the first column to be wavelength and the second flux.
SPAN can load and process a list of n 1D spectra, where n must be greater than 1. In order to do this, you need to create and load a txt file containing the relative path of the spectra (with respect to the location of the main SPAN program) and the complete spectra names. The first row of this list file must be commented with # and usually contains something like that: #Spectrum. You can put any type of 1D spectra in this file list, but I strongly suggest to insert spectra with at least the SAME wavelength unit scale.
It seems difficult, but don't worry: the button 'Create spectra list' will help you to create a spectra list file by selecting a folder containing the spectra you want to process.
You can find example file lists in the example_files directory. They are:
1) xshooter_vis_sample_list_spectra.dat , already preloaded in the main application (you just need to click "Load Spectra"), contains 5 spectra of the central regions of nearby galaxies observed with the VIS arm of ESO XShooter spectrograph at resolution of R = 5000. Wavelength units are in nm. Sampling is linear and the wavelength units are 'nm';
2) NGC5320_bins.dat, ngc5806_bins.dat and ic3392_bins.dat contain the spatial bins of three spiral galaxies observed with the TNG telescope at resolution FWHM of 3.5 A from 470 to 670 nm. Sampling is logaritmic and wavelengths are in log(A). SPAN will take care of weverything; you just need to set "A" in the "Open Spectra" frame of the main application before clicking "Load";
3) irtf_K-M_giants_list_spectra.dat contains a sample of 31 giant K and M stellar spectra of the IRTF library from 800 to 5000 nm. Wavelength units are in mu, so you need to check the option "mu" in the Open Spectra" frame befoce clicking "Load spectra".


****Quick start****
If you want to compile the source code, type in the terminal: python3 span_5.0.py, then press the "Load spectra" button to load the example files.
The spectra loaded will appear in the upper central frame (the white window). 
Just select one spectrum with the mouse, then click "Plot" to see the spectrum. Close the plot to activate again the main panel. 
You can select one of the many operations possible in the other frames, for example the "Add noise", then press "Preview spectrum" to see the result. If you like it, you can click to "Process selected" button to save the work (but first you need to close the plot window!). If you press "Process all" you will apply the operation selected to all the loaded spectra. The terminal is your friend: it will tell you all the things that the program is doing.

Now, let's try something in the "Spectral analysis" frame. Let's activate the "EW measurement" task and this time click on the "Preview result analysis" to see the result of the operation. 
The spectrum does look strange to you? Did you unactivate the "Add noise" option above? If not, the "EW measurement" task is analysing this noisy spectrum and not the original one. Congratulations, you just discovered why the "Spectral Analysis" frame is separated from the others: it will consider the spectrum (or the spectra) processed by the other upper frames. If you activate 10 tasks, the spectrum processed in this frame will be the sum of all the activated tasks above. If you don't activate any task, the spectrum processed will be the original loaded. 

Don't be shy: the program has been made to not (almost) crash and if you perform an invalid operation a pop-up window will tell you what's the problem and how to solve it. 


#################################################################################################
################################### General description and usage ###############################

SPAN is composed by a main graphical window that shows all the operations that the user can perform on 1D spectra.
In this window you will find three panels separated by a horizontal line: the top one, the middle and the bottom.
The top and the middle panels are composed by three frames each. 
The program is provided with example spectra and preloaded values that allow the user to run it without any modification, to see how it works.

###The upper panel###
Any operation begins within the upper-left frame (Open Spectra), where you need to insert a valid ASCII file containing the list of spectra you want to process. If you have only one spectrum to analyze, you can upload it by checking the option "I have just one spectrum." By selecting the wavelength unit scale for the spectra (assuming they all have the same wavelength scale) and clicking "Load spectra," the spectra (or the single spectrum) are loaded into the system and displayed as a list in the middle frame. Don't worry: if there is any problem with a spectrum (either it does not exist or is not readable), either on the list or in the "one spectrum" option, the program will alert you with a popup message. If everything goes smoothly, the spectra (or the single spectrum) will be displayed in the middle frame; just select one spectrum with the mouse and perform the operations you want to do.

The right frame (Utilities) is a standalone frame that allows you to find out information about the selected spectrum, such as the header, the sampling (in nm), the Signal-to-Noise Ratio (SNR) on a user-defined window with a user-defined interval (in nm), or simply convert the spectrum to ASCII or binary fits.


###The middle panel###
This section contains some basic tasks that can be performed on spectra, grouped into the 'Spectra pre-processing,' 'Spectra processing,' and 'Spectra math' frames. Any task executed within these frames modifies the selected spectrum and affects the 'Spectral analysis' frame. You can choose multiple tasks (e.g., rebinning, dopcor, adding noise...) without limitations. The 'Preview spec.' button allows you to observe the effect of the task(s) performed. The spectrum displayed and used in the bottom panel ('Spectral analysis') will be the one resulting from the selected operations. No intermediate graphical information is available, so if you perform three tasks, you will see the combined effect of all three. If you don't perform any task, don't worry; the original spectrum will be visible and used for spectral analysis.

The four math tasks in the 'Spectra math' panel that involve all the spectra ('Average all,' 'Normalize and average all,' 'Sum all,' 'Norm. and sum all') act on all the original spectra loaded (and don't work if you have loaded just one spectrum) and remain insensitive to other tasks performed. By activating the 'Use for spec. an.' option, you force the program to utilize the results of these operations for spectral analysis, disregarding any other tasks performed on individual spectra. Be cautious in managing this option. In any case, a message in the terminal window will appear, indicating that you will use the combined spectra for spectral analysis.



###The bottom panel###
This panel is composed of two frames. The left one contains basic and more advanced Spectral analysis tools. The Spectral analysis frame contains the following task: 1) Blackbody fitting, 2) Cross-correlation, 3) Velocity dispersion measurement, 4) Equivalent width line measurement, 5) Line fitting, 6) Kinematics with ppxf, 7) Stellar populations with ppxf, 8) Sigma coefficients determination, and 9) Correction of the EWs for the sigma coefficients. Each task is independent and does not modify the spectrum.
The input spectra can be modified by the operations performed in the upper panel. If you want to perform the analysis on the original spectra, ensure that all previous operations are deactivated (they are all deactivated by default!).
The "Preview" button will display the task(s) result on the selected spectrum in a graphic window and in the output frame on the right, except for the "Sigma coeff determination" and "Correct EWs for sigma" tasks, which do not depend on the spectra loaded in the list. If no task is selected, a warning message will pop up when clicking on the "Preview" button.

The right frame displays the text output of the software. This is how SPAN communicates with you. This panel reproduces the computer terminal and shows the output of the operations performed, including errors and warnings.


###Apply the tasks###
Once you are satisfied with your work, you can process the spectra or the single spectrum. Clicking on the 'Process selected' button will perform all the tasks activated in the middle frame and the "Spectral analysis" frame, saving the new processed spectrum to a text file. By default, the program will save intermediate files, one for each task performed, in the order they appear. For example, if you have selected rebinning, sigma broadening, and add noise, the program will save a spectrum with rebinning done, a second spectrum with rebinning + sigma broadening applied, and a third with rebinning + sigma broadening + add noise applied.
If you are not interested in saving all the intermediate spectra files, you can deselect the 'Save intermediate files' option, and only the spectrum at the end of the selected operations will be saved.
The results of the spectral analysis frame will be written in the output frame.

By clicking "Process all," you will apply all the tasks to all the spectra in your list. This is the only way to save the results of the Spectral analysis frame in a text ACII file.


### The sub windows###
In the bottom part of the GUI, along with 'Process selected' and 'Process all' buttons, you will find 4 light blue buttons. These are sub-program utilities that will allow you to: open, create and modify an ASCII file (Text editor), add and modify the keywords in the header of fits files (FITS header editor), plot the data generated by SPAN and, in general, any data stored in an ASCII text file (Plot data), and extract 1D spectra from a 2D reduced fits image (2D spec extraction).


********************************************************************************************************
*************************************** The input files ************************************************

In order to work properly, the program needs some text files containing some informations about your data. To see how they must be formatted, please take a look at those coming with SPAN and already set by default in the graphic interface.

IMPORTANT: 1)The text files MUST always have the first line as header, identified by # (e.g. #spectrum)
           2)The spectra list files MUST always have the path to their location, if different than the location of SPAN. If it is in a subfolder you can put the local path. If it is in a different directory, you need to put the absolute path.
           
1) Spectra file list task:
    It is essential. If you don't believe it, try to perform any task without upload the spectra and you will see the effects! It is just an ASCII file containing the names of the spectra you want to process and [optionally] the path to find them. You can use any spectra you want, with different format (fits, ASCII...) and resolutions, but I strongly suggest to use spectra with the same wavelength units and, more important, that cover roughly the same wavelength interval. Don't mix up, for example, spectra in the visible and NIR domain. This might cause the program to crash when the 'Process all' is performed, but you should still be able to process and preview the tasks performed on the single spectra. If you just want to play with one spectrum, then load the ASCII or fits 1D spectrum and activate the option "I have just one spectrum" before clicking the button "Load spectra".

                                            example_list.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]spectrum1.fits
                                            [path/]spectrum2.fits
                                            [path/]spectrum3.fits
                                            [path/]spectrum4.fits

                                            
2) Doppler correction file for the Dopcor task and the 'I have file' option selected:
    It's an ASCII file containing two columns: 1) Name of the spectrum and 2) Radial velocity to correct the spectrum. This file has the same format of the output text file generated by the Cross-correlation task, so you can directly use it. 
                    
                                            example_dopcor.dat
                                        
                                        #spectrum       RV(km/s) ---> header: always necessary!
                                        [path/]spectrum1.fits  1000
                                        [path/]spectrum2.fits  1001
                                        [path/]spectrum3.fits  1002
                                        [path/]spectrum4.fits  1003
                                            
3) Heliocentric correction file for the Helio correction task and the 'I have file' option selected: 
    It's an ASCII file containing three columns, separated by a space: 1) Name of the location, 2) Date of the observation (just year, month, day, not the hour), 3) RA of the object (format: degree.decimal), 4) Dec of the object (format: degree.decimal). IMPORTANT: the file must be in the same directory of SPAN.
    
                                            example_heliocorr.dat
                                        
                                #where  date        RA          Dec
                                paranal 2016-6-4    4.88375     35.0436389
                                paranal 2016-6-30   10.555      1.11121
                                aao     2011-12-24  -50.034     55.3232
                                aao     2018-2-13   -11.443     11.2323
                                SRT     2020-7-31   70.234      55.32432

                            
4) Cross-correlation and velocity dispersion tasks:
    These task require a template, in fits or ASCII format (ie. just a spectrum!)
    
    
5) Equivalent width measurement task and 'I have an index list file' option selected:
    It's an ASCII text file containing the index definitions. One index per column. Don't mess it uo with the index file, otherwise you will obtain inconstistent results! Luckily, you can always test a single index and see the graphical preview before running the wrong indices on 240913352 spectra and waste one year of your life.
    
                                            example_idx_list_file.dat

                                    #Idx1    Idx2  ---> header: always necessary!
                                    847.4   847.4 ---> row2: left blue continuum band, in nm
                                    848.4   848.4 ---> row3: right blue continuum band, in nm
                                    856.3   856.3 ---> row4: left red continuum band, in nm
                                    857.7   857.7 ---> row5: right red continuum band, in nm
                                    846.1   848.4 ---> row6: left line limits, in nm
                                    847.4   851.3 ---> row7: right line limits, in nm


6) Sigma coeff determination task: 
    It determines 4 spline correction coefficients in order to correct the EW of galactic spectra modified by the broadening due to the velocity dispersion. It needs a sample of unbroadened spectra that are a good match of the expected stellar populations of the galactic spectra you want to correct to the zero velocity dispersion level. The input file is just an ASCII file containing the list of the spectra used as sample. by default, the program contains the spectra of 31 giant K and early M (<5) stars of the IRTF catalog. This means that this sampe is suitable only for the NIR band (850-2400 nm). 

                                            example_coeff_determ.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]stellar_spectrum1.fits
                                            [path/]stellar_spectrum2.fits
                                            [path/]stellar_spectrum3.fits
                                            [path/]stellar_spectrum4.fits

                                    
7) Correct the EWs for the sigma task: 
    To apply the velocity dispersion coefficients for correcting the EW to a zero velocity dispersion value you'll need this task and three files: 
        1) Sigma list file: a file containing the name of the spectra, the velocity dispersion and their error. It has the same format of the output file generated by the Velocity dispersion measurement task. 
        
                                            example_sigma_vel.dat
                                            
                                            #Spectrum       Sigma(km/s) err ---> Header: always necessary
                                            spectrum_name1  166.2       3.0
                                            spectrum_name2  241.5       3.1
                                            spectrum_name3  335.1       6.2
                                            spectrum_name4  241.5       3.2
        
        2) EW file list to correct: the text file containing the EW values you want to correct. It has the same format of the output file generated by the EW measurement task. BE CAREFULL to check that the indices are in the EXACT same order of those you used in the Sigma coeff determination task for the correction coefficient determination.
        
                                            
                                            example_uncorrected_ew.dat
                                            
                            #Spectrum       idx1    idx2    idx3   idx1err idx2err idx3err
                            spectrum_name1  0.27    1.38    3.56    0.01     0.01    0.02
                            spectrum_name2  0.15    1.32    3.43    0.01     0.02    0.02
                            spectrum_name3  0.08    0.75    2.81    0.01     0.02    0.02
                            spectrum_name4  0.14    1.25    3.18    0.01     0.01    0.01

        
        3) Correction coefficients file: it's the output file generated by the Sigma coeff determination task. 
                                    
                                            example_correction_coeff.dat
                        Pa1          Ca1           Ca2          Pa1e         Ca1e         Ca2e
                        4.3282e-08   1.06712e-08  -2.7344e-09  -5.7463e-09   2.2911e-09   2.8072e-10
                       -2.9602e-05  -1.2012e-05   -3.5782e-07   3.9353e-06  -1.9246e-06  -2.9293e-07
                        0.0017       0.0021        8.5793e-05  -0.0001       0.0004       9.9212e-05
                       -0.0029      -0.0085       -0.0016       0.0053      -0.0003      -0.0002

                                    


In general, be always sure that the wavelength interval you put as parameter in some tasks, matches the wavelength interval of the spectrum and the wavelength units. If not, the result will be a general crash. I will work to fix that and give just a mild warning from the program instead a catastrophic crash, but I'm sure that a crash is more effective to let you check all the wavelength intervals!


##################################################################################################
################################## List of operations you can perform ############################
SPAN can perform many operations on the spectra and

WARNING: Except where clearly stated by the program, all the wavelengths of SPAN are given in nm, in air, and all the velocities in km/s.

Here is a description of the functions:

1) Utilities frame: stand alone frame with action buttons on the right.
    a) HDR = shows the header of the selected spectrum, both fits and ASCII;
    b) Step = shows the step of the selected spectrum;
    c) Res. = shows the resolution of the selected spectrum by trying to fit an emission sky line. In The W1 and W2 you should put a small wavelength interval containing a sky line: it's up to you!
    d) Convert spectrum to = converts the selected spectrum to ASCII of Fits;
    e) Compare with = Compares the selected spectrum with another one selected by the user;
    f) Convert flux = converts the flux from frequency to lambda and vice-versa. The buttons "see plot", "save one" and "save all" are active to see and save the results for one or all the spectra;
    g) SNR = measures the Signal to Noise in the selected spectrum, in the W. central wavelength selected by the user. The buttons "save one" and "save all" are active to save one or all the SNR calculated for the spectra.

2) Spectra pre-processing frame
    a) Cropping = performs a simple cropping of the spectra. If the wavelength window to crop is outside the spectrum, SPAN will ignore the task and will not perform the crop;
    b) Dynamic cleaning = performs a sigma clipping on the spectra. The sigma clip factor, the resolving power of the spectrum and the velocity dispersion (instrumental and/or intrinsic) of the selected spectrum is required in order to perform a better cleaning. For the "Process all" the option "R and sigma vel file" is available in order to have R and sigma values for all the spectra to be processed. Be careful to use it with undersampled spectra;
    c) Wavelet cleaning = performs a wavelet based denoise of the spectra. The mean standard deviation of the spectra continuum (sigma) and the number of wavelet layers to consider are required. You don't need to measure it, just try different values. Be careful to not delete the signal;
    d) Filtering and denoising = smooths the spectra by performing some denoising filter: box window moving average, gaussian kernel moving average, low-pass Butterworth filter and band-pass Butterworth filter;
    e) Dopcor = performs the doppler correction of the spectrum. Single shot option with user input value of recession velocity (in km/s) is available both for one or all the spectra. "I have a file" option" only works with the "Process all": it contains a text file with the spectra name and the recession velocities. This files is generated by the "Cross-correlation" task in "Process all" mode;
    f) Helio corr = performs the heliocentric correction on the spectra. The "Single" option require a location, that can be selected from the "loc.list" button (it requires an internet connection). The other fields are the date in the format YYYY-MM-DD and the RA and Dec. of the observed object (in decimals). In the "I have a file" option, available only for "Process all" a list file with location, date, RA and Dec. coordinates for each object is provided.

3) Spectra processing frame
    a) Rebin = performs a rebin/resample of the spectra in linear wavelength step ("pix.lin" option, with the step in nm) and in sigma linear step ("sigma lin." option, with the sigma step in km/s);
    b) DegradeR = degrades the resolution of the spectra in terms of resolving power R. The original and the final R must be provided;
    c) DegradeL = degrades the resolution of the spectra in terms of delta lambda (useful, for example, to match the resolution of the Lick/IDS indices);
    d) Normalize to = normalises the spectra to the wavelength provided by the user (in nm);
    e) Continuum model = models the continuum shape with two options: 1) Simple filtering of the continuum by reducing the spectrum to a very small resolution (R = 50), and 2) polynomial fitting, with the possibility to masks emission/contaminated regions. Both the continuum models can be divided or subtracted to the original spectrum;
    f) Sigma broadening = broads the spectra by convolving with a gaussian function with the standard deviation provided by the user, in km/s. Remember that the real broadening of the spectra will be the quadrature sum between the broadening and the instrumental sigma of the spectra;
    g) Add noise = adds a random poisson noise to the spectra with a SNR defined by the user. Remember that the final SNR of the spectra will the the sum in quadrature between the added noise and the intrinsic SNR of the spectra;


4) Spectra math 
    a) Average all = averages all the spectra (only available in "Process selected" mode);
    b) Norm. and average all = normalizes to a common wavelength and average all the spectra (only available in "Process selected" mode);
    c) Sum all = sums all the spectra (only available in "Process selected");
    d) Norm. and sum all = Normalizes and sum all the spectra (only available in "Process selected"). The option "Use for spec. an." forces the program to use the result of one of these 4 operations for the following spectral analysis;
    e) Subtract normalized average = subtracts to the spectra the normalized average made from all the spectra loaded;
    f) Subtract norm. spec. = subtracts to the spectra a normalized spectrum selected by the user;
    g) Add pedestal = add a constant to the spectra;
    i) Multiply by a constant = multiplies the spectra by a user defined constant value.
    
5) Spectral analysis. This is the core of the program. Each task can be fine-tuned by clicking on the relative parameters button on the right:
    a) Blackbody fitting = performs a fit of the spectrum with Plack's blackbody equation and gives the temperature estimation. It works with any type of spectra but it performs better for stellar spectra, with wide (at least 500 nm) wavelength range.
    b)  Cross-correlation = performs a cross-correlation of the spectra with a user selected template. This uses the function crosscorrRV of pyastronomy. The user can smooth the template to a velocity dispersion value in order to improve the cross-correlation and should identify a narrow region of the spectrum to be cross-correlated (tip: the Calcium triplet lines are the best features);
    c) Velocity dispersion measurement = performs the measurement of the velocity dispersion of the spectra by fitting with a user provided template. Some pre-loaded bands in the visible and NIR are shown but the user can select an independent band. The routine succeeds with strong features (the CaT is the best);
    d) Equivalent width = performs the equivalent width measurement of the spectra, with a single user provided index, with a list of indices or the Lick/IDS system. The results are provided in Angstrom. Montecarlo simulations are run for the uncertainties estimation. The calculation of the Lick/IDS indices can be personalized in many ways: you can correct for the emission, for the velocity dispersion and the recession velocity. You can also perform the linear interpolation with the SSP models of Thomas et al. 2010 to retrieve the age, metallicity and alpha enhancement of the stellar populations;
    e) Line(s) fitting = performs the fitting of an user defined line with user defined parameters and a combination of gaussian model with straight line. If "CaT lines" is selected, the task will perform an automatic fitting of the Calcium Triplet lines, assuming they have been corrected to the rest frame velocity;
    f) Kinematics with ppxf = uses the known ppxf algorith of Cappellari to fit a user defined wavelength region of the spectra with a combination of templates. You can select the template library you want between the EMILES, GALAXEV e FSPS. It returns the radial velocity, the velocity dispersion and the higher moments H3 and H4 (and a nice plot courtesy of Cappellari);
    g) Stellar populations with ppxf = uses the known ppxf algorith of Cappellari to fit a user defined wavelength region of the spectra with a combination of templates. You can select the template library you want between the EMILES, GALAXEV e FSPS. The user can decide whether include the gas emission or not, the reddening and the order of multiplicative and additive polynomials of the fit. It returns a beautiful plot, the kinematics, the weighted age (in luminosity and mass), metallicity (weighted in luminosity and mass), the M/L and save the best fit template and the emission corrected spectra (if any). Works great in the visible and in the NIR, but that depends on the quality of the template library used;
    h) Sigma coeff determination = standalone procedure to determine the coefficients necessary to correct the EW to zero velocity dispersion value. It broadens a sample of K and early M stars of the IRTF library up to 400 km/s and calculates the deviation of the equivalent width of the index/index file provided in the EW measurement task. It works only by pressing the "Compute!" button and in "Preview result analysis" mode and creates a text file with a third order polynomial curve that fits the behaviour of any broadened index;
    i) Correct EWs for the sigma = performs the correction of the EW based on the coefficients estimated with the Sigma coeff determination task. It works only by pressing the "Correct!" button mode and require an EW measurement files with the same indices, in the same order, to that considered in the Sigma coeff determination. The output files of the "EW measurement", "Sigma coeff determination" and "Velocity dispersion measurement" are ready to be used for this task, if we are considering the same spectra and indices.
    

##################################################################################################
################################## The sub-programs ############################

The four light-blue buttons at the bottom of the GUI are sub-programs that might help you in the difficult task of analysing and processing astronomical spectra. They work independently from the main program, so you can also not load spectra if you don't need to perform tasks on them. Here is how they works:

1) Text editor: a simple ASCII file editor where you can create, read or modify ASCII files, included those generated by the SPAN tasks. Some basics operations are available, such find, replace and merge rows;

2) FITS header editor: an header editor to add, remove and save the keywords of fits header files. You can choice between: 'Single fits header editor' to work with the keywords of one fits file, 'List of fits header editor' to modify the keywords of a list of fits files, 'Extract keyword from list' to extract and save in an ASCII file one or more keywords from the headers of a list of fits files;

3) Plot data: a sub-program to plot the data generated by the Spectral analysis task and, in general, all the data stored in ASCII space separated data. Once you browse for the text file and click the 'Load' button, the program will automatically recognise the column names. Select a name for the x and y axis and plot the data to see them in an IDL style plot.
You can personalise the plot by adding the error bars, set the log scale, add a linear fit (simple fit without considering the uncertainties), set the labels, the range, the size of characters, size and colours of the markers and decide if visualise the legend or not. You can also save the plot in .eps format, in the same main directory where resides SPAN.
If any error occur, the program will warn you;

4) 2D spec extraction: allows the extraction of a single 1D spectrum or a series of 1D spectra from a reduced and wavelength calibrated 2D fits image containing the long-slit spectrum of a source, with dispersion axis along the X-axis and the spatial axis along the Y-axis.
Before proceed to the extraction, you need to load a valid 2D fits image, then you need to:
    1) Open the spectrum and see if everything is ok;
    2) Fit the trace in order to find the maximum along the dispersion axis. You need to set the degree of polynomial curve that will be used to fit the trace and correct the distortion and slope of the spectrum;
    3) Correct the spectrum for distortion and slope using the model trace obtained in the previous step.
Then, you can choose:
    a) extract and save only one 1D spectrum within the selected Y range (useful for point sources);
    b) extract and save a series of n 1D spectrum covering all the spatial axis and obtained by binning contiguous rows in order to reach the desired Signal to Noise Ratio (SNR).
    The SNR threshold that you must insert is just a rough estimation of the real SNR. A good starting value to produce 1D spectra with bins with realistic SNR > 30 is 35. Adjust the SNR Threshold to your preference by looking at the real SNR of the bins.
    The pixel scale parameter is optional. If you set to zero it will not be considered. This option is useful if you have the spectrum of an extended source (e.g. a galaxy) and want to sample different regions.


##################################################################################################
################################## Utilities in the menu bar ############################


The menu bar was introduced in version 4.5 of SPAN, offering several helpful options to enhance your experience with spectral analysis. Here's a detailed overview of some options that you won't find in the main panel:

1) File --> Save Parameters...: Allows you to save all parameters and values from the main panel in a .json file.
This feature is invaluable as it enables you to preserve any modifications made to parameters, facilitating the restoration of your session each time you reopen SPAN. Multiple configuration files can be saved and loaded into the program;
2) File --> Load Parameters...: Permits the loading of parameters saved in a .json file. This functionality allows you to resume your work with personalized parameters instead of modifying default ones every time. Default parameters are loaded each time the program starts, making it convenient to use saved configurations;
3) File --> Restore Default Parameters: Provides the ability to reset all parameters to their default values. Useful if numerous parameter modifications during a lengthy session have resulted in issues, allowing you to start fresh;
4) Edit --> Clear All Tasks: Immediately deactivates all tasks activated during the session, enabling a clean restart;
5) Edit --> Clean Output: Deletes all content in the output window. Particularly useful during extended sessions where the generated output may become extensive.


Please, report any bug or comment to daniele.gasparri@gmail.com
Have fun!

Daniele Gasparri
01/02/2024
Greetings from the Atacama desert!
