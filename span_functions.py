#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#
#This file contains all the functions in order to run SPAN

"""
    Copyright (C) 2018-2024, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyze 1D long-slit reduced spectra.
    It is provided "as is" without any warranty whatsoever.
    Permission to use for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided that this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.

"""

######################################################################################
#Version histories

#22/01/2024:
#Updated the functions for kinematics and stellar populations to use the new release of ppxf (9.1.1).
#Inserted consistency checks for the templates' wavelength range.

#14/01/2024:
#Modified the functions involving the EW measurement to align with what spans 4.6 does with the Lick/IDS indices measurements.
#Fixed a bug in the error determination of stellar populations with ppxf.

#21/12/2023:
#Added the Lick/IDS index determination and correction for the sigma.

#07/12/2023:
#Added functions for the plotting window, fits header modifications, and matching rows from the text editor.

#02/12/2023:
#Added the blackbody fitting functions.

#12/11/2023:
#Added Monte Carlo simulation to measure errors in ppxf stellar populations.

#08/11/2023:
#Changed the Pandas .at method to .iloc to assign flux values to be averaged and summed.

#03/11/2023:
#In the read_spec function, added the possibility to correctly load and read 1dfits IRAF style spectra with log lambda.

#20/11/2022:
#Added ppxf routines for stellar kinematics and populations.
#Incorporated the crosscorr function of pyasl and the modified version of the emission line templates of ppxf.

#16/11/2022:
#Now the read_spec can also read SDSS spectra-type.

#14/11/2022:
#Fixed some bugs.

#09/22/2022:
#Modified the Cross-correlation function. Now it works.

#28/01/2022:
#Added the calculation of the EW in magnitudes.

#19/01/2022:
#Added the SNR calculation in the EW measurement function.

#22/11/2021:
#Added the degrade function in lambda (FWHM) to degrade the spectrum to a certain FWHM in Angstrom.

#27/10/2021:
#Fixed the error calculation from the sigma coeff correction.

#19/10/2021:
#Fixed a bug in the normalize and average function.

#07/11/2020:
#Added resolution function to calculate the resolution of the spectrum considering the sky emission lines (the user must provide the wavelength interval of a well-visible sky line to measure its FWHM).

#15/06/2020:
#Changed the read_spec function for the 1dFit spec, excluding PyAstronomy.pyasl since it is not working with the standalone version of the program.
####################################################################################


import numpy as np
import PySimpleGUI as sg   
import math as mt
from scipy import interpolate
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from scipy.optimize import curve_fit
import scipy as scipy
from scipy import ndimage, misc # for smooting like IDL smooth function
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, filtfilt

from scipy.interpolate import interp1d
from scipy.constants import h,k,c
from scipy.integrate import quad
from scipy.optimize import leastsq
from scipy.stats import pearsonr

from matplotlib.backends.backend_ps import FigureCanvasPS
import time
import scipy.stats

#for ppxf
import glob
from os import path
from time import perf_counter as clock

from scipy.signal import correlate2d
from matplotlib.colors import LogNorm
import os

from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from urllib import request

from skimage.restoration import denoise_wavelet
#*****************************************************************************************************
# 1) READ THE SPECTRA: INPUT: SPECTRUM NAME, LAMBDA UNITS. OUTPUT: WAVELENGTH AND FLUX
# Reads ASCII, 1D (flux and delta lambda) and 2D fits spectra (flux, lambda), automatically recognizing them.
def read_spec(spec_name, lambda_units):
    fmt_spec1 = '.txt' in spec_name or '.dat' in spec_name
    fmt_spec2 = '.fits' in spec_name

    # If I have an ASCII spectrum with lambda
    if (fmt_spec1 == True and fmt_spec2 == False or (fmt_spec1 == True and fmt_spec2 == True)):
        spec_type = '2d ASCII table'
        wavelength, flux = np.loadtxt(spec_name, usecols = (0,1)).T
        start_lambda = wavelength[0]
        if (start_lambda < 12. and start_lambda > 5 and lambda_units != 'mu'):
            print ('I think you have ln lambda, try to convert to lambda...')
            wavelength_log = wavelength
            #ln_wavelength= wavelength_log*np.log(10)         # Convert lg --> ln
            wavelength = np.exp(wavelength_log)
        print(spec_type, 'spec with lambda in', lambda_units)
        obj_name = spec_name
        
    # if I have fits files, they can be of different type
    else:
        hdu = fits.open(spec_name)
        hdr_fits = hdu[0].header
        oned_key = 'CDELT1' in hdr_fits # ho assunto che CELT1 sia la keyword che identifica UNICAMENTE i fits 1D: è corretto???????

        #if fits table are 1d (IRAF style, with flux and delta lamba)
        if (oned_key == True):
            spec_type = '1d fits table IRAF style'
            print (spec_type, 'spec with lambda in', lambda_units)
            points = hdr_fits['NAXIS1']
            start_lambda = hdr_fits['CRVAL1']
            step=hdr_fits['CDELT1']
            wavelength = np.arange(points)*step+start_lambda
            flux_tmp = hdu[0].data
            flux = np.array(flux_tmp)

            #reading 1dfits IRAF style with logaritmic wavelength
            if (start_lambda < 5. and start_lambda > 2.5 and lambda_units != 'mu'):
                print ('I think you have log lambda, try to convert to lambda...')
                wavelength_log = wavelength
                ln_wavelength= wavelength_log*np.log(10)         # Convert lg --> ln
                wavelength = np.exp(ln_wavelength)
            if (start_lambda < 12. and start_lambda > 5 and lambda_units != 'mu'):
                print ('I think you have ln lambda, try to convert to lambda...')
                wavelength_log = wavelength
                wavelength = np.exp(wavelength_log)

            star_name = 'HNAME' in hdr_fits
            if (star_name == True):
                star = ''
                star = hdr_fits['HNAME']
                obj_name = star.replace(" ", "")
            else:
                obj_name = spec_name

        # if the fits table are 2d:
        elif (oned_key == False): 
            # Define the columns
            flux_tmp = 0
            spec_type = '2d fits table'
            print (spec_type, 'spec with lambda in', lambda_units)
            

            # trying a connon sense fits table with wavelength and flux, like the ESO does
            try:
                flux_tmp = hdu[1].data['FLUX']
                waves_tmp = hdu[1].data['WAVE']
                eso_spec = True
            except KeyError:
                eso_spec = False #bad_spec = 1
            
            if eso_spec == True:
                wavelength = np.array(waves_tmp)
                flux = np.array(flux_tmp)
                flux_tmp = 0
                # per riconoscere il nome dell'oggetto associato allo spetto. Serve??
                star_name = 'HNAME' in hdr_fits
                if (star_name == True):
                    star = ''
                    star = hdr_fits['HNAME']
                    obj_name = star.replace(" ", "")

                elif (star_name == False):
                    obj_name = spec_name
            
            #americans are different: I try to see if the spectra are in the SDSS weird format, where I have flux and loglam instead of flux and wave
            elif(eso_spec == False):
                
                try:
                    t = hdu['COADD'].data
                    flux_tmp = t['flux']
                    wavelength_tmp = t['loglam']
                    sdss_new_spec = True
                except KeyError:
                    sdss_new_spec = False

                if (sdss_new_spec == True):
                    #the new sdss spectra have the log_wave instead of wave!
                    print ('with log lambda')

                    wavelength_log = np.array(wavelength_tmp)
                    flux = np.array(flux_tmp)
                    #since the lambda are in log, I transform to real numbers, maintaining the log spaced values
                    ln_wavelength= wavelength_log*np.log(10)         # Convert lg --> ln
                    wavelength = np.exp(ln_wavelength)
                    #wavelength = 10**(wavelength_log)
                    obj_name = spec_name
                elif(sdss_new_spec == False):

                    #trying the older release of sdss
                    t = hdu[1].data
                    flux_tmp = t['flux']
                    wavelength_tmp = t['wavelength']
                    wavelength = np.array(wavelength_tmp)
                    flux = np.array(flux_tmp)
                    #wavelength = 10**(wavelength_log)
                    obj_name = spec_name
                    
                
    #remove NaN values
    not_nan_mask = ~np.isnan(flux)
    nan_mask = np.isnan(flux)

    len_flux1 = len(flux)
    len_wave1 = len(wavelength)
    wavelength = wavelength[not_nan_mask]
    flux = flux[not_nan_mask]
    len_flux2 = len(flux)
    len_wave2 = len(wavelength)
    

    if len_flux1 != len_flux2:
        print ('NaN values found and deleted')
    
    spec_components = len(flux)
    #convert all to nm
    if(lambda_units == 'mu'):
        wavelength = wavelength *1000.
    elif(lambda_units == 'A' or lambda_units == 'a'):
        wavelength = wavelength/10.
    
    #calculating the step
    original_step = wavelength[1]-wavelength[0]
    
    #check on the sampling
    step_chk1 = wavelength[1]-wavelength[0]
    step_chk2 = wavelength[spec_components-1]-wavelength[spec_components-2]
    step_diff = abs(step_chk2-step_chk1)
    eps_step = 2.e-4

    if (step_diff > eps_step):
        print ('Warning: step not constant!')
        #force the linear rebinning because SPAN works better with that
        wavelength, flux, points_spec = resample(wavelength, flux, original_step)
        print('Resampled to linear step')
    return wavelength, flux, original_step, obj_name

#*****************************************************************************************************
# 2) RESAMPLE FUNZIONA. Rebin a spectrum to a useer defined, linear step. INPUT: wavelength, flux, new_step to resample, in the units given in the spectra. OUTPUT: resampled wavelength, resampled flux, number of points of the resampled quantities.
def resample(wavelength, flux, new_step):
    wave_components = len(wavelength)
    initial_wave = wavelength[0]
    final_wave = wavelength[wave_components-1]
    npoint_resampled = int(round((final_wave-initial_wave)/new_step))
    res_wave = np.zeros(npoint_resampled)
    value = 0
    
    for i in range(npoint_resampled):
        res_wave[i] = (value*new_step + initial_wave)
        value = value + 1

    #Now calculate the interpolated flux at the res_wave points
    interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear')
    res_flux = (interpfunc(res_wave))
    spec_components = npoint_resampled
    return res_wave, res_flux, npoint_resampled
    
    
#*****************************************************************************************************
# 3) NORM SPEC. FUNZIONA!
def norm_spec(wavelength, flux, wave_norm, epsilon_avg, flux_to_normalize):
    flux_tmp = 0.
    points_wavelength = len(wavelength) 
    norm_flux = np.zeros(points_wavelength)
    points_flux_to_normalize = len(flux_to_normalize)
    points_flux_to_normalize = len(flux_to_normalize)
    npoints = 0
    
    for i in range(points_wavelength):
        if (wavelength[i] >= (wave_norm - epsilon_avg) and wavelength[i] <= (wave_norm + epsilon_avg)):
                flux_tmp = flux_tmp + flux[i]
                npoints = npoints + 1
 
    avg_flux_ref = flux_tmp/npoints
    for i in range(points_flux_to_normalize):
        norm_flux[i] = (flux_to_normalize[i]/avg_flux_ref)
    return norm_flux
      
      
#*****************************************************************************************************                
# 4) EXTRACT INDEX. SEMBRA FUNZIONARE
def extract_index(wavelength, flux, index):
    left_band = index[0]
    right_band = index[3]
    left_index = index[4]
    right_index = index[5]
    epsilon = 5. #arbitrary values in nm to add to the index window
    index_flux = []
    index_wave = []
    wave_components = len(wavelength)
    
    #;case 1): Line band in the middle
    if (left_band <= left_index and left_index <= right_band): 
        left_band_tmp = left_band-epsilon
        right_band_tmp = right_band+epsilon 
        left_index_tmp = left_index

    #;case 2): Continuum bands on the left
    if (left_band <= left_index and left_index >= right_band):
        left_band_tmp = left_band -epsilon
        right_band_tmp = right_index+epsilon
        right_index_tmp = right_band

#;case 3): Continuum bands on the right
    if (left_band >= left_index and left_band <= right_band):
        left_band_tmp = left_index-epsilon
        right_band_tmp = right_band+epsilon
        left_index_tmp = left_band
    
    for i in range (wave_components): 
        if (wavelength[i] >= left_band_tmp and wavelength[i] <= right_band_tmp):
            index_flux.append(flux[i])
            index_wave.append(wavelength[i])
            #counter = counter +1
    return index_wave, index_flux
  

#*****************************************************************************************************
# 5) Generatind a pseudo continuum, useful for the Equivalent Width (EW) calculation
def idx_cont(index, wavelength, flux):
    left_wave_a = index[0]
    left_wave_b = index[1]
    right_wave_a = index[2]
    right_wave_b = index[3]
    index_left_band = index[4]
    index_right_band = index[5]
    index_central_band = (index_left_band+index_right_band) /2.
    nx = len(wavelength)
    
    cont_flux_left = 0.
    points_left = 0.
    cont_flux_right = 0.
    points_right = 0.

    #;find the total flux on the blue and red continuum band
    for i in range(nx):
        if (wavelength[i] >= left_wave_a and wavelength[i] <= left_wave_b):
            cont_flux_left = cont_flux_left + flux[i]
            points_left += 1
        if (wavelength[i] >= right_wave_a and wavelength[i] <= right_wave_b):
            cont_flux_right = cont_flux_right + flux[i]
            points_right += 1


    #;finding the average blue and red flux
    avg_left_flux = cont_flux_left/points_left
    avg_right_flux = cont_flux_right/points_right

    #;finding the central blue and red band of the continuum
    central_left_wave = (left_wave_a + left_wave_b) /2.
    central_right_wave = (right_wave_a + right_wave_b ) /2.

    #;add the data to arrays
    flux_ref_cont = [avg_left_flux, avg_right_flux]
    wave_ref_cont = [central_left_wave, central_right_wave]

    #;finding the x points (lambda) where interpolate the continuum 
    #points = 0.
    wave_pseudo_cont = []
    line_flux = []

    for i in range(nx): 
        if (wavelength[i] >= index_left_band and wavelength[i] <= index_right_band): 
            wave_pseudo_cont.append(wavelength[i])
            line_flux.append(flux[i])

    #doing the interpolation here
    interpfunc = interpolate.interp1d(wave_ref_cont, flux_ref_cont, kind = 'linear',fill_value='extrapolate')
    flux_pseudo_cont = (interpfunc(wave_pseudo_cont))
        
    return flux_pseudo_cont, line_flux, wave_pseudo_cont, flux_ref_cont, wave_ref_cont, central_left_wave, central_right_wave, avg_left_flux, avg_right_flux


#*****************************************************************************************************
#6) Calculating the EW of an index in A
def eq_width(flux_pseudo_cont, line_flux, step):
    ew = 0.
    points = len(line_flux)
    for i in range(points):
        ew = ew + (1- line_flux[i]/flux_pseudo_cont[i])*step
    return ew

#6) Calculating the EW of an index in mag
def eq_width_mag(flux_pseudo_cont, line_flux, step,lambda_blue_line, lambda_red_line):
    ew_mag = 0.
    points = len(line_flux)
    for i in range(points):
        ew_mag = ew_mag + (line_flux[i]/flux_pseudo_cont[i]*step)
    #Check for negative argument of the log!
    if ew_mag <= 0:
        ew_mag = 999
    else:
        ew_mag = -2.5*mt.log10(1/(lambda_red_line-lambda_blue_line)*ew_mag)
    return ew_mag

#*****************************************************************************************************
#7) Uncertainties of the EW values
def ew_err (index, wavelength, flux, step, flux_pseudo_cont, wave_pseudo_cont, flux_ref_cont, wave_ref_cont):
    nx = len(wavelength)
    nx_check = len(flux)
    if (nx != nx_check):
        print ('ERROR: wavelength componens different from flux. Something went wrong!')
    
    #;definying useful arrays and variables
    flux_red_band = []
    flux_blue_band = []
    red_wave = []
    blue_wave = []
    h = 0

    #;extract the flux and lambda arrays from the continuum blue and red bands. 
    for i in range(nx):
        if (wavelength[i] >= index[0] and wavelength[i] <= index[1]):
            flux_blue_band.append(flux[i])
            blue_wave.append(wavelength[i])
        
        if (wavelength[i] >= index[2] and wavelength[i] <= index[3]):
            flux_red_band.append(flux[i])
            red_wave.append(wavelength[i])
                    
    wave_all = []
    flux_all = []
    
    for i in range(len(blue_wave)):
        wave_all.append(blue_wave[i])
        flux_all.append(flux_blue_band[i])
        
    for j in range(len(red_wave)):
        wave_all.append(red_wave[j])
        flux_all.append(flux_red_band[j])
    
    #interpolate to create the a fake continuum where measure the EW 
    interpfunc = interpolate.interp1d(wave_ref_cont, flux_ref_cont, kind = 'linear', fill_value='extrapolate')
    interp_all = (interpfunc(wave_all))
    
    #calculating the residual continuum
    components = len(flux_all)
    residuals = []
    for i in range(components):
        residuals.append(flux_all[i] - interp_all[i])

    sigma_cont_real = np.std(residuals)
    
    scale = sigma_cont_real
    line_components = len(flux_pseudo_cont)
    
    number_noisy_cont = 1000 #how many syntethics? a lot!
    ews= []
    
    
    for k in range(number_noisy_cont):
        noise_array = np.random.standard_normal((line_components,))
        noise_array_scaled = noise_array * scale
        noisy_cont = []
        
        for i in range(line_components):
            noisy_cont.append(flux_pseudo_cont[i] + noise_array_scaled[i])
    
        #;measuring the EW
        ew = eq_width(flux_pseudo_cont, noisy_cont, step)
        ews.append(ew)
    
    error = np.std(ews)
    return error, sigma_cont_real

#*****************************************************************************************************
#8) reading index file
def read_idx(index_file):
    indices = []
    idx_names = np.loadtxt(index_file, dtype = 'str', delimiter = ' ', max_rows = 1) #reading the first line only, with the header
    indices = np.loadtxt(index_file, comments = '#', skiprows = 1) #reading the other lines
    return idx_names, indices


#*****************************************************************************************************
# 9) Degradation of the spectrum
def degrade(wavelength, flux, original_resolution, final_resolution, verbose):
    c = 299792.458
    wave_components = len(wavelength)
    log_wave = []
    initial_wave = wavelength[0]
    final_wave = wavelength[wave_components-1]
    initial_step = wavelength[1]- wavelength[0]
    log_wave.append(wavelength[0])
    tmp_wave = wavelength[0]
    i = 0
    vel_pix = c/wavelength[0]*initial_step
    
    # 1) converto tutto a campionamento log, così che sigma = c/R = cost
    while tmp_wave < final_wave:
        new_step = (vel_pix/c)*log_wave[i]
        log_wave.append(log_wave[i]+new_step)
        tmp_wave = tmp_wave +new_step
        i = i+1

    if verbose == True:
        print ('Resampling to log...')
        print ('New step at lambda initial ', log_wave[0],' : ', (log_wave[1]-log_wave[0]) ,' nm' )
        print ('New step at lambda final ' , log_wave[len(log_wave)-1], ' : ', (log_wave[len(log_wave)-1]-log_wave[len(log_wave)-2]), ' nm')
    
    #Now calculate the interpolated flux at the res_wave points
    interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear',fill_value='extrapolate')
    log_flux = (interpfunc(log_wave))
    original_resolution_vel = c/original_resolution
    
    #converting from FHWN values to stdev since the gaussian broadening wants this value
    fwhm_to_sigma = 2.3548
    sigma_original_resolution = original_resolution_vel/fwhm_to_sigma
    final_resolution_vel = c/final_resolution
    sigma_final_resolution = final_resolution_vel/fwhm_to_sigma

    if verbose == True:
        print ('Original resolution in sigma vel: ', sigma_original_resolution)
        print ('Final resolution in sigma vel: ', sigma_final_resolution)
        print ('Velocity per pixel: ', vel_pix )

    
    #2) Degrading the resolution only if the resolution selected is smaller than the initial
    if ((sigma_final_resolution-sigma_original_resolution)> 0):
        sigma_to_broad = mt.sqrt(sigma_final_resolution**2-sigma_original_resolution**2)
    
        gauss_stdev_pix = sigma_to_broad/vel_pix
        
        if verbose == True:
            print ('Sigma to broad the spectra: ', sigma_to_broad)
            print ('Gaussian kernel sigma to apply (in pixels): ', gauss_stdev_pix)
        
        kernel = Gaussian1DKernel(gauss_stdev_pix)
        log_degraded_flux = convolve(log_flux, kernel)
        
        #testing an alternative: works the same as the previous!
        #log_degraded_flux = gaussian_filter1d(log_flux,gauss_stdev_pix)
    
    else:
        print('Final resolution is less than the original: cannot degrade the spectrum. Doing nothing.')
        return wavelength, flux

    
    #3) Resampling the spectrum to the original step
    interpfunc1 = interpolate.interp1d(log_wave, log_degraded_flux, kind = 'linear')
    degraded_flux = (interpfunc1(wavelength))

    return wavelength, degraded_flux


# 9) BIS: DEGRADATION OF THE SPECTRA TO DELTA LAMBDA
def degrade_lambda(wavelength, flux, original_resolution_lambda, final_resolution_lambda):
    fwhm_to_sigma = 2.3548
    wave_components = len(wavelength)
    step = wavelength[1]-wavelength[0] # STEP SUPPOSED LINEAR!
    
    original_resolution_lambda_nm = original_resolution_lambda/10. #converting to nm!
    final_resolution_lambda_nm = final_resolution_lambda/10. #converting to nm!
    
    
    print ('Original resolution: ', original_resolution_lambda, 'A')
    print ('Final resolution in sigma vel: ', final_resolution_lambda , 'A')
    

    #2) Degrading the resolution only if the resolution selected is smaller than the initial
    if ((final_resolution_lambda_nm - original_resolution_lambda_nm)> 0):
        real_value_to_broad = mt.sqrt(final_resolution_lambda_nm**2-original_resolution_lambda_nm**2)
    
        gauss_fwhm_pix = real_value_to_broad/step
        gauss_stdev_pix = gauss_fwhm_pix/fwhm_to_sigma
        
        print ('Sigma to broad the spectra: ', real_value_to_broad*10, 'A')
        print ('Gaussian kernel sigma to apply (in pixels): ', gauss_stdev_pix)
        
        kernel = Gaussian1DKernel(gauss_stdev_pix)
        degraded_flux = convolve(flux, kernel)
        
        #testing and alternative: works the same as the previous!
        #log_degraded_flux = gaussian_filter1d(log_flux,gauss_stdev_pix)
    
    else:
        print('Final resolution is less than the original: cannot degrade the spectrum. Doing nothing.')
        return wavelength, flux

    return wavelength, degraded_flux



# 9) TRIS: DEGRADATION OF THE SPECTRA FOR LICK INDICES MEASUREMENT
def degrade_to_lick(wavelength, flux, original_resolution, res_delta_lambda):
    fwhm_to_sigma = 2.3548
    wave_components = len(wavelength)
    step = wavelength[1]-wavelength[0] # STEP SUPPOSED LINEAR!

    #IF THE RESOLUTION GIVEN IS IN FWHM!
    if res_delta_lambda == True:
        original_resolution_lambda_nm = original_resolution/10. #converting to nm!
        final_resolution_lambda_nm_lick = 8.4/10. #converting to nm!


        print ('Original resolution: ', original_resolution, 'A')

        #2) Degrading the resolution only if the resolution selected is smaller than the initial
        if ((final_resolution_lambda_nm_lick - original_resolution_lambda_nm)> 0):

            real_value_to_broad = mt.sqrt(final_resolution_lambda_nm_lick**2-original_resolution_lambda_nm**2)

            gauss_fwhm_pix = real_value_to_broad/step
            gauss_stdev_pix = gauss_fwhm_pix/fwhm_to_sigma

            print ('Sigma to broad the spectra: ', real_value_to_broad*10, 'A')
            print ('Gaussian kernel sigma to apply (in pixels): ', gauss_stdev_pix)

            kernel = Gaussian1DKernel(gauss_stdev_pix)
            degraded_flux = convolve(flux, kernel)

            #testing and alternative: works the same as the previous!
            #log_degraded_flux = gaussian_filter1d(log_flux,gauss_stdev_pix)

        else:
            print('Final resolution is less than the original: cannot degrade the spectrum. Doing nothing.')
            degraded_flux = flux
            return wavelength, degraded_flux

        return wavelength, degraded_flux



    #IF THE RESOLUTION GIVEN IS IN R!
    if res_delta_lambda == False:
        #original_resolution_lambda_nm = original_resolution/10. #converting to nm!
        final_resolution_lambda_nm_lick = 8.4/10. #converting to nm!

        original_resolution_lambda_nm = wavelength/original_resolution #array contenente le risoluzioni in FWHM
        final_resolution_lambda_nm = np.full_like(original_resolution_lambda_nm, final_resolution_lambda_nm_lick, dtype=float)#Array con le stesse dimensioni contenente la risoluzione degli indici di Lick, in FWHM
        real_value_to_broad = np.zeros_like(wavelength)
        degraded_flux = np.zeros_like(flux)

        for i in range (len(wavelength)):

            real_value_to_broad[i] = mt.sqrt(final_resolution_lambda_nm[i]**2-original_resolution_lambda_nm[i]**2)

        gauss_sigma = real_value_to_broad/fwhm_to_sigma

        #using the varsmooth function of ppxf.util for convolution with variable sigma. Works great!
        degraded_flux = util.varsmooth(wavelength, flux, gauss_sigma, xout=None, oversample=1)
        return wavelength, degraded_flux


#*****************************************************************************************************
# 10) Rough continuum subtraction by degrading the spectrum to a very low resolution and considering this as continuum
def sub_cont(wavelength, flux, operation):
    cont_wave, cont_flux = degrade(wavelength, flux, 5000., 50., False) # degrade the spectra to a fixed value
    if operation == 'divide':
        norm_flux = flux/cont_flux #I normalise the continuun, not subtract
    if operation == 'subtract':
        norm_flux = flux - cont_flux #I normalise the continuun, not subtract
    return norm_flux, cont_flux

#*****************************************************************************************************
# 11) Correction EW for the sigma broadening
def corr_ew(ew_file, corr_file, sigma_file):
    #reading the files
    data_ew = pd.read_csv(ew_file, header=None, sep = ' ')
    data_corr = pd.read_csv(corr_file, header=None, sep = ' ')
    
    data_number = len(data_ew.index)

    index_number = round((len(data_ew.columns)-1)/2)
    total_column = len(data_ew.columns)
    ew_values_starting_at = 1
    ew_values_end_at = index_number+1
    err_ew_starting_at = ew_values_end_at
    err_ew_ends_at = total_column
    n_column_ew_to_correct = index_number
    
    column_corr_file = len(data_corr.columns)
    corr_ew_ends_at = round(column_corr_file/2)
    coeff_number = len(data_corr.index)
 
    ew_values = pd.DataFrame(data_ew.iloc[1:data_number, ew_values_starting_at:ew_values_end_at])
    err_values = pd.DataFrame(data_ew.iloc[1:data_number, err_ew_starting_at:err_ew_ends_at ])
    index_names = pd.DataFrame(data_ew.iloc[0:1, ew_values_starting_at:err_ew_ends_at ])
    all_idx = pd.DataFrame(data_ew.iloc[0:1, 0:err_ew_ends_at ])
    sigma_values = np.loadtxt(sigma_file, usecols = [1])
    
    spec_names = pd.DataFrame(data_ew.iloc[1:data_number, 0:1])
    
    corr_ew_values = pd.DataFrame(data_corr.iloc[1:coeff_number, 0:corr_ew_ends_at])
    corr_err_values = pd.DataFrame(data_corr.iloc[1:coeff_number, corr_ew_ends_at:column_corr_file])
    corr_index_names = pd.DataFrame(data_corr.iloc[0:1, 0:column_corr_file])

    ew_values_np = ew_values.to_numpy(dtype = float)
    err_values_np = err_values.to_numpy(dtype = float)
    sigma_values_np = sigma_values
    #rec_vel_values_np = rec_vel_values.to_numpy(dtype = float)
    corr_ew_values_np = corr_ew_values.to_numpy(dtype = float)
    corr_err_values_np = corr_err_values.to_numpy(dtype = float)
    corr_index_names_np = corr_index_names.to_numpy(dtype = str)
    spec_names_np = spec_names.to_numpy()
    index_names_np = index_names.to_numpy(dtype = str)    
    all_idx_np = np.loadtxt(ew_file, dtype = 'str', delimiter = ' ' , max_rows = 1, comments = '##')

    spectra_number = data_number -1
    new_ew = np.zeros((spectra_number, n_column_ew_to_correct))
    ew_correction = np.zeros((spectra_number, n_column_ew_to_correct))
    sigma_correction = np.zeros((spectra_number, n_column_ew_to_correct))
    new_err_tot = np.zeros((spectra_number, n_column_ew_to_correct))
    
    #wiping the zeros
    epsilon = 1e-5
    ew_values_np[ew_values_np == 0] = epsilon

    #checking correspondence between the files
    is_equal = np.char.equal(corr_index_names_np, index_names_np)
    if False in is_equal:
        print ('Ops: seems that theere is no correspondence between the ew file and the correction one!')
        print ('Doing nothing')
        return 0
    
    for i in range  (n_column_ew_to_correct):
        for j in range (spectra_number):
            
            #correcting ews
            ew_correction[j,i] = (corr_ew_values_np[3,i] + corr_ew_values_np[2,i]*sigma_values_np[j]+corr_ew_values_np[1,i]*sigma_values_np[j]**2+corr_ew_values_np[0,i]*sigma_values_np[j]**3)
            new_ew[j,i] = ew_values_np[j,i]/ (ew_correction[j,i]+1)
            
            #correcting uncertainties
            sigma_correction[j,i] = abs(corr_err_values_np[3,i] + corr_err_values_np[2,i]*sigma_values_np[j]+corr_err_values_np[1,i]*sigma_values_np[j]**2+corr_err_values_np[0,i]*sigma_values_np[j]**3)
            
            ##total uncertainties
            new_err_tot[j,i] =  mt.sqrt(sigma_correction[j,i]**2 + (err_values_np[j,i]/abs(ew_values_np[j,i]))**2)*abs(new_ew[j,i]) #c'è uno zero o qualcosa che non va nella righa 9
            
    #Ok, putting together
    new_ew_data = np.column_stack((spec_names_np, new_ew, new_err_tot))

    return all_idx_np, new_ew_data

#*************************************************************************************************
#12 Sigma clipping
def sigma_clip(wavelength, flux, clipping, resolution, sigma_vel):
    c = 299792.458
    npoints = len(wavelength)
    base_smooth_width = 41
    counter = 1
    original_flux = flux
    clipped_flux = np.zeros(npoints)
    
    if sigma_vel == 0:
        dynamic_size = 41
        sigma_coeff = 0
    else:
        sigma_instrum = c/resolution
        sigma_coeff = sigma_vel/sigma_instrum
        dynamic_size = round(base_smooth_width*sigma_coeff)
        if sigma_coeff < 1:
            dynamic_size = int(41)
            window_width_small = 1
        else:
            small_frequency_size = 51
            window_width_small = small_frequency_size
    
    window_width_dyn = int(dynamic_size)
    
    #se la window width è pari, falla dispari:
    if (2*(window_width_small/2) == window_width_small):
        window_width_small =  window_width_small + 1
    
    clipping_dyn = clipping
    clipping_small = 2.5
    
    elements_flux = len(flux)
    clipped_flux = np.zeros(elements_flux)
    n_clipped = 0
    sigma_flux_ext_small = np.zeros(window_width_small)
    sigma_flux_ext_dyn = np.zeros(window_width_dyn)
    
    print ('Velocity dispersion: ', sigma_vel)
    print ('Smooth windows: ', window_width_small, window_width_dyn)

    original_flux = flux
    counter = 1
    
    
    step1 = wavelength[1]-wavelength[0]
    step2 = wavelength[npoints-1]-wavelength[npoints-2]
    epsilon = 1e-4
    
    #If the step is not linear, resampling to linear
    #if abs(step2-step1) > epsilon:
        #new_step = step1
        #res_wave, res_flux, npoint_resampled = resample (wavelength, flux, new_step)
        #wavelength = res_wave
        #flux = res_flux
    

    #pre cleaning at small frequencies, only for sigma_vel greater than 100 or more
    if (sigma_coeff > 1):
        threshold = 50
        cond = 1
        print ('Performing pre-cleaning with smooth size: ', window_width_small)
        while (threshold > 0):
            n_clipped = 0
            smooth_flux = scipy.ndimage.filters.uniform_filter(flux, size = window_width_small)
            
            #For all the flux values
            for k in range (elements_flux):
                
                #clipping the edges
                if (k >= window_width_small and k <= (elements_flux-1-window_width_small)):
                    
                    #calculating sigma for the smoothed window
                    for i in range (window_width_small):
                        sigma_flux_ext_small[i] = flux[int(k-(window_width_small-1)/2+i)]
                    
                    sigma = np.std(sigma_flux_ext_small)
                    
                    #clipping
                    if (abs(smooth_flux[k]-flux[k]) <= sigma*clipping_small):
                        clipped_flux[k] = flux[k]
                    else:
                        clipped_flux[k] = smooth_flux[k]
                        n_clipped = n_clipped + 1
        
                #On the edges I do nothing
                if (k < window_width_small or k > (elements_flux-1-window_width_small)):
                    clipped_flux[k] = flux[k]

            print ('Iteration N.', counter, ' Clipped points: ', n_clipped)
            flux = clipped_flux
            counter = counter + 1
            threshold = n_clipped
            
            
    # dinamycal clipping with windows size based on the sigma_vel or sigma instrumental
    print ('')
    print ('Performing cleaning with dynamical smooth size: ', window_width_dyn)
    counter = 1
    
    
    threshold = 50
    while (threshold > 0):
        
        #if true the cycle will stop
        n_clipped = 0
    
        #Smoothing
        smooth_flux = scipy.ndimage.filters.uniform_filter(flux, size = window_width_dyn)
        
        #ciclo su tutti i valori di flusso
        for k in range (elements_flux):
            
            #condizione per escludere i bordi, brutalmente li taglio via dalla valutazione del sigma clipping
            if (k >= window_width_dyn and k <= (elements_flux-window_width_dyn-1)):
                
                #calculating sigma for the smoothed window
                for i in range (window_width_dyn):
                    sigma_flux_ext_dyn[i] = flux[int(k-(window_width_dyn-1)/2+i)]
                
                sigma = np.std(sigma_flux_ext_dyn)
    
                #clipping
                if (abs(smooth_flux[k]-flux[k]) <= sigma*clipping_dyn):
                    clipped_flux[k] = flux[k]
                else:
                    clipped_flux[k] = smooth_flux[k]
                    n_clipped = n_clipped + 1
                    #cond = 0 #until this is zero the cycle will proceed!
    
            #se sono sui bordi non faccio nulla!
            if (k < window_width_dyn or k > (elements_flux-1-window_width_dyn)):
                clipped_flux[k] = flux[k]
        
        print ('Iteration N.', counter, ' Clipped points: ', n_clipped)
        flux = clipped_flux
        counter = counter + 1
        threshold = n_clipped
    
         
    new_flux = flux
    return wavelength, new_flux

#*************************************************************************************************
#13 Log rebin
def log_rebin(wavelength, flux, sigma_pix):
    c = 299792.458
    wave_components = len(wavelength)
    log_wave = []
    initial_wave = wavelength[0]
    final_wave = wavelength[wave_components-1]
    initial_step = wavelength[1]- wavelength[0]
    log_wave.append(wavelength[0])
    tmp_wave = wavelength[0]
    i = 0
    
    #if sigma_pix input value == 0, I use the initial step
    if sigma_pix == 0:
        sigma_pix = c/wavelength[0]*initial_step
    
        
    # 1) converting to log resample so sigma = c/R = cost
    while tmp_wave < final_wave:
        new_step = (sigma_pix/c)*log_wave[i]
        log_wave.append(log_wave[i]+new_step)
        tmp_wave = tmp_wave +new_step
        i = i+1
        
     #Now calculate the interpolated flux at the res_wave points
    log_wave = np.array(log_wave)
    interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear',fill_value='extrapolate')
    log_flux = (interpfunc(log_wave))
    
    return log_wave, log_flux

#*************************************************************************************************
#14 Doppler correction
def dopcor(wavelength, flux, rec_vel):
    c = 299792.458
    z = rec_vel/c
    npoints = len(wavelength)
    corr_wave = np.zeros(npoints)
    for i in range(npoints): 
        corr_wave[i] = wavelength[i]/(1+z)
    #rebinning to the initial step
    in_step = corr_wave[1]-corr_wave[0]
    res_wave, res_flux, npoint_resampled = resample(corr_wave, flux, in_step)
    return res_wave, res_flux
    
#*************************************************************************************************
#15 Average spectra
def average(lambda_units, spectra_number, spec_names):
    #df_list_spec = pd.read_csv(list_spec, delimiter = ' ')
    #spectra_number = len(df_list_spec.index)
    #spec_names = np.loadtxt(list_spec, dtype = 'str', delimiter = ' ', usecols=[0])
    
    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = read_spec(spec_names[i], lambda_units)
        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = flux
        if i > 0:

            #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))
            df_to_avg.iloc[:,i] = interp_flux
    average = df_to_avg.mean(axis=1) #average all over the columns
    average_flux = average.to_numpy(dtype = float)
    average_spec = np.column_stack((wavelength_grid, average_flux))
    return average_spec


#*************************************************************************************************
#16 Normalize and average spectra
def average_norm(lambda_units, wavelength, flux, spectra_number, spec_names):
    #df_list_spec = pd.read_csv(list_spec, delimiter = ' ')
    #spectra_number = len(df_list_spec.index)
    #spec_names = np.loadtxt(list_spec, dtype = 'str', delimiter = ' ', usecols=[0])
    
    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0
    
    wavelength_grid = wavelength
    npoints = len(wavelength_grid)
    step = wavelength_grid[1]-wavelength_grid[0]
    wave_norm = (wavelength_grid[0] + wavelength_grid[npoints-1])/2.
    epsilon_norm = 10*step

    #reference spectrum is the one selected in span:
    norm_flux_reference_spec = norm_spec(wavelength, flux, wave_norm, epsilon_norm, flux)

    shape = (len(wavelength_grid), spectra_number)
    data_to_df = np.zeros(shape)
    df_to_avg = pd.DataFrame(data_to_df)
            
    for i in range (spectra_number):
        wavelength, flux, step, name = read_spec(spec_names[i], lambda_units)

        #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
        interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
        interp_flux = (interpfunc(wavelength_grid))
        norm_interp_flux = norm_spec(wavelength_grid, interp_flux, wave_norm, epsilon_norm, interp_flux)
        df_to_avg.iloc[:,i] = norm_interp_flux
        
    average = df_to_avg.mean(axis=1) #average all over the columns
    average_flux = average.to_numpy(dtype = float)
    average_norm_spec = np.column_stack((wavelength_grid, average_flux))
    return average_norm_spec


#*************************************************************************************************
#17 Sum the spectra
def sum_spec(lambda_units, spectra_number, spec_names):
    #df_list_spec = pd.read_csv(list_spec, delimiter = ' ')
    #spectra_number = len(df_list_spec.index)
    #spec_names = np.loadtxt(list_spec, dtype = 'str', delimiter = ' ', usecols=[0])
    
    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = read_spec(spec_names[i], lambda_units)
        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = flux

        if i > 0:
            #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))
            df_to_avg.iloc[:,i] = interp_flux
            
    sum_all = df_to_avg.sum(axis=1) #average all over the columns
    sum_flux = sum_all.to_numpy(dtype = float)
    sum_spec = np.column_stack((wavelength_grid, sum_flux))
    return sum_spec



#*************************************************************************************************
#18 Normalize and sum the spectra
def sum_norm_spec(lambda_units, spectra_number, spec_names):
    #df_list_spec = pd.read_csv(list_spec, delimiter = ' ')
    #spectra_number = len(df_list_spec.index)
    #spec_names = np.loadtxt(list_spec, dtype = 'str', delimiter = ' ', usecols=[0])
    
    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0
    
    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = read_spec(spec_names[i], lambda_units)
        epsilon_norm = 10*step
        npoints = len(wavelength)
        
        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)
            
            #find the normalization wavelength = half way to the beginning and the end of the spectrum
            wave_norm = (wavelength_grid[0] + wavelength_grid[npoints-1])/2. 
            norm_flux = norm_spec(wavelength_grid, flux, wave_norm, epsilon_norm, flux)
            
            #storing data
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = norm_flux


        if i > 0:
            #finding the flux at the wavelength_grid location
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))
            
            #normalizing the flux
            norm_flux = norm_spec(wavelength, interp_flux, wave_norm, epsilon_norm, interp_flux)
            
            #storing data
            df_to_avg.iloc[:,i] = norm_flux
            
    sum_norm_all = df_to_avg.sum(axis=1) #average all over the columns
    sum_norm_flux = sum_norm_all.to_numpy(dtype = float)
    sum_norm_spec = np.column_stack((wavelength_grid, sum_norm_flux))
    return sum_norm_spec

#*************************************************************************************************
#19 Subtract normalized average

def sub_norm_avg(wavelength, flux, lambda_units, spectra_number, spec_names):
        

    #averaging all spectra
    average_norm_spec = average_norm(lambda_units, wavelength, flux, spectra_number, spec_names)
    average_wavelength = average_norm_spec[:,0]
    average_flux = average_norm_spec[:,1]
        
    # interpolating the average flux to the wavelength range of the spectrum
    interpfunc = interpolate.interp1d(average_wavelength, average_flux, kind = 'linear', fill_value='extrapolate')
    interp_average_flux = (interpfunc(wavelength))
        
    #normalize the spectra
    npoints = len(wavelength)
    wave_normalization = (wavelength[0] + wavelength[npoints-1])/2. 
    epsilon_norm = (wavelength[1]-wavelength[0])*10
    norm_flux = norm_spec(wavelength, flux, wave_normalization, epsilon_norm, flux)
    norm_interp_average = norm_spec(wavelength, interp_average_flux, wave_normalization,epsilon_norm, interp_average_flux)
    subtracted_flux = norm_flux - norm_interp_average
    return subtracted_flux

#*************************************************************************************************
#20 Subtract normalised spec
def sub_norm_single(wavelength, flux, spectrum_to_subtract, lambda_units):
    wave_to_sub, flux_to_sub, step, name = read_spec(spectrum_to_subtract, lambda_units)

    # interpolating the average flux to the wavelength range of the spectrum
    interpfunc = interpolate.interp1d(wave_to_sub, flux_to_sub, kind = 'linear', fill_value='extrapolate')
    interp_flux_to_sub = (interpfunc(wavelength))
        
    #normalise the spectra
    npoints = len(wavelength)
    wave_normalization = (wavelength[0] + wavelength[npoints-1])/2. 
    epsilon_norm = (wavelength[1]-wavelength[0])*10

    #normalisation
    norm_flux = norm_spec(wavelength, flux, wave_normalization, epsilon_norm, flux)
    norm_interp_sub = norm_spec(wavelength, interp_flux_to_sub, wave_normalization, epsilon_norm, interp_flux_to_sub)
    subtracted_flux = norm_flux - norm_interp_sub
    return subtracted_flux


#*************************************************************************************************
#21 Sigma broadening: broad the spectra to a user defined sigma value (km/s)
def sigma_broad(wavelength, flux, sigma_to_broad):
    c = 299792.458
    step = wavelength[1]-wavelength[0]
    sigma_pix = c/wavelength[0]*step
    
    if sigma_to_broad == 0:
        flux_gauss = flux
        return flux_gauss
    
    #rebinning to sigma = cost
    log_wave, log_flux = log_rebin(wavelength,flux, sigma_pix)

    gauss_stdev_pix = sigma_to_broad/sigma_pix
    kernel = Gaussian1DKernel(gauss_stdev_pix)
    log_flux_gauss = convolve(log_flux, kernel)
    
    #back to the original, linear sample
    interpfunc1 = interpolate.interp1d(log_wave, log_flux_gauss, kind = 'linear')
    flux_gauss = (interpfunc1(wavelength))
    
    return flux_gauss

#*************************************************************************************************
#22 Add noise
def add_noise(wavelength, flux, snr):
    
    #normalise the continuum. This gives me the possibility to have the new spectrum at level = 1, usefull for adding noise
    norm_flux, cont_flux = sub_cont(wavelength, flux, 'divide')
    npoints = len(wavelength)
    
    #generate the noise array
    noise_array = np.random.standard_normal((npoints,))
    scale = 1./snr
    noise_array_scaled = noise_array * scale
    noisy_norm_flux = norm_flux + noise_array_scaled
    
    #give back the continuum shape
    noisy_flux = noisy_norm_flux * cont_flux
    return noisy_flux
    
#*************************************************************************************************
#23 Cross-correlation, use with cautions for galaxies. Works well for stellar spectra
def crosscorr (wavelength_spec, flux_spec, template, lambda_units_template, wave_interval, smooth_vel, vel_interval):
    
    #reading the template
    wavelength_template_orig, flux_template_orig, step, name = read_spec(template, lambda_units_template)
    
    #Sampling to the common, smallest, delta lambda
    step_spec = wavelength_spec[1]-wavelength_spec[0]
    step_template = step
    
    if step_spec > step_template:
        wavelength, flux, npoint_resampled = resample(wavelength_spec, flux_spec, step_template)
        
        #resampling to linear also the template
        wavelength_template, flux_template, npoint_template = resample(wavelength_template_orig, flux_template_orig, step_template)
    
    elif step_spec <= step_template:
        wavelength_template, flux_template, npoint_template = resample(wavelength_template_orig, flux_template_orig, step_spec)
        #resampling to linear also the spec
        wavelength, flux, npoint_resampled = resample(wavelength_spec, flux_spec, step_spec)
        
    #Extracting the arrays values, considering that they can be inverted, so I look for max and min values)
    low_lim_vel =  np.min(vel_interval)
    high_lim_vel = np.max(vel_interval)
    low_wave_interval = np.min(wave_interval)
    high_wave_interval = np.max(wave_interval)
    step_vel = 2.
    
    #normalising to a mid wavelength
    npoints_tot = len(wavelength)
    
    if low_wave_interval == 0:
        wave_norm_spec = (wavelength[0]+wavelength[npoints_tot-1])/2.
    else: 
        wave_norm_spec = (low_wave_interval+high_wave_interval)/2.
        
    epsilon_norm = step*10.
    norm_flux = norm_spec(wavelength, flux, wave_norm_spec, epsilon_norm, flux)
    norm_flux_template = norm_spec(wavelength_template, flux_template, wave_norm_spec, epsilon_norm, flux_template)
    flux = norm_flux
    flux_template = norm_flux_template
    
    #smoothing the template
    if(smooth_vel > 0):
        flux_gauss_template = sigma_broad(wavelength_template, flux_template, smooth_vel)
        flux_template = flux_gauss_template
        
    
    #if thew interval is zero:
    if low_wave_interval == 0:
        rv, cc = crosscorrRV(wavelength, flux, wavelength_template, flux_template, low_lim_vel, high_lim_vel, step_vel)
        
        #estrapolating the most probable rv value
        max_corr_fcn = np.argmax(cc)
        rv_at_max = rv[max_corr_fcn]
        
    #otherwise I need to extract lambda and flux both for the spectrum and the template.
    else:
        corr_wave = []
        corr_flux = []
        corr_wave_template = []
        corr_flux_template = []
        
        npoints_wave = len(wavelength)
        for i in range(npoints_wave):
            if wavelength[i]>= low_wave_interval and wavelength[i]<= high_wave_interval:
                corr_wave.append(wavelength[i])
                corr_flux.append(flux[i])

        #converting to numpy because Python sucks!
        corr_wave_np = np.asarray(corr_wave, dtype = float)
        corr_flux_np = np.asarray(corr_flux, dtype = float)
        wavelength = corr_wave_np
        flux = corr_flux
        

        ##calculating the crosscorrelation
        rv, cc = crosscorrRV(wavelength, flux, wavelength_template, flux_template, low_lim_vel, high_lim_vel, step_vel)
        
        #estrapolating the most probable rv value
        max_corr_fcn = np.argmax(cc)
        rv_at_max = rv[max_corr_fcn]
        cc_at_max = cc[max_corr_fcn]
        
    return rv, cc, rv_at_max, cc_at_max, wavelength, flux, wavelength_template, flux_template  #rv_at_max




#*************************************************************************************************
#24 Velocity dispersion measurement
def sigma_measurement(wavelength, flux, spec_test_template, lambda_units_template, resolutionR_spec, resolution_template, banda1, banda1_cont, err_calc):
    c = 299792.458

        #************************************ select the bands/lines ******************************
    wave_range = banda1
    snr_range = banda1_cont
    line_name = 'banda1'
    band_cont = banda1_cont
    wave_norm = np.mean(banda1_cont)
    
    sigma_instrumental = (c/resolutionR_spec)/2.355
    resolution_lambda_fwhm_spec = (np.mean(wave_range)/resolutionR_spec)

    #if the resolution_template variable is zero, I assume they are (E)MILES spectra, otherwise I transform the resolution from R to sigma.
    if resolution_template != 0:
        resolution_lambda_fwhm_temp = (np.mean(wave_range)/resolution_template)
        resolution_temp = (c/resolution_template)/2.355

    #preparing variables for the determination of sigma
    initial_sigma = 0.
    final_sigma = 360. # massima sigma che ha senso fare
    step_sigma = 2. #incremento sigma
    number_values_sigma = round(final_sigma/step_sigma+1)
    sigma_values = np.zeros(number_values_sigma)

    #filling the sigma vector
    for i in range(1, number_values_sigma):
        sigma_values[i] = sigma_values[i-1]+step_sigma

    chisquare_fit = [] #vector for the fit chi squared values



    #********************************* Let's rock! *********************************

    #read the template
    wavelength_template, flux_template, step_template, name = read_spec(spec_test_template, lambda_units_template)


    #test: subtract the continuum. I risultati migliorano. Sembrano più stabili
    flux, cont_flux = sub_cont(wavelength, flux, 'divide')
    flux_template, cont_flux = sub_cont(wavelength_template, flux_template, 'divide')

    #rebin to smallest # better to rebin all to a constant step, in case it's not
    optimal_step = 2*np.mean(wave_range)*step_sigma/c # I tested that this sampling is the optimal one: smaller values don't change the result, while greater values start to feel the quantization effect.

    wavelength_template, flux_template, points = resample(wavelength_template, flux_template, optimal_step)
    wavelength, flux, point_spec = resample(wavelength, flux, optimal_step)

    #uniform the wavelength grid
    interpfunc = interpolate.interp1d(wavelength_template, flux_template, kind = 'linear', fill_value='extrapolate')
    flux_template = (interpfunc(wavelength))

    #storing the original template flux
    flux_template_original = flux_template
    
    #extract the line flux and wavelength arrays
    line_wave = wavelength[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]  
    line_flux_spec = flux[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]
    
    #resolution for the EMILES models: if lambda < 895, the resolution is constant with wavelength, which means that the resolution in sigma diminishes with the increasing wavelength.
    if resolution_template == 0:
        if np.mean(line_wave) < 895.:
            resolution_lambda_fwhm_temp = 0.251 #fwhm
            #converting to sigma
            resolution_temp = resolution_lambda_fwhm_temp/np.mean(line_wave)*c/2.3548 # in sigma velocity. it's an approximation, since the resolution in sigma changes with the wavelength, but if the band is small, the approximation is good.
            print (resolution_temp)
        else:
            resolution_vel_fwhm_temp = 60. #fwhm, in km/s
            resolution_lambda_fwhm_temp = (resolution_vel_fwhm_temp*np.mean(wave_range))/c
            resolution_temp = resolution_vel_fwhm_temp/2.3548
            print (resolution_temp)


    #normalise the spectrum
    epsilon_norm = optimal_step*200
    line_flux_spec_norm = norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #calculating the SNR
    snr_flux = line_flux_spec_norm[(line_wave >= snr_range[0]) & (line_wave <= snr_range[1])]    
    snr = np.mean(snr_flux)/np.std(snr_flux)


    #TEST FITTING TAMPLATE
    #The idea is: select the template, broad to 200 km/s, normalise, extract the working band, extract the working continuum, fitting, finding the displacement with respect to the spectrum, then apply to the original, not broadened and normalized template.

    #preparing the template
    #broadening the template to 200 km/s
    flux_template_broad = sigma_broad(wavelength, flux_template_original, 200.)

    #extract the line flux template
    line_flux_template_orig = flux_template_broad[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]    

    #normalize the template in the wavelength range
    line_flux_temp_test = norm_spec(line_wave, line_flux_template_orig, wave_norm, epsilon_norm, line_flux_template_orig)

    #normalize the whole, original template
    flux_temp_test = norm_spec(wavelength, flux_template_original, wave_norm, epsilon_norm, flux_template_original)

    ##estraggo flusso e lunghezza d'onda del continuo
    cont_spec_banda = line_flux_spec_norm[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]
    cont_temp_banda = line_flux_temp_test[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]
    wave_banda = line_wave[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]
    
    #calculate the rms of the banda1
    rms_banda = np.std(cont_spec_banda)
    #the shift_step will be half or one rms
    shift_step = rms_banda/2.

    starting_flux_level = 0.95
    delta_starting_flux = 1-starting_flux_level
    cont_temp_banda = cont_temp_banda - delta_starting_flux
    end_flux_level = 1.1
    actual_level = starting_flux_level
    chisquare_test = []
    shift_array = []

    #fitting procedure
    while actual_level <= end_flux_level:
        #chi square
        chisqr_test = 0. #it is a chisquared
        for i in range(len(wave_banda)):
            chisqr_test = chisqr_test + (cont_spec_banda[i]-cont_temp_banda[i])**2/cont_temp_banda[i]

        chisquare_test.append(chisqr_test) #chi squared vector!
        shift_array.append(actual_level)
    #test to stop the cycle when reached the minimum chi square, without exploring other values
        actual_level = actual_level + shift_step
        cont_temp_banda = cont_temp_banda + shift_step
        
    min_index_test = np.argmin(chisquare_test)
    best_fit_level_temp = shift_array[min_index_test]
    print ('Adjusting the template to new level:', best_fit_level_temp)

    #applying the corrections
    flux_template_shifted = flux_temp_test + (best_fit_level_temp-1)
    line_flux_temp_test_shifted = line_flux_temp_test + (best_fit_level_temp-1)

    print ('line selected: ', line_name)
    print ('SNR line:', round(snr))
    print('Processing...')


    #********************** FITTING PROCEDURE, BOTH SPECTRAL REGION AND LINES **************************
    shape = len(line_wave)
    line_flux_template_norm = np.zeros(shape)

    for j in range(number_values_sigma):
        
        previous_line_flux_template_norm = line_flux_template_norm
        
        #broadening the template
        flux_template_broad = sigma_broad(wavelength, flux_template_shifted, sigma_values[j])

        #extract the wavelength range for the template. I use it only to know the initial guesses
        line_flux_template_orig = flux_template_broad[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

        #normalise the template
        line_flux_template_norm = line_flux_template_orig #norm_spec(line_wave, line_flux_template_orig, wave_norm, epsilon_norm, line_flux_template_orig)
    
    
        #resudials
        chisqr = 0. #it is a chisquared
        for i in range(len(line_wave)):
            chisqr = chisqr + (line_flux_spec_norm[i]-line_flux_template_norm[i])**2/line_flux_template_norm[i]

        chisquare_fit.append(chisqr) #chi squared vector!

        #test to stop the cycle when reached the minimum chi square, without exploring other values
        min_index = np.argmin(chisquare_fit)
        if j > 0:
            if min_index < j:
                line_flux_template_norm = previous_line_flux_template_norm # if I enter this cycle, I use the previous stored value of the broadened template because its the best
                break

    if (j == 1): # that means that I found the minumum at 0, with no broadening of the template
        print ('Warning: the resolution of the template is greater than the sigma you want to measure. The value obtained will be just a superior limit')

    #finding the minimum chi square and the respective value of sigma
    min_value = np.argmin(chisquare_fit)
    
    sigma_vel = sigma_values[min_value] # this is the holy grail: the most probable velocity dispersion value
    min_residual = chisquare_fit[min_value] # this is the chi square of the best fit


    #***************************** Uncertainties with MonteCarlo simulations ************
    # I perturb the line_flux_temp_fit_norm by adding noise equal to the SNR of my line, then I do the Gaussian fit with the same non-noisy template and see how the sigma that I obtain from the Gaussian fit fluctuates.
    if err_calc == True:
        scale_noise = np.mean(snr_flux)/snr
        number_noisy_cont = 20 #how many syntethics? a lot, but you need to consider also computation time

        sigma_vel_err_tot_ = []

        print ('Calculating the error...')
        for k in range (number_noisy_cont):        
            #generate the noisy line
            noise_array = np.random.standard_normal((len(line_wave),))
            noise_array_scaled = noise_array * scale_noise
            noisy_template = []
            for i in range(len(line_wave)):
                noisy_template.append(line_flux_template_norm[i] + noise_array_scaled[i])
            
            #maximum error on the sigma estimation
            max_error_sigma = int(round(700./snr)) #empirical value. Checked and seems ok
            if max_error_sigma > 70 and k == 0:
                max_error_sigma = 70
                print ('Warning: SNR < 10, large error and long processing time!')
            step_error = 2.
            initial_sigma_err = sigma_vel - max_error_sigma
            
            #if I reach negative values in the interval of possible sigma values:
            if initial_sigma_err < 0:
                initial_sigma_err = 0
            final_sigma_err = sigma_vel + max_error_sigma
            
            number_values_sigma_err = int(round((final_sigma_err-initial_sigma_err)/step_error+1))
            #filling the sigma vector
            sigma_values_err = []
            sigma_value_err = initial_sigma_err
            
            for i in range(1, number_values_sigma_err):
                sigma_values_err.append(sigma_value_err)
                sigma_value_err =  sigma_value_err + step_error
                
            residuals_err = []
            
            # fitting procedure
            for h in range(len(sigma_values_err)):
                #broadening the template
                #sigma_to_broad = 250. # questo dovrà essere un vettore
                flux_template_new = sigma_broad(wavelength, flux_template_original, sigma_values_err[h])

                #extract the wavelength range for the template. I use it only to know the initial guesses
                line_flux_template_new = flux_template_new[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

                #normalize the template
                line_flux_template_new_norm = norm_spec(line_wave, line_flux_template_new, wave_norm, epsilon_norm, line_flux_template_new)

                #resudials
                residual_err = 0. #it is a chisquared
                for i in range(len(line_wave)):
                    residual_err = residual_err + (noisy_template[i]-line_flux_template_new_norm[i])**2/line_flux_template_new_norm[i]
                    
                residuals_err.append(residual_err) #chi squared vector!
                
                min_index = np.argmin(residuals_err)
                if h > 1:
                    if min_index < h:
                        break

            #finding the minimum
            min_value_err = np.argmin(residuals_err)
            sigma_vel_err = sigma_values_err[min_value_err]
            sigma_vel_err_tot_.append(sigma_vel_err)

        #storing the data
        sigma_vel_err_tot = np.asarray(sigma_vel_err_tot_, dtype = float)
        error_sigma_fit = np.std(sigma_vel_err_tot)

        #the total error is the quadratic sum of the error above + the quantization error due to the step of sigma values selected.
        total_error = mt.sqrt(error_sigma_fit**2+step_error**2 +step_sigma**2)
    else:
        total_error = 0.

    ####################################################################################################

    sigma_real_broadened = mt.sqrt(sigma_vel**2+resolution_temp**2) # This is the 
    
    if (sigma_real_broadened == resolution_temp):
        sigma_real = sigma_real_broadened
    elif (sigma_real_broadened < sigma_instrumental):    
        sigma_real = sigma_instrumental
        print ('WARNING: The real velocity dispersion is lower than the instrumental sigma of the spectrum. Do not trust the result!')
    else:
        sigma_real = mt.sqrt(sigma_real_broadened**2 - sigma_instrumental**2)
    
    print ('Resolution template in A (FWHM): ', resolution_lambda_fwhm_temp*10)
    print ('Resolution spectrum in A (FWHM): ', resolution_lambda_fwhm_spec*10)
    print ('Resolution sigma template (km/s): ', resolution_temp)
    print ('Resolution sigma spectrum: (km/s)', sigma_instrumental)
    print ('Best sigma gaussian broadening: ', sigma_vel)
    print ('Template best real total broadening: ', sigma_real_broadened , 'km/s')
    print ('Sigma spectrum = sqrt(best broadening^2- resolution sigma spectrum^2): ', sigma_real , 'km/s')

    return sigma_real, total_error, min_residual, line_wave, line_flux_spec_norm, line_flux_template_norm, sigma_instrumental


#*************************************************************************************************
#25 Show sampling and identify linear or log spectrum (hopefully!)
def show_sampling(wavelength):
    step1 = wavelength[1]-wavelength[0]
    npoints = len(wavelength)
    step2 = wavelength[npoints-1]-wavelength[npoints-2]
    epsilon = 1e-4
    if (abs(step1-step2) >= epsilon):
        linear_step = False
        print ('Not linear spectrum')
        return step1, linear_step
    else:
        linear_step = True
        print ('Linear spectrum')
        return step1, linear_step
    
#*************************************************************************************************
#26 Show the SNR of a selected window
def show_snr(wavelength, flux, wave_snr, epsilon_wave_snr):
    step = wavelength[1]-wavelength[0]
    
    #check for constant step and if not resample
    step2 = wavelength[len(wavelength)-1]- wavelength[len(wavelength)-2]
    epsilon = 1e-4
    if abs(step-step2) > epsilon:
        wavelength, flux, npoint_resampled = resample(wavelength, flux, step)
        print ('Spectrum resampled to a linear step')
                   
    initial_wave = wave_snr - epsilon_wave_snr
    final_wave = wave_snr + epsilon_wave_snr
    
    #extracting the wavelength and the flux within the selected band
    wave_snr = wavelength[(wavelength >= initial_wave) & (wavelength <= final_wave)]
    flux_snr = flux[(wavelength >= initial_wave) & (wavelength <= final_wave)]
    mean_flux = np.mean(flux_snr)
    snr_pix = mean_flux/np.std(flux_snr)
    
    snr_ang = snr_pix * mt.sqrt(1/(step*10)) # supposing all the units in nm!
    
    return snr_pix, snr_ang
    

#*************************************************************************************************
#27 Show the header of a fits file
def show_hdr (spec_name):
    fmt_spec2 = '.fits' in spec_name
    fmt_spec1 = '.txt' in spec_name or '.dat' in spec_name
    
    #if I have an ASCII file, well, no good
    if (fmt_spec1 == True and fmt_spec2 == False or (fmt_spec1 == True and fmt_spec2 == True)):
        header ='ASCII files do not have header!'
        
    #if I have a fits file, reading the header stored in the hdu[0]
    else:
        hdu = fits.open(spec_name)
        hdr = hdu[0].header
        header = (repr(hdr))
    return header


#*************************************************************************************************
#28 Convert a spectrum to ASCII or binary fits file
def convert_spec(wavelength, flux, spec_name, type_spec_to_convert):
    if type_spec_to_convert == 'ASCII':
        filename = spec_name + '.txt'
        np.savetxt(filename, np.column_stack([wavelength, flux]), header="wavelength \t flux")
    
    elif type_spec_to_convert == 'FITS':
        filename = spec_name +'.fits'
        t = Table([wavelength, flux], names=('wavelength', 'flux'))
        t.write(filename, format='fits')
        
    return 0
        
        
#*************************************************************************************************
#29 Simple box window moving average
def mov_avg(flux, window_size):
    avg_flux = np.convolve(flux, np.ones((window_size,))/window_size, mode='same') # for now I decided to not handle the edges with the keyword 'same'
    return avg_flux


#*************************************************************************************************
#31 Heliocentric calculation and correction on the spectrum
def helio_corr(wavelength, flux, epoch, where, ra_obj, dec_obj):
    location = EarthLocation.of_site(where)
    sc = SkyCoord(ra=ra_obj*u.deg, dec=dec_obj*u.deg)
    heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(epoch), location = location)
    heliocorr = heliocorr.to(u.km/u.s)  
    location_list = EarthLocation.get_site_names()
    
    heliocorr = (str(heliocorr))
    
    #extract the value from heliocorr
    helio_number = 0.
    for heliocorr in heliocorr.split():
        try:
            helio_number = (float(heliocorr))
        except ValueError:
            pass
    
    #correct_the_spectrum
    new_wavelength, new_flux = dopcor(wavelength, flux, helio_number)
    
    return helio_number, new_wavelength, new_flux


#*************************************************************************************************
#32 Equivalenth width measurements
def ew_measurement(wavelength, flux, index_file, is_usr_idx, want_plot, verbose, calculate_error, save_plot, spec_name, normalize_spec):
    new_step = wavelength[1]-wavelength[0]

    #reading the index file
    if is_usr_idx == False: #if index_file is a file containing the index definitions
        id_array, index = read_idx(index_file)
        num_indices = len(id_array)
        
    else: # if index_file is a numpy array already containing the definition of a usr_idx
        id_array = 'usr_idx'
        num_indices = 1
        index = index_file
    
    if verbose == True:
        print ('Number of indices to measure = ', num_indices)

    if verbose == True:
        if is_usr_idx == False:
            index_transpose = np.transpose(index) #only for the following visualization
            for i in range(num_indices):
                print (id_array[i], index_transpose[i,:]) #printing the indices names and their wavelengths
        else: 
            print(id_array)
            print(index)
        
    #defining the arrays containing the data
    ew_array = np.zeros(num_indices)
    ew_array_mag = np.zeros(num_indices)
    snr_ew_array = np.zeros(num_indices)
    err_array = np.zeros(num_indices)
    err_array_mag = np.zeros(num_indices)
    wave_limits_spec = np.array([wavelength[0], wavelength[len(wavelength)-1]])
    #************************************* Cycling all over the indices **********************************
  

        #if I have a list of indices
    if is_usr_idx == False and normalize_spec == True:
        for t in range(num_indices): 


            #let's normalise the spectrum to one of the band index. Easier and with no difference with respect to before
            lambda_ref_norm = index[4,t]
     
            #check if the index interval is in the wavelength range of the spectrum. Need to insert also the normalization lambda because it can be outside the range!
            min_idx = np.min(index[:,t])
            max_idx = np.max(index[:, t])
            if (min_idx < wave_limits_spec[0] or max_idx > wave_limits_spec[1] or lambda_ref_norm > np.max(wave_limits_spec) or lambda_ref_norm < np.min(wave_limits_spec)):
                print ('Warning: index not in the spectrum limits. Skipping')
                continue
            
            # 1)Normalizing the spectra
            epsilon_wave = new_step*10. #just an epsilon value to average the flux for the normalization value
            if normalize_spec == True:
                norm_flux = norm_spec(wavelength, flux, lambda_ref_norm, epsilon_wave, flux)
            elif normalize_spec == False:
                norm_flux = flux
                
            # 2) Extract the spectral region around the index (+/- 5 nm) to speed up the process
            index_wave, index_flux = extract_index(wavelength, norm_flux, index[:,t])

            #3) Extract the pseudo continuum
            interp_cont_flux, line_flux, interp_lambda, flux_ref_cont, lambda_ref_cont, central_left_lambda, central_right_lambda, avg_left_flux, avg_right_flux = idx_cont(index[:,t], index_wave, index_flux)
            
            #4) Determining the EW of the index
            ew = eq_width(interp_cont_flux, line_flux, new_step)
            
            #4) BIS Determining the EW of the index in mag
            ew_mag = eq_width_mag(interp_cont_flux, line_flux, new_step, index[4,t],index[5,t])
            
            #5) Calculate the error via MonteCarlo simulation, if you want
            if calculate_error == True:
                err, sigma_cont = ew_err(index[:,t], index_wave, index_flux, new_step, interp_cont_flux, interp_lambda, flux_ref_cont, lambda_ref_cont)
            else:
                err = 0.
                
            #trasform in A
            ew = ew*10
            err = err*10
            
            #calculating the errors in magnitudes
            err_mag= 0.434*abs(err/ew)

            #fill the arrays
            ew_array[t] = ew
            err_array[t] = err
            
            ew_array_mag[t] = ew_mag
            err_array_mag[t] = err_mag
            
            # Calculate the snr
            if calculate_error == True:
                snr_per_pix = (avg_left_flux + avg_right_flux)/(2*sigma_cont)
                pix_per_a = 1/(new_step*10) #new_step*10 because the step is given in nm
                snr_per_a = snr_per_pix*mt.sqrt(pix_per_a)
                
                #fill the array
                snr_ew_array[t] = snr_per_pix
            else:
                snr_ew_array[t] = 0

            # doing the plots for all the indices
            if(want_plot == True and save_plot == False):
                band = []
                for i in range (len(wavelength)):
                    if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                        band.append(wavelength[i])
                interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
                interp_flux_bands = (interpfunc(band))
                
                #set the y limits for the plots
                ylim_low = np.mean(interp_flux_bands)-0.6
                ylim_high = np.mean(interp_flux_bands)+0.8
            
                yeps = 0.1
                ew_string = str(round(ew,2))
                err_string = str(round(err,2))
                snr_ew_string = str(round(snr_per_pix,0))
                
                if (index[3, t] < index[5, t]):
                    plt.title('EW ' + spec_name + ' ' + id_array[t] + ' ' + ew_string + '+/-' + err_string + ' A.  SNR =' + snr_ew_string)
                    plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                    plt.xlim(index[0, t]-5., index[5, t]+5.)
                    plt.ylim (ylim_low, ylim_high)  #(0.4, 1.8)
                    plt.xlabel('Wavelength nm', fontsize = 9)
                    plt.ylabel('Flux', fontsize = 9)
                    plt.tick_params(axis = 'both', labelsize = 9)

                else:
                    plt.title('EW ' + spec_name + ' ' + id_array[t] + ' ' + ew_string + '+/-' + err_string + ' A.  SNR =' + snr_ew_string)
                    plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                    plt.xlim(index[0, t]-5., index[3, t]+5.)
                    plt.ylim(ylim_low, ylim_high)
                    plt.xlabel('Wavelength nm', fontsize = 9)
                    plt.ylabel('Flux', fontsize = 9)
                    plt.tick_params(axis = 'both', labelsize = 9)
                
                #pseudo continuum
                plt.plot(band, interp_flux_bands, linewidth = 0.5, color = 'black')
                plt.plot(interp_lambda, interp_cont_flux, linewidth = 1, color = 'blue')
                
                #polygons
                x_polygon_bband = [index[0, t], index[1,t], index[1,t], index[0,t], index[0,t]]
                y_polygon_bband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                plt.fill(x_polygon_bband,y_polygon_bband, 'blue')
                
                x_polygon_rband = [index[2,t], index[3,t], index[3,t], index[2,t], index[2,t]]
                y_polygon_rband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                plt.fill(x_polygon_rband,y_polygon_rband, 'red')
                
                x_polygon_line = [index[4,t], index[5,t], index[5,t], index[4,t], index[4,t]]
                y_polygon_line = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                plt.fill(x_polygon_line,y_polygon_line, 'gray')
                plt.show()
            
            if(want_plot == False and save_plot == True):
                band = []
                for i in range (len(wavelength)):
                    if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                        band.append(wavelength[i])
                interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
                interp_flux_bands = (interpfunc(band))
                
                #set the y limits for the plots
                ylim_low = np.mean(interp_flux_bands)-0.6
                ylim_high = np.mean(interp_flux_bands)+0.8
                
                yeps = 0.1
                ew_string = str(round(ew,3))
                err_string = str(round(err,3))
                
                #how mani subplot?
                fig, axs = plt.subplots(3)
                fig.set_size_inches(7.5, 10.)
            
                #space between subplots
                fig.subplots_adjust(hspace=0.3, top = 0.95) #wspace=0.4)
            
                #first plot
                axs[0].plot(wavelength, flux, linewidth = 0.5, color = 'black')
                #axs[0].set_title(id_array[t])
                axs[0].set_xlabel('Wavelength nm', fontsize = 9)
                axs[0].set_ylabel('Flux', fontsize = 9)
                axs[0].tick_params(axis = 'both', labelsize = 9)
                #axs[0].plt.title('Scores by group and gender')
                
                #second plot
                axs[1].plot(wavelength, norm_flux, label = 'Original', linewidth=0.5, color = 'black' )
                axs[1].plot(index_wave, index_flux, label = 'Index', color = 'green')
                axs[1].set_xlabel('Wavelength nm', fontsize = 9)
                axs[1].set_ylabel('Flux', fontsize = 9)
                axs[1].set_xlim(index[0,t]-40., index[3,t]+40.)
                axs[1].set_ylim(ylim_low, ylim_high)
                axs[1].legend(fontsize = 6)
                axs[1].tick_params(axis = 'both', labelsize = 9)
                
                #third plot
                yeps = 0.1
                if (index[3,t] < index[5,t]):
                    axs[2].plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                    axs[2].set_xlim(index[0,t]-5., index[5,t]+5.)
                    axs[2].set_ylim(ylim_low, ylim_high)
                    axs[2].set_xlabel('Wavelength nm', fontsize = 9)
                    axs[2].set_ylabel('Flux', fontsize = 9)
                    axs[2].tick_params(axis = 'both', labelsize = 9)

                else:
                    axs[2].plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                    axs[2].set_xlim(index[0,t]-5., index[3,t]+5.)
                    axs[2].set_ylim(ylim_low, ylim_high)
                    axs[2].set_xlabel('Wavelength nm', fontsize = 9)
                    axs[2].set_ylabel('Flux', fontsize = 9)
                    axs[2].tick_params(axis = 'both', labelsize = 9)
                
                #pseudo continuum
                axs[2].plot(band, interp_flux_bands, linewidth = 0.5, color = 'black')
                axs[2].plot(interp_lambda, interp_cont_flux, linewidth = 1, color = 'blue')
                
                #polygons
                x_polygon_bband = [index[0,t], index[1,t], index[1,t], index[0,t], index[0,t]]
                y_polygon_bband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
                axs[2].fill(x_polygon_bband,y_polygon_bband, 'blue')
                
                x_polygon_rband = [index[2,t], index[3,t], index[3,t], index[2,t], index[2,t]]
                y_polygon_rband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
                axs[2].fill(x_polygon_rband,y_polygon_rband, 'red')
                
                x_polygon_line = [index[4,t], index[5,t], index[5,t], index[4,t], index[4,t]]
                y_polygon_line = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
                axs[2].fill(x_polygon_line,y_polygon_line, 'gray')
                
                
                result_plot_dir = 'results/plots'
                #pathlib.Path(result_plot_dir).mkdir(parents=True, exist_ok=True)
                os.makedirs(result_plot_dir, exist_ok=True)
                
                plt.savefig(result_plot_dir + '/'+spec_name + '_' + id_array[t]+'.eps', format='eps', dpi=100)
                plt.close() # se non chiudo il plot va in overflow di memoria

                # FINE NUOVO #
            
        return id_array, ew_array, err_array, snr_ew_array, ew_array_mag, err_array_mag

        
        # if I have just a single index
    else: 
        id_array = 'usr_idx'
        if normalize_spec == True:

            #normalization wavelength
            lambda_ref_norm = index[4]

            # 1)Normalising the spectra
            epsilon_wave = new_step*10. #just an epsilon value to average the flux for the normalization value
            norm_flux = norm_spec(wavelength, flux, lambda_ref_norm, epsilon_wave, flux)
            
        if normalize_spec == False:
            norm_flux = flux


        # 2) Extract the spectral region around the index (+/- 5 nm) to speed up the process
        index_wave, index_flux = extract_index(wavelength, norm_flux, index)

        #3) Extract the pseudo continuum
        interp_cont_flux, line_flux, interp_lambda, flux_ref_cont, lambda_ref_cont, central_left_lambda, central_right_lambda, avg_left_flux, avg_right_flux = idx_cont(index, index_wave, index_flux)
            
        #4) Determining the EW of the index
        ew = eq_width(interp_cont_flux, line_flux, new_step)
        
                
        #4) BIS Determining the EW of the index in MAG
        ew_mag = eq_width_mag(interp_cont_flux, line_flux, new_step, index[4],index[5])
        
        #5) Calculate the error via MonteCarlo simulation
        if calculate_error == True:
            err, sigma_cont = ew_err(index, index_wave, index_flux, new_step, interp_cont_flux, interp_lambda, flux_ref_cont, lambda_ref_cont)
        else:
            err = 0.
            
        #fill the vectors
        ew = ew*10
        err = err*10
        err_mag= 0.434*abs(err/ew)


        
        # Calculate the snr
        if calculate_error == True:
            snr_per_pix = (avg_left_flux + avg_right_flux)/(2*sigma_cont)
            pix_per_a = 1/(new_step*10) #new_step*10 because the step is given in nm
            snr_per_a = snr_per_pix*mt.sqrt(pix_per_a)
        else: 
            snr_per_pix = 0
            snr_per_a = 0
            
        # doing the plots for the single index and only if I want
        if(want_plot == True and save_plot == False):
            band = []
            for i in range (len(wavelength)):
                if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                    band.append(wavelength[i])
            interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
            interp_flux_bands = (interpfunc(band))
                
            #set the y limits for the plots
            ylim_low = np.mean(interp_flux_bands)-0.6
            ylim_high = np.mean(interp_flux_bands)+0.8
                
            yeps = 0.1
            ew_string = str(round(ew,3))
            err_string = str(round(err,3))
            snr_ew_string = str(round(snr_per_pix,0))
                
            if (index[3] < index[5]):
                plt.title('EW usr index ' + spec_name + ' ' + ew_string + '+/-' + err_string + ' A.  SNR =' + snr_ew_string)
                plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                plt.xlim(index[0]-5., index[5]+5.)
                plt.ylim(ylim_low, ylim_high)
                plt.xlabel('Wavelength nm', fontsize = 9)
                plt.ylabel('Flux', fontsize = 9)
                plt.tick_params(axis = 'both', labelsize = 9)
                    
            else:
                plt.title('EW usr index ' + spec_name + ' ' + ew_string + '+/-' + err_string + ' A.  SNR =' + snr_ew_string)
                plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                plt.xlim(index[0]-5., index[3]+5.)
                plt.ylim(ylim_low, ylim_high)
                plt.xlabel('Wavelength nm', fontsize = 9)
                plt.ylabel('Flux', fontsize = 9)
                plt.tick_params(axis = 'both', labelsize = 9)
            
            #pseudo continuum
            plt.plot(band, interp_flux_bands, linewidth = 0.5, color = 'black')
            plt.plot(interp_lambda, interp_cont_flux, linewidth = 1, color = 'blue')
                
            #polygons
            x_polygon_bband = [index[0], index[1], index[1], index[0], index[0]]
            y_polygon_bband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_bband,y_polygon_bband, 'blue')
                
            x_polygon_rband = [index[2], index[3], index[3], index[2], index[2]]
            y_polygon_rband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_rband,y_polygon_rband, 'red')
                
            x_polygon_line = [index[4], index[5], index[5], index[4], index[4]]
            y_polygon_line = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_line,y_polygon_line, 'gray')
            plt.show()

            
            
        if(want_plot == False and save_plot == True):
            band = []
            for i in range (len(wavelength)):
                if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                    band.append(wavelength[i])
            interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
            interp_flux_bands = (interpfunc(band))
                
            yeps = 0.1
            ew_string = str(round(ew,3))
            err_string = str(round(err,3))
                
            #how mani subplot?
            fig, axs = plt.subplots(3)
            fig.set_size_inches(7.5, 10.)
            
            #space between subplots
            fig.subplots_adjust(hspace=0.3, top = 0.95) #wspace=0.4)
            
            #first plot
            axs[0].plot(wavelength, flux, linewidth = 0.5, color = 'black')
            #axs[0].set_title(id_array[t])
            axs[0].set_xlabel('Wavelength nm', fontsize = 9)
            axs[0].set_ylabel('Flux', fontsize = 9)
            axs[0].tick_params(axis = 'both', labelsize = 9)
            #axs[0].plt.title('Scores by group and gender')
                
            #second plot
            axs[1].plot(wavelength, norm_flux, label = 'Original', linewidth=0.5, color = 'black' )
            axs[1].plot(index_wave, index_flux, label = 'Index', color = 'green')
            axs[1].set_xlabel('Wavelength nm', fontsize = 9)
            axs[1].set_ylabel('Flux', fontsize = 9)
            axs[1].set_xlim(index[0]-40., index[3]+40.)
            axs[1].set_ylim(0.3, 1.8)
            axs[1].legend(fontsize = 6)
            axs[1].tick_params(axis = 'both', labelsize = 9)
                
            #third plot
            yeps = 0.1
            if (index[3] < index[5]):
                axs[2].plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                axs[2].set_xlim(index[0]-5., index[5]+5.)
                axs[2].set_ylim(0.4, 1.8)
                axs[2].set_xlabel('Wavelength nm', fontsize = 9)
                axs[2].set_ylabel('Flux', fontsize = 9)
                axs[2].tick_params(axis = 'both', labelsize = 9)

            else:
                axs[2].plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                axs[2].set_xlim(index[0]-5., index[3]+5.)
                axs[2].set_ylim(0.4, 1.8)
                axs[2].set_xlabel('Wavelength nm', fontsize = 9)
                axs[2].set_ylabel('Flux', fontsize = 9)
                axs[2].tick_params(axis = 'both', labelsize = 9)
                
            #pseudo continuum
            axs[2].plot(band, interp_flux_bands, linewidth = 0.5, color = 'black')
            axs[2].plot(interp_lambda, interp_cont_flux, linewidth = 1, color = 'blue')
                
            #polygons
            x_polygon_bband = [index[0], index[1], index[1], index[0], index[0]]
            y_polygon_bband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
            axs[2].fill(x_polygon_bband,y_polygon_bband, 'blue')
                
            x_polygon_rband = [index[2], index[3], index[3], index[2], index[2]]
            y_polygon_rband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
            axs[2].fill(x_polygon_rband,y_polygon_rband, 'red')
                
            x_polygon_line = [index[4], index[5], index[5], index[4], index[4]]
            y_polygon_line = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
            axs[2].fill(x_polygon_line,y_polygon_line, 'gray')
                
                
            result_plot_dir = 'results/plots'
            #pathlib.Path(result_plot_dir).mkdir(parents=True, exist_ok=True)
            os.makedirs(result_plot_dir, exist_ok=True)
            
            plt.savefig(result_plot_dir + '/'+spec_name + '_' + id_array +'.eps', format='eps', dpi=100)
            plt.close()
            
    return id_array, ew, err, snr_per_pix, ew_mag, err_mag
    
# 33) Calculating the velocity dispersion coefficients to correct the EW due to velocity dispersion broadening
def sigma_coeff (spectra_file, index_file, lambda_units, is_usr_idx, want_plot, smooth_value, save_plot):
    
    spectra_list = np.loadtxt(spectra_file, dtype = 'str')
    number_spec = len(spectra_list)
    print ('Number of spectra', number_spec)
    #array of sigmas to be tested
    sigma_array = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400])
    sigma_values = len (sigma_array)

    
    
    if is_usr_idx == True:
        print ('User index')
        ew_sigma = np.zeros((sigma_values, number_spec))
        ew_mean = np.zeros(sigma_values)
        ew_std = np.zeros(sigma_values)
        idx_array = 'usr_idx'
        # for every spectrum in my list:
        stop_condition = 0
        for i in range(number_spec):
            if not sg.OneLineProgressMeter('Task progress', i+1, number_spec,  'single', 'Processing spectra:', orientation='h',button_color=('white','red')):
                print ('***CANCELLED***')
                print ('')
                stop_condition = 1
                break
            
            #reading the spectrum
            wavelength, flux, original_step, obj_name = read_spec(spectra_list[i], lambda_units)
            
            #resample if step not constant
            step1 = wavelength[1]-wavelength[0]
            step2 = wavelength[len(wavelength)-1]- wavelength[len(wavelength)-2]
            epsilon = 1e-4
            if abs(step1-step2) > epsilon:
                wavelength, flux, npoint_resampled = resample(wavelength, flux, original_step)
                print('Spectrum resampled to a linear step')
                
            if smooth_value != 0:
                flux = sigma_broad(wavelength, flux, smooth_value)
            
            #for every sigma value:
            for j in range(sigma_values):
                
                # broadening the spectra
                flux_broadened = sigma_broad(wavelength, flux, sigma_array[j])
                
                #measuring the index/indices, for one index
                if j == 0:
                    id_array, ew_orig, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, index_file, True, False, False, False, False, 'fake_name', True)
                    ew_sigma[j,i] = 0.
                    
                #measuring the EW
                else:
                    id_array, ew, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, index_file, True, False, False, False, False, 'fake_name', True)
                    
                    #storing the ew in the array containing row = sigma values; columns = spectra number
                    ew_sigma[j,i] = (ew-ew_orig)/ew_orig
            
            if want_plot == True:
                if i == 0:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                    fig.suptitle('Broadening coefficients for ' + idx_array)
                    
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                else:
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)
            
            #NUOVO
            if save_plot == True: 
                if i == 0:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                    fig.suptitle('Broadening coefficients for ' + idx_array)
                    
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                else:
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)
            
        #filling the vectors with mean and std
        for h in range (sigma_values):
            ew_mean[h] = np.mean(ew_sigma[h,:])
            ew_std[h] = np.std(ew_sigma[h,:])
        
        #plotting
        if want_plot == True:
            ax1.errorbar(sigma_array, ew_mean, yerr = ew_std, color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
            ax1.set_xlabel('Broadening (km/s)')
            ax1.set_ylabel('Relative variation')
            ax1.legend(fontsize = 10)
            
            ax2.plot(sigma_array, ew_std, color = 'red')
            ax2.set_xlabel('Broadening (km/s)')
            ax2.set_ylabel('Error (1 sigma)')
            #plt.title(idx_array)
            plt.show()
            
        #NUOVO
        if save_plot == True:
            ax1.errorbar(sigma_array, ew_mean, yerr = ew_std, color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
            ax1.set_xlabel('Broadening (km/s)')
            ax1.set_ylabel('Relative variation')
            ax1.legend(fontsize = 10)
            
            ax2.plot(sigma_array, ew_std, color = 'red')
            ax2.set_xlabel('Broadening (km/s)')
            ax2.set_ylabel('Error (1 sigma)')
            
            result_plot_dir = 'results/plots'
            os.makedirs(result_plot_dir, exist_ok=True)
            
            plt.savefig(result_plot_dir + '/'+ 'sigma_coeff_' + id_array +'.eps', format='eps', dpi=100)
            plt.close() # remember to always close the plots!

        
        ew_coeff = np.polyfit(sigma_array, ew_mean, 3)
        err_coeff = np.polyfit(sigma_array, ew_std, 3)
        
        return idx_array, ew_coeff, err_coeff, ew_mean, ew_std, stop_condition
        
        
        
    #if I have an index file
    if is_usr_idx == False:
        idx_array, indices = read_idx(index_file)
        num_indices = len(idx_array)
        ew_sigma = np.zeros((sigma_values, number_spec))
        ew_mean = np.zeros((sigma_values, num_indices))
        ew_std = np.zeros((sigma_values, num_indices))
        num_coeff = 4
        ew_coeff_array = np.zeros((num_coeff, num_indices))
        err_coeff_array = np.zeros((num_coeff, num_indices))
        print ('Number of indices: ', num_indices)
        
        # for every index in my list:
        stop_condition = 0
        for k in range (num_indices):
            
            if not sg.OneLineProgressMeter('Task progress', k+1, num_indices,  'single', 'Processing index:', orientation='h',button_color=('white','red')):
                print ('***CANCELLED***')
                print ('')
                stop_condition = 1
                break

            # for every spectrum in my list:
            for i in range(number_spec):
                
                #reading the spectrum
                wavelength, flux, original_step, obj_name = read_spec(spectra_list[i], lambda_units)
                
                #resample if step not constant
                step1 = wavelength[1]-wavelength[0]
                step2 = wavelength[len(wavelength)-1]- wavelength[len(wavelength)-2]
                epsilon = 1e-4
                if abs(step1-step2) > epsilon:
                    wavelength, flux, npoint_resampled = resample(wavelength, flux, original_step)
                    print('Spectrum resampled to a linear step')
                
                if smooth_value != 0:
                    flux = sigma_broad(wavelength, flux, smooth_value)
                
                #for every sigma value:
                for j in range(sigma_values):
                    
                    # broadening the spectra
                    flux_broadened = sigma_broad(wavelength, flux, sigma_array[j])
                    
                    #measuring the index/indices, for one index
                    if j == 0:
                        id_array, ew_orig, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, indices[:,k], True, False, False, False, False, 'fake_name', True)
                        ew_sigma[j,i] = 0.
                    #measuring the EW
                    else:
                        id_array, ew, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, indices[:,k], True, False, False, False, False, 'fake_name', True)
                        
                        #storing the ew in the array containing row = sigma values; columns = spectra number
                        ew_sigma[j,i] = (ew-ew_orig)/ew_orig
                
                # plotting the single values
                if want_plot == True:
                    if i == 0:
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                        fig.suptitle('Broadening coefficients for ' + idx_array[k])
                        
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                    else:
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)
                        
                #NUOVO
                if save_plot == True:                       
                    if i == 0:
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                        fig.suptitle('Broadening coefficients for ' + idx_array[k])
                        
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                    else:
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)
                        
            #filling the vectors with mean and std
            for h in range (sigma_values):
                ew_mean[h,k] = np.mean(ew_sigma[h,:])
                ew_std[h,k] = np.std(ew_sigma[h,:])
            
            #plotting
            if want_plot == True:
                ax1.errorbar(sigma_array,ew_mean[:,k], yerr = ew_std[:,k], color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
                ax1.set_xlabel('Broadening (km/s)')
                ax1.set_ylabel('Relative variation')
                ax1.legend(fontsize = 10)
                
                ax2.plot(sigma_array, ew_std[:,k], color = 'red')
                ax2.set_xlabel('Broadening (km/s)')
                ax2.set_ylabel('Relative error')
                plt.show()
                
               
            #NUOVO
            if save_plot == True:
                ax1.errorbar(sigma_array,ew_mean[:,k], yerr = ew_std[:,k], color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
                ax1.set_xlabel('Broadening (km/s)')
                ax1.set_ylabel('Relative variation')
                ax1.legend(fontsize = 10)
                
                ax2.plot(sigma_array, ew_std[:,k], color = 'red')
                ax2.set_xlabel('Broadening (km/s)')
                ax2.set_ylabel('Relative error')
            
                result_plot_dir = 'results/plots'
                os.makedirs(result_plot_dir, exist_ok=True)
                
                plt.savefig(result_plot_dir + '/'+ 'sigma_coeff_' + idx_array[k] +'.eps', format='eps', dpi=100)
                plt.close() #Remember to always close the plots!

            
            ew_coeff = np.polyfit(sigma_array, ew_mean[:,k], 3)
            err_coeff = np.polyfit(sigma_array, ew_std[:,k], 3)
            
            ew_coeff_array[:,k] = ew_coeff
            err_coeff_array[:,k] = err_coeff

        return idx_array, ew_coeff_array, err_coeff_array, ew_mean, ew_std, stop_condition
            
            
# 34) Calculating the resolution R
def Gauss(x, y0, x0, a, sigma):
    return y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Gaussian model with slope
def Gauss_slope(x, y0, x0, a, sigma, m, c):
    y = y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))+m*x + c
    return y

#Multiple gaussians
def multiple_gauss(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 6):
        y0 = params[i]
        x0 = params[i+1]
        a = params[i+2]
        sigma = params[i+3]
        m = params[i+4]
        c = params[i+5]
        y = y + y0 + a * np.exp(-(x - x0)**2/(2*sigma**2)) +m*x+c
    return y

#Measure the resolution from an emission (sky) line
def resolution (wavelength, flux, wave1, wave2):
    step = wavelength[1]-wavelength[0]
    
    #extract the line flux and wavelength arrays
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]  
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    #normalize the spectrum
    wave_norm = line_wave[10] # guessing for now. fix it later
    epsilon_norm = step*10
    line_flux_spec_norm = norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)


    ############### FITTING PROCEDURE ###########
    #initial data guess for the spectra; looking for EMISSION patterns
    max_flux_spec = np.argmax(line_flux_spec_norm) # searching for a peak
    mean_spec = line_wave[max_flux_spec] # lambra corresponding to the flux peak found
    sigma_spec = 0.1
    offset = 1. #wave_norm # offset sull'asse y, somaro! Se lo spettro è già normalizzato è = 1!
        
    #fitting to the spectra
    popt_spec, pconv_spec = curve_fit(Gauss, line_wave, line_flux_spec_norm, p0=[offset, mean_spec, max_flux_spec, sigma_spec])
    
    #assigning values of the gaussian to new variables
    sigma_fit_spec = popt_spec[3]
    mean_fit_spec = popt_spec[2]
   
    #assigning the fit value
    line_flux_spec_fit = Gauss(line_wave, *popt_spec)

    #calculating the resolution
    resolution_fwhm = sigma_fit_spec*2.35
    resolution_R = int(line_wave[0]/resolution_fwhm)
   
    print ('Resolution in A (FWHM): ', round(resolution_fwhm*10,2))
    return resolution_R, line_wave, line_flux_spec_norm, line_flux_spec_fit
    
      
# 35) Convert flux from Jansky to f_lambda or f_nu
def convert_flux(wavelength, flux, spec_name, type_to_convert, lambda_units):
    flux_points = len(flux)
    converted_flux = np.zeros(flux_points)
 
    #re-convert to the original wavelength, so the flux units are consistent!
    if lambda_units == 'mu':
        wavelength = wavelength/1000.
    if lambda_units == 'A':
        wavelength == wavelength/10.
        
        
    if type_to_convert == 'to_flambda':
        conversion_factor = 2.999792458e12
        
        #converting
        for i in range (flux_points):
            converted_flux[i] = flux[i]*conversion_factor/(wavelength[i])**2

    elif type_to_convert == 'to_fnu':
        conversion_factor = 1e26
        #wave_mu = wavelength*1000
        for i in range (flux_points):
            converted_flux[i] = flux[i]*conversion_factor

    return converted_flux


#fitting line
def line_fitting (wavelength, flux, wave_interval, guess_param):
    step = wavelength[1]-wavelength[0]
    wave1 = min(wave_interval)
    wave2 = max(wave_interval)
    
    #isolating the region of interest
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]  
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    #normalize the spectrum
    wave_norm = line_wave[10] # guessing for now. fix it later
    epsilon_norm = step*10
    line_flux_spec_norm = norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #fitting to the spectra
    popt_spec, pconv_spec = curve_fit(Gauss_slope, line_wave, line_flux_spec_norm, p0=guess_param)

    fit = Gauss_slope(line_wave, *popt_spec)
    return line_wave, line_flux_spec_norm, fit, popt_spec



#fitting threee gaussians to the CaT lines
def cat_fitting (wavelength, flux):
    step = wavelength[1]-wavelength[0]
    
    wave1 = 844 #845
    wave2 = 872 #870 
    
    #extract the line flux and wavelength arrays
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]  
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    #normalize the spectrum
    wave_norm = line_wave[10] # guessing for now. fix it later
    epsilon_norm = step*10
    line_flux_spec_norm = norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #initial guesses
    y0 = 1
    x0 = 850
    a = -0.8
    sigma = 0.1
    m = 0.1
    c = 1

    guess = [y0,x0,a,sigma,m,c]
    for i in range(3):
        if i == 0:
            guess = guess
        if i == 1:
            guess+= [y0,854,-0.6, 0.2, m, c]
        if i == 2:
            guess+= [y0,866,-0.6, 0.2, m, c]

    #fitting to the spectra
    popt_spec, pconv_spec = curve_fit(multiple_gauss, line_wave, line_flux_spec_norm, p0=guess)
    fit = multiple_gauss(line_wave, *popt_spec)

    return line_wave, line_flux_spec_norm, fit, popt_spec
        

#kinematics with ppxf and EMILES SSP models
def ppxf_kinematics(wavelength, flux, wave1, wave2, FWHM_gal, is_resolution_gal_constant, R, redshift_guess, sigma_guess, stellar_library, additive_degree):
    ppxf_dir = path.dirname(path.realpath(lib.__file__))
    #converting the wavelengths to angstrom because this is what ppxf wants
    wavelength = wavelength*10
    wave1 = wave1*10
    wave2 = wave2*10
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]  
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]
    flux = line_flux_spec
    wavelength = line_wave
    
    #select the range in wavelength
    lamRange1 = [wave1,wave2]

    redshift_0 = 0                  # Ignore cosmological redshift for local galaxies
    redshift = redshift_guess               # Initial redshift estimate of the galaxy

    #log rebin
    galaxy, ln_lam1, velscale = util.log_rebin(lamRange1, flux)


    galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues
    noise = np.full_like(galaxy, 0.0163)           # Assume constant noise per pixel here

    sps_name = stellar_library

    lam_range_temp = [lamRange1[0]/1.02, lamRange1[1]*1.02]

    # Read SPS models file from my GitHub if not already in the ppxf package dir.
    # The SPS model files are also available here https://github.com/micappe/ppxf_data
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = path.join(ppxf_dir, 'sps_models', basename)
    if not path.isfile(filename):
        url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
        request.urlretrieve(url, filename)


    if is_resolution_gal_constant == True:
        sps = lib.sps_lib(filename, velscale, FWHM_gal, wave_range=lam_range_temp)
    if is_resolution_gal_constant == False:
        FWHM_gal = wavelength/R
        FWHM_gal = {"lam": ln_lam1, "fwhm": FWHM_gal}
        sps = lib.sps_lib(filename, velscale, FWHM_gal, wave_range=lam_range_temp)

    # Compute a mask for gas emission lines
    goodPixels = util.determine_goodpixels(ln_lam1, lam_range_temp, redshift)


    c = 299792.458
    vel = c*np.log(1 + redshift)   # eq.(8) of Cappellari (2017, MNRAS)
    start = [vel, sigma_guess]  # (km/s), starting guess for [V, sigma]
    t = clock()

    try:
        pp = ppxf(sps.templates, galaxy, noise, velscale, start,
                goodpixels=goodPixels, plot=True, moments=4, lam=np.exp(ln_lam1),
                lam_temp=sps.lam_temp, degree=additive_degree)


        # The updated best-fitting redshift is given by the following
        # lines (using equations 5 of Cappellari 2022, arXiv, C22)
        errors = pp.error*np.sqrt(pp.chi2)  # Assume the fit is good chi2/DOF=1
        redshift_fit = (1 + redshift_0)*np.exp(pp.sol[0]/c) - 1  # eq. (5c) C22
        redshift_err = (1 + redshift_fit)*errors[0]/c            # eq. (5d) C22

        print("Formal errors:")
        print("     dV    dsigma   dh3      dh4")
        print("".join("%8.2g" % f for f in errors))
        print('Elapsed time in pPXF: %.2f s' % (clock() - t))
        prec = int(1 - np.floor(np.log10(redshift_err)))  # two digits of uncertainty
        print(f"Best-fitting redshift z = {redshift_fit:#.{prec}f} "
            f"+/- {redshift_err:#.{prec}f}")

        #output kinematics parameters
        kinematics = pp.sol
        error_kinematics = errors

        #output fit_model
        bestfit_flux = pp.bestfit
        bestfit_wavelength = ln_lam1

        return kinematics, error_kinematics, bestfit_flux, bestfit_wavelength
    except:
        print ('The selected template does not cover the wavelength range you want to fit')
    #kinematics=info_pop=info_pop_mass= mass_light= errors= bestfit_flux= bestfit_wave= bestfit_gas_flux= chi_square= error_age_abs= error_met= error_mass_age_abs= error_mass_met= emission_corrected_flux = 0
            

    
    
#stellar populations with ppxf
def ppxf_pop(wave, flux, wave1, wave2, FWHM_gal, z, sigma_guess, fit_components, with_plots, with_errors, save_plot, spec_name, regul_err, additive_degree, multiplicative_degree, tied_balmer, stellar_library, gas_reddening):
    
    ppxf_dir = path.dirname(path.realpath(lib.__file__))
    #print(ppxf_dir)
    #fit_components: 'with_gas' or 'without_gas' or 'balmer' or 'balmer_forbidden'
    #converting wavelength to angstrom
    wave = wave*10
    galaxy = flux
    galaxy = galaxy/np.median(galaxy)

    #selecting the input range
    wave1 = wave1*10
    wave2 = wave2*10
    line_wave = wave[(wave >= wave1) & (wave <= wave2)]  
    line_flux_spec = galaxy[(wave >= wave1) & (wave <= wave2)]
    
    #updating the variables
    galaxy = line_flux_spec
    wave = line_wave

    print('Rebinning to log')
    galaxy, ln_lam1, velscale = util.log_rebin(wave, galaxy)
    lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1 + z)
    wave = np.exp(ln_lam1)

    noise = np.full_like(galaxy, 0.0163) #noise per pixel
    c = 299792.458
    velscale = c*np.log(wave[1]/wave[0])

    #loading the templates EMILES
    #sps_name = 'emiles'
    sps_name = stellar_library
    #if sps_name == 'emiles':
        #lam_range_temp = [3500, 1e4]
    basename = f"spectra_{sps_name}_9.0.npz"
    filename = path.join(ppxf_dir, 'sps_models', basename)
    if not path.isfile(filename):
        url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
        request.urlretrieve(url, filename)

    #if sps_name == 'emiles':

    sps = lib.sps_lib(filename, velscale, norm_range=[5070, 5950])
    #else:
        #sps = lib.sps_lib(filename, velscale, norm_range=[5070, 5950])

    reg_dim = sps.templates.shape[1:]
    stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)

    #gas_reddening = reddening
    gas_reddening = 0.1
    #errors
    error_age = 0
    error_met = 0
    error_age_abs = 0
    error_mass_age = 0
    error_mass_met = 0
    error_mass_age_abs = 0

  ###################### Now without gas ##################  
    if fit_components == 'without_gas':
        print ('Fitting without gas component')

        try:
            templates = stars_templates
            vel = c*np.log(1 + z)
            start = [vel, sigma_guess]
            n_temps = stars_templates.shape[1]
            component = 0
            gas_component = np.array(component) > 0
            moments = 4
            start = start

            t = clock()
            #do the fit!
            pp = ppxf(templates, galaxy, noise, velscale, start,
                moments=moments, degree= additive_degree, mdegree=multiplicative_degree,
                lam=wave, lam_temp=sps.lam_temp,
                regul=1/regul_err, reg_dim=reg_dim,
                component=component)

            #setting up the result parameters
            light_weights = pp.weights[~gas_component]
            light_weights = light_weights.reshape(reg_dim)
            light_weights /= light_weights.sum()
            print('')
            print('Luminosity weighted stellar populations:')
            info_pop = sps.mean_age_metal(light_weights)
            mass_weights = light_weights/sps.flux
            mass_weights /= mass_weights.sum()              # Normalize to mass fractions
            mass_light = sps.mass_to_light(mass_weights, band="v")
            print('')
            print('Mass weighted stellar populations:')
            info_pop_mass = sps.mean_age_metal(mass_weights)

            #printing the output infos
            print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
            print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
            print(f"Chi^2: {(pp.chi2):#.4g}")
            print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

            #extracting the output parameters
            kinematics = pp.sol
            bestfit_flux = pp.bestfit
            bestfit_wave = wave
            bestfit_gas_flux = 0.
            chi_square = pp.chi2
            emission_corrected_flux = galaxy
            errors = pp.error*np.sqrt(pp.chi2)

            if with_plots == True:
                #setting up the plot
                plt.clf()
                pp.plot()
                plt.title('Template fitting and residuals')
                plt.show()
                # Plot stellar population luminosity-fraction distribution
                plt.subplot(211)
                sps.plot(light_weights,cmap='BuPu', title = 'Luminosity fraction')
                plt.tight_layout()

                plt.subplot(212)
                sps.plot(mass_weights,cmap='BuPu', title = 'Mass fraction' )
                plt.tight_layout()


                plt.show()

            if save_plot == True:
                            #setting up the plot
                result_plot_dir = 'results/plots'
                #pathlib.Path(result_plot_dir).mkdir(parents=True, exist_ok=True)
                os.makedirs(result_plot_dir, exist_ok=True)

                plt.clf()
                pp.plot()
                plt.title('Template fitting and residuals')
                plt.savefig(result_plot_dir + '/'+ 'pop_fit_ppxf_'+ spec_name + '.eps', format='eps', dpi=100)

                # Plot stellar population luminosity-fraction distribution
                plt.subplot(211)
                sps.plot(light_weights,cmap='BuPu', title = 'Luminosity fraction')
                plt.tight_layout()

                plt.subplot(212)
                sps.plot(mass_weights,cmap='BuPu', title = 'Mass fraction' )
                plt.tight_layout()
                plt.savefig(result_plot_dir + '/'+ 'pop_pop_ppxf_'+ spec_name + '.eps', format='eps', dpi=100)
                plt.close()



            if with_errors == True: #calculating the errors of age and metallicity with MonteCarlo simulations
                residual = galaxy - bestfit_flux
                snr = 1/np.std(residual)
                print ('snr spectrum:', snr)

                n_sim = 10 #how many simulated templates I want to create. Watch out for the computation time!
                age_dist = []
                met_dist = []
                mass_age_dist = []
                mass_met_dist = []
                #as done by Cappellari: https://github.com/micappe/ppxf_examples/blob/main/ppxf_example_population_bootstrap.ipynb
                for i in range(n_sim):
                    noisy_template = add_noise(bestfit_wave, bestfit_flux, snr)

                    #no regularization!
                    pp = ppxf(templates, noisy_template, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=bestfit_wave, lam_temp=sps.lam_temp,
                    component=component)

                        #setting up the result parameters
                    light_weights_err = pp.weights[~gas_component]
                    light_weights_err = light_weights_err.reshape(reg_dim)
                    light_weights_err /= light_weights_err.sum()
                    print('')
                    print('Luminosity weighted stellar populations:')
                    info_pop_err = sps.mean_age_metal(light_weights_err)
                    mass_weights_err = light_weights_err/sps.flux
                    mass_weights_err /= mass_weights_err.sum()              # Normalize to mass fractions
                    mass_light_err = sps.mass_to_light(mass_weights_err, band="v")
                    print('')
                    print('Mass weighted stellar populations:')
                    info_pop_mass_err = sps.mean_age_metal(mass_weights_err)


                    #printing the output infos
                    print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                    print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                    print(f"Chi^2: {(pp.chi2):#.4g}")
                    print(f"Elapsed time in pPXF: {(clock() - t):.2f}")
                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])

                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])
                    mass_age_dist.append(info_pop_mass_err[0])
                    mass_met_dist.append(info_pop_mass_err[1])


                error_age = np.std(age_dist)
                error_age_abs = np.log(10)*error_age*(10**info_pop_err[0])/1e9
                error_met = np.std(met_dist)
                error_mass_age = np.std(mass_age_dist)
                error_mass_age_abs = np.log(10)*error_mass_age*(10**info_pop_mass_err[0])/1e9
                error_mass_met = np.std(mass_met_dist)

                print(f"Error luminosity age (dex): {(error_age):#.4g}")
                print(f"Error luminosity age (Gyr): {(error_age_abs):#.4g}")
                print(f"Error luminosity met (dex): {(error_met):#.4g}")
                print(f"Error mass age (dex): {(error_mass_age):#.4g}")
                print(f"Error mass age (Gyr): {(error_mass_age_abs):#.4g}")
                print(f"Error mass met (dex): {(error_mass_met):#.4g}")


            return kinematics, info_pop, info_pop_mass, mass_light, errors, bestfit_flux, bestfit_wave, bestfit_gas_flux, chi_square, error_age_abs, error_met, error_mass_age_abs, error_mass_met, emission_corrected_flux

        except AssertionError:
            print ('The selected template does not cover the wavelength range you want to fit')
            kinematics=info_pop=info_pop_mass= mass_light= errors= bestfit_flux= bestfit_wave= bestfit_gas_flux= chi_square= error_age_abs= error_met= error_mass_age_abs= error_mass_met= emission_corrected_flux = 0

#################### WITH GAS #########################


    if fit_components == 'with_gas':
        
        print ('Fitting with at least one gas component')
        try:
            tie_balmer=tied_balmer
            limit_doublets=False
            check_gas_cond = 1
            #retrieving the emission lines in the wavelength range
            gas_templates, gas_names, line_wave = emission_lines(
            sps.ln_lam_temp, lam_range_gal, FWHM_gal,
            tie_balmer=tie_balmer, limit_doublets=limit_doublets)

            # grouping the emission lines: 1) balmer, 2) forbidden, 3) others
            n_forbidden = np.sum(["[" in a for a in gas_names])
            if tie_balmer == False:
                n_balmer = np.sum(["(" in a for a in gas_names])
            else:
                n_balmer = np.sum(["Balmer" in a for a in gas_names])
                print ('Tied Balmer lines')

            n_others = np.sum(["-" in a for a in gas_names])


            #looking for the existence of at least one line of each group in the selected spectral window
            if n_forbidden !=0 and n_balmer !=0 and n_others !=0:
                ##### THREE GAS COMPONETS
                print('Balmer, forbidden and other lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden +[3]*n_others
                gas_component = np.array(component) > 0
                moments = [4, 2, 2, 2]
                start = [start, start, start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden !=0 and n_balmer !=0 and n_others == 0:
                #####
                print ('Forbidden and Balmer lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden !=0 and n_balmer == 0 and n_others !=0:
                #####
                print ('Forbidden and other lines')
                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_others + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)
                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden !=0 and n_balmer == 0 and n_others ==0:
                #######
                print ('Only forbidden lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden ==0 and n_balmer != 0 and n_others ==0:
                ######
                print('Only balmer lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_balmer
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden ==0 and n_balmer != 0 and n_others !=0:
                #######
                print ('Balmer and other lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_balmer [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden ==0 and n_balmer == 0 and n_others !=0:
                ########
                print ('Only other lines')

                templates = np.column_stack([stars_templates, gas_templates])
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]

                component = [0]*n_temps + [1]*n_others
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]
                #gas_reddening = None
                t = clock()

                #fitting
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, gas_reddening=gas_reddening)

                        #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = pp.gas_bestfit
                chi_square = pp.chi2

            if n_forbidden ==0 and n_balmer == 0 and n_others ==0:
                ########### NO GAS COMPONENT
                print ('No gas lines found')
                check_gas_cond = 0
                templates = stars_templates
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]
                component = 0
                gas_component = np.array(component) > 0
                moments = 4
                start = start

                t = clock()
                #do the fit!
                pp = ppxf(templates, galaxy, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=1/regul_err, reg_dim=reg_dim,
                    component=component)

                #setting up the result parameters
                light_weights = pp.weights[~gas_component]
                light_weights = light_weights.reshape(reg_dim)
                light_weights /= light_weights.sum()
                print('')
                print('Luminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights)
                mass_weights = light_weights/sps.flux
                mass_weights /= mass_weights.sum()              # Normalize to mass fractions
                mass_light = sps.mass_to_light(mass_weights, band="v")
                print('')
                print('Mass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights)

                #printing the output infos
                print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                print(f"Chi^2: {(pp.chi2):#.4g}")
                print(f"Elapsed time in pPXF: {(clock() - t):.2f}")

                #extracting the output parameters
                kinematics = pp.sol
                bestfit_flux = pp.bestfit
                bestfit_wave = wave
                bestfit_gas_flux = 0.
                chi_square = pp.chi2

            errors = pp.error[0]*np.sqrt(pp.chi2)
            try:
                emission_corrected_flux = galaxy - pp.gas_bestfit #pp.bestfit - pp.gas_bestfit
            except TypeError:
                emission_corrected_flux = galaxy


            if with_plots == True:

                #setting up the plot
                plt.clf()
                pp.plot()
                plt.title('Template fitting and residuals')
                plt.show()

                # Plot stellar population luminosity-fraction distribution
                plt.subplot(211)
                sps.plot(light_weights,cmap='BuPu', title = 'Luminosity fraction')
                plt.tight_layout()

                plt.subplot(212)
                sps.plot(mass_weights,cmap='BuPu', title = 'Mass fraction' )
                plt.tight_layout()
                plt.show()

            if save_plot == True:
                            #setting up the plot
                result_plot_dir = 'results/plots'
                #pathlib.Path(result_plot_dir).mkdir(parents=True, exist_ok=True)
                os.makedirs(result_plot_dir, exist_ok=True)

                plt.clf()
                pp.plot()
                plt.title('Template fitting and residuals')
                plt.savefig(result_plot_dir + '/'+ 'pop_fit_ppxf_'+ spec_name + '.eps', format='eps', dpi=100)

                # Plot stellar population luminosity-fraction distribution
                plt.subplot(211)
                sps.plot(light_weights,cmap='BuPu', title = 'Luminosity fraction')
                plt.tight_layout()

                plt.subplot(212)
                sps.plot(mass_weights,cmap='BuPu', title = 'Mass fraction' )
                plt.tight_layout()
                plt.savefig(result_plot_dir + '/'+ 'pop_pop_ppxf_'+ spec_name + '.eps', format='eps', dpi=100)
                plt.close()

            if with_errors == True: #uncertainties for age and metallicity with MonteCarlo simulations
                residual = galaxy - bestfit_flux
                snr = 1/np.std(residual)
                print ('snr spectrum:', snr)

                n_sim = 10 #number of simulated noisy templates to create
                age_dist = []
                met_dist = []
                mass_age_dist = []
                mass_met_dist = []
                for i in range(n_sim):
                    noisy_template = add_noise(bestfit_wave, bestfit_flux, snr)

                    pp = ppxf(templates, noisy_template, noise, velscale, start,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=bestfit_wave, lam_temp=sps.lam_temp,
                    component=component)

                        #setting up the result parameters
                    light_weights_err = pp.weights[~gas_component]
                    light_weights_err = light_weights_err.reshape(reg_dim)
                    light_weights_err /= light_weights_err.sum()
                    print('')
                    print('Luminosity weighted stellar populations:')
                    info_pop_err = sps.mean_age_metal(light_weights_err)
                    mass_weights_err = light_weights_err/sps.flux
                    mass_weights_err /= mass_weights_err.sum()              # Normalize to mass fractions
                    mass_light_err = sps.mass_to_light(mass_weights_err, band="v")
                    print('')
                    print('Mass weighted stellar populations:')
                    info_pop_mass_err = sps.mean_age_metal(mass_weights_err)

                    #printing the output infos
                    print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                    print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                    print(f"Chi^2: {(pp.chi2):#.4g}")
                    print(f"Elapsed time in pPXF: {(clock() - t):.2f}")
                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])

                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])
                    mass_age_dist.append(info_pop_mass_err[0])
                    mass_met_dist.append(info_pop_mass_err[1])

                error_age = np.std(age_dist)
                error_age_abs = np.log(10)*error_age*(10**info_pop_err[0])/1e9
                error_met = np.std(met_dist)
                error_mass_age = np.std(mass_age_dist)
                error_mass_age_abs = np.log(10)*error_mass_age*(10**info_pop_mass_err[0])/1e9
                error_mass_met = np.std(mass_met_dist)

                print(f"Error luminosity age (dex): {(error_age):#.4g}")
                print(f"Error luminosity age (Gyr): {(error_age_abs):#.4g}")
                print(f"Error luminosity met (dex): {(error_met):#.4g}")
                print(f"Error mass age (dex): {(error_mass_age):#.4g}")
                print(f"Error mass age (Gyr): {(error_mass_age_abs):#.4g}")
                print(f"Error mass met (dex): {(error_mass_met):#.4g}")

            plt.close()

            return kinematics, info_pop, info_pop_mass, mass_light, errors, bestfit_flux, bestfit_wave, bestfit_gas_flux, chi_square, error_age_abs, error_met, error_mass_age_abs, error_mass_met, emission_corrected_flux
        except AssertionError:
            print ('The selected template does not cover the wavelength range you want to fit')
            kinematics=info_pop=info_pop_mass= mass_light= errors= bestfit_flux= bestfit_wave= bestfit_gas_flux= chi_square= error_age_abs= error_met= error_mass_age_abs= error_mass_met= emission_corrected_flux = 0
    
    
#Crosscorr, adapted from  pyasl, but working
def crosscorrRV(w, f, tw, tf, rvmin, rvmax, drv, mode="doppler", skipedge=0, edgeTapering=None):


  # Copy and cut wavelength and flux arrays
  w, f = w.copy(), f.copy()
  if skipedge > 0:
    w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]
  
  if edgeTapering is not None:
    # Smooth the edges using a sine
    if isinstance(edgeTapering, float):
      edgeTapering = [edgeTapering, edgeTapering]

    indi = np.where(w < w[0]+edgeTapering[0])[0]
    f[indi] *= np.sin((w[indi] - w[0])/edgeTapering[0]*np.pi/2.0)
    # Carry out edge tapering (right edge)
    indi = np.where(w > (w[-1]-edgeTapering[1]))[0]
    f[indi] *= np.sin((w[indi] - w[indi[0]])/edgeTapering[1]*np.pi/2.0 + np.pi/2.0)
  
  # Speed of light in km/s
  c = 299792.458

  # Calculate the cross correlation
  drvs = np.arange(rvmin, rvmax, drv)
  cc = np.zeros(len(drvs))
  for i, rv in enumerate(drvs):
    if mode == "lin":
      # Shift the template linearly
      fi = interpolate.interp1d(tw+meanWl*(rv/c), tf)
    elif mode == "doppler":
      # Apply the Doppler shift
      fi = interpolate.interp1d(tw*(1.0 + rv/c), tf)
    # Shifted template evaluated at location of spectrum
    cc[i] = np.sum(f * fi(w))

  return drvs, cc





#Modified emission line list from the ppxf_util.py of the ppxf package.

def emission_lines(ln_lam_temp, lam_range_gal, FWHM_gal, pixel=True,
                   tie_balmer=False, limit_doublets=False, vacuum=False):
    """
    Generates an array of Gaussian emission lines to be used as gas templates in PPXF.

    ****************************************************************************
    ADDITIONAL LINES CAN BE ADDED BY EDITING THE CODE OF THIS PROCEDURE, WHICH
    IS MEANT AS A TEMPLATE TO BE COPIED AND MODIFIED BY THE USERS AS NEEDED.
    ****************************************************************************


    Output Parameters
    -----------------

    emission_lines: ndarray
        Array of dimensions ``[ln_lam_temp.size, line_wave.size]`` containing
        the gas templates, one per array column.

    line_names: ndarray
        Array of strings with the name of each line, or group of lines'

    line_wave: ndarray
        Central wavelength of the lines, one for each gas template'

    """
    #        Balmer:     H10       H9         H8        Heps    Hdelta    Hgamma    Hbeta     Halpha
    balmer = np.array([3798.983, 3836.479, 3890.158, 3971.202, 4102.899, 4341.691, 4862.691, 6564.632])  # vacuum wavelengths

    if tie_balmer:

        # Balmer decrement for Case B recombination (T=1e4 K, ne=100 cm^-3)
        # from Storey & Hummer (1995) https://ui.adsabs.harvard.edu/abs/1995MNRAS.272...41S
        # In electronic form https://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/64
        # See Table B.7 of Dopita & Sutherland 2003 https://www.amazon.com/dp/3540433627
        # Also see Table 4.2 of Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/
        wave = balmer
        if not vacuum:
            wave = util.vac_to_air(wave)
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
        ratios = np.array([0.0530, 0.0731, 0.105, 0.159, 0.259, 0.468, 1, 2.86])
        ratios *= wave[-2]/wave  # Account for varying log-sampled pixel size in Angstrom
        emission_lines = gauss @ ratios
        line_names = ['Balmer']
        w = (lam_range_gal[0] < wave) & (wave < lam_range_gal[1])
        line_wave = np.mean(wave[w]) if np.any(w) else np.mean(wave)

    else:

        line_wave = balmer
        if not vacuum:
            line_wave = util.vac_to_air(line_wave)
        line_names = ['(H10)', '(H9)', '(H8)', '(Heps)', '(Hdelta)', '(Hgamma)', '(Hbeta)', '(Halpha)']
        emission_lines = util.gaussian(ln_lam_temp, line_wave, FWHM_gal, pixel)

    if limit_doublets:

        # The line ratio of this doublet lam3727/lam3729 is constrained by
        # atomic physics to lie in the range 0.28--1.47 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[OII]-----
        wave = [3727.092, 3729.875]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[OII]3726_d1', '[OII]3726_d2']
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
        doublets = gauss @ [[1, 1], [0.28, 1.47]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

        # The line ratio of this doublet lam6717/lam6731 is constrained by
        # atomic physics to lie in the range 0.44--1.43 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #        -----[SII]-----
        wave = [6718.294, 6732.674]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[SII]6731_d1', '[SII]6731_d2']
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
        doublets = gauss @ [[0.44, 1.43], [1, 1]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    else:

        # Here the two doublets are free to have any ratio
        #         -----[OII]-----     -----[SII]-----
        wave = [3727.092, 3729.875, 6718.294, 6732.674]  # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
        emission_lines = np.column_stack([emission_lines, gauss])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    # Here the lines are free to have any ratio
    #       -----[NeIII]-----    HeII      HeI
    wave = [3968.59, 3869.86, 4687.015, 5877.243]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['[NeIII]3968', '[NeIII]3869', '-HeII4687-', '-HeI5876-']
    gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)

    # NIR H lines
    #       paeps      pad      pab
    wave = [10052.1, 10941.1, 12821.6]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['-PaEps-', '-Pad-', '-Pab-']
    gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)
    

    
    ######### Doublets with fixed ratios #########

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OIII]-----
    wave = [4960.295, 5008.240]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OI]-----
    wave = [6302.040, 6365.535]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [1, 0.33]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[NII]-----
    wave = [6549.860, 6585.271]    # air wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # Only include lines falling within the estimated fitted wavelength range.
    #
    w = (lam_range_gal[0] < line_wave) & (line_wave < lam_range_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave




# BLACK BODY FITTING FUNCTIONS. Adapted from: https://github.com/Professor-G/BlackbodyFit

def select_data(spect_class):
    """This function selects the flux data used to calibrate and
    convert the instrumental flux to flux. If there's no data for the
    specific spectral class, the function picks the closest spectral class
    to the input one.

    :param spect_class: Spectral class. Must be a string. Ex) 'A5'

    :rtype: array
    """
    flux_data = np.loadtxt('all_spectra.txt')

    if spect_class == 'O' or spect_class == 'O5' or spect_class == 'O4' or spect_class == 'O3' or spect_class == 'O2' or spect_class == 'O1':
        data = flux_data[:,1]
        print('Intensity data for star O5 V HD 242908 selected')

    elif spect_class == 'O6':
        data = flux_data[:,3]
        print('Intensity data for star O6.5 V HD 12993 selected')

    elif spect_class == 'O7':
        data = flux_data[:,4]
        print('Intensity data for star O7 V HD 35619 selected')

    elif spect_class == 'O8':
        data = flux_data[:,6]
        print('Intensity data for star O8 V HD 242935 selected')

    elif spect_class == 'O9':
        data = flux_data[:,8]
        print('Intensity data for star O9 V HD 17520 selected')

    elif spect_class == 'B' or spect_class == 'B4':
        data = flux_data[:,15]
        print('Intensity data for star B4 V FEIGE 40 selected')

    elif spect_class == 'B0':
        data = flux_data[:,11]
        print('Intensity data for star B0 V HD 158659 selected')

    elif spect_class == 'B1' or spect_class == 'B2':
        data = flux_data[:,13]
        print('Intensity data for star  B1.5  V HD  35215 selected')

    elif spect_class == 'B3' or spect_class == 'B4':
        data = flux_data[:,14]
        print('Intensity data for star   B3    V HD  37767 selected')

    elif spect_class == 'B5' or spect_class == 'B6':
        data = flux_data[:,17]
        print('Intensity data for star   B6    V HD  30584 selected')

    elif spect_class == 'B7' or spect_class == 'B8' or spect_class == 'B9':
        data = flux_data[:,18]
        print('Intensity data for star    B8    VO    1015 selected')

    elif spect_class == 'A0' or spect_class == 'A1':
        data = flux_data[:,19]
        print('Intensity data for star   A1    V HD 116608 selected')

    elif spect_class == 'A2':
        data = flux_data[:,21]
        print('Intensity data for star    A2    V HD 124320 selected')

    elif spect_class == 'A3':
        data = flux_data[:,24]
        print('Intensity data for star    A3    V HD 221741 selected')

    elif spect_class == 'A4' or spect_class == 'A5' or spect_class == 'A':
        data = flux_data[:,25]
        print('Intensity data for star   A5    V HD   9547 selected')

    elif spect_class == 'A6':
        data = flux_data[:,26]
        print('Intensity data for star   A6    V HD  21619 selected')

    elif spect_class == 'A7':
        data = flux_data[:,27]
        print('Intensity data for star  A7    V HD  23863 selected')

    elif spect_class == 'A8':
        data = flux_data[:,29]
        print('Intensity data for star  A8    V HD   9972 selected')

    elif spect_class == 'A9':
        data = flux_data[:,30]
        print('Intensity data for star  A9    V HD  23733 selected')

    elif spect_class == 'F0' or spect_class == 'F1':
        data = flux_data[:,31]
        print('Intensity data for star   F0    V HD  10032 selected')

    elif spect_class == 'F2' or spect_class == 'F3':
        data = flux_data[:,32]
        print('Intensity data for star   F3    V HZ    948 selected')

    elif spect_class == 'F4':
        data = flux_data[:,33]
        print('Intensity data for star  F4    V HD  23511 selected')

    elif spect_class == 'F5' or spect_class == 'F':
        data = flux_data[:,34]
        print('Intensity data for star   F5    V HZ    227 selected')

    elif spect_class == 'F6':
        data = flux_data[:,35]
        print('Intensity data for star   F6    V SAO 57199 selected')

    elif spect_class == 'F7':
        data = flux_data[:,37]
        print('Intensity data for star  F7    V HD   5702 selected')

    elif spect_class == 'F8':
        data = flux_data[:,40]
        print('Intensity data for star   F8    V HD   6111 selected')

    elif spect_class == 'F9':
        data = flux_data[:,41]
        print('Intensity data for star    F9    V HD  31084 selected')

    elif spect_class == 'F9':
        data = flux_data[:,41]
        print('Intensity data for star    F9    V HD  31084 selected')

    elif spect_class == 'G0':
        data = flux_data[:,43]
        print('Intensity data for star    G0    V HD  28099 selected')

    elif spect_class == 'G1':
        data = flux_data[:,44]
        print('Intensity data for star     G1    V HD  17647 selected')

    elif spect_class == 'G2':
        data = flux_data[:,45]
        print('Intensity data for star   G2    V HD  66171  selected')

    elif spect_class == 'G3':
        data = flux_data[:,46]
        print('Intensity data for star    G3    V BD+581199 selected')

    elif spect_class == 'G4':
        data = flux_data[:,47]
        print('Intensity data for star   G4    V TR A   14 selected')

    elif spect_class == 'G5' or spect_class == 'G6' or spect_class == 'G':
        data = flux_data[:,48]
        print('Intensity data for star  G6    V HD  22193 selected')

    elif spect_class == 'G7':
        data = flux_data[:,49]
        print('Intensity data for star  G7    V HD  27685 selected')

    elif spect_class == 'G8' or spect_class == 'G9':
        data = flux_data[:,50]
        print('Intensity data for star   G9    V HD  33278 selected')

    elif spect_class == 'K0' or spect_class == 'K1' or spect_class == 'K2':
        data = flux_data[:,52]
        print('Intensity data for star  K0    V HD  23524 selected')

    elif spect_class == 'K3' or spect_class == 'K4':
        data = flux_data[:,53]
        print('Intensity data for star  K4    V HD   5351 selected')

    elif spect_class == 'K' or spect_class == 'K5' or spect_class == 'K6' or spect_class == 'K7' or spect_class == 'K8' or spect_class == 'K9':
        data = flux_data[:,54]
        print('Intensity data for star  K5    V SAO 76803 selected')

    elif spect_class == 'M0':
        data = flux_data[:,55]
        print('Intensity data for star   M0    V HD 260655 selected')

    elif spect_class == 'M1' or spect_class == 'M2' or spect_class == 'M3':
        data = flux_data[:,56]
        print('Intensity data for star  M1    V BD+63 137 selected')

    elif spect_class == 'M' or spect_class == 'M4' or spect_class == 'M5' or spect_class == 'M6' or spect_class == 'M7' or spect_class == 'M8' or spect_class == 'M9':
        data = flux_data[:,57]
        print('Intensity data for star   M5    V YALE 1755 selected')

    else:
        data = flux_data[:,31]
        print('Spectral class not identified. Intensity data for star   F0    V HD  10032 selected')

    return data

def blackbody(wavelength, T):
    """
    Planck's law, which describes the black body radiation
    of a source in thermal equilibrium at a given temperature T.
    """
    return 2*h*c**2 / (wavelength**5 * (np.e**(h*c / (wavelength*k*T)) - 1))

#Blackbody fitting
def blackbody_fit(wavelength, flux, initial_wave, final_wave, t_guess, with_plots):

    #extracting the wavelength and the flux within the selected band
    wave_bb = wavelength[(wavelength >= initial_wave) & (wavelength <= final_wave)]
    flux_bb = flux[(wavelength >= initial_wave) & (wavelength <= final_wave)]

    wavelength = wave_bb
    flux = flux_bb

    #normalising the flux to the median wavelength
    median_wave = np.median(wavelength)
    epsilon_wave = 5.

    norm_flux = norm_spec(wavelength, flux, median_wave, epsilon_wave, flux)
    flux = norm_flux

    wavelength = wavelength*1e-9 #convert to meters for SI consistency

    def blackbody_flux(wave, flux, T):
        blackbody_flux = blackbody(wave, T)

        flux_integral = np.abs(np.trapz(flux, wave))
        planck_integral = np.abs(quad(blackbody, np.min(wave), np.max(wave), args = T ))[0]
        scale_factor = flux_integral / planck_integral

        return scale_factor*blackbody_flux

    def residuals(T, y, lam):
        return y - blackbody_flux(wavelength, flux, T)

    t0 = np.array([t_guess]) #the initial temperature guess for the optimization
    T = leastsq(residuals, t0, args = (flux, wavelength))[0].astype('float32')
    bbody_fit = blackbody(wavelength, T)

    instrumental_integral = np.abs(np.trapz(flux, wavelength))
    planck_integral = np.abs(quad(blackbody, np.min(wavelength), np.max(wavelength), args = T))[0]
    scale_factor = instrumental_integral / planck_integral
    y = scale_factor*bbody_fit
    wavelength = wavelength*1e9
    temperature = int(round(T[0]))
    residual_bb = flux-y

    if with_plots == True:
        fig = plt.figure()
        ax = fig.add_subplot(211)
        #plt.subplot(211)
        plt.plot(wavelength, flux, label = "DATA")
        plt.plot(wavelength, y, 'r-',  label = "FIT")
        plt.title("Blackbody fitting, T = " + str(temperature) +' K' )
        plt.xlabel("$\lambda$ (nm)")
        plt.ylabel("Relative flux")
        plt.legend(loc = 1)
        plt.tight_layout()

        ax2 = fig.add_subplot(212)
        plt.title("Residuals" )
        plt.ylabel("Relative flux")
        plt.xlabel("$\lambda$ (nm)")
        plt.plot(wavelength, residual_bb, 'g.', label = "residuals")
        plt.tight_layout()

        plt.show()

    return temperature, residual_bb



#FUCTIONS FOR THE FITS HEADER MANIPULATION

#for the modification of just one file
def read_fits_header(file_path):
    try:
        with fits.open(file_path) as hdul:
            header = hdul[0].header
        return header
    except Exception as e:
        return str(e)

def save_fits_header(file_path, header):
    try:
        with fits.open(file_path, mode='update') as hdul:
            hdul[0].header = header
            hdul.flush()
        return True
    except Exception as e:
        return str(e)

def delete_keyword(header, key):
    try:
        del header[key]
        return True
    except KeyError:
        return f"Keyword '{key}' not found in the header."


#for the manipulation of a list of fits with a list of keywords:
def read_keyword_values_from_file(file_path):
    key_name = []
    key_value = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                key, value, value_type = map(str.strip, line.split('='))
                # Recognise the value type and convert
                if value_type.lower() == 'int':
                    value = int(value)
                elif value_type.lower() == 'float':
                    value = float(value)
                elif value_type.lower() == 'string':
                    value = str(value)

                # Add the value to the list
                key_name.append(key)
                key_value.append(value)

    return key_name, key_value


def read_file_list(file_path):
    try:
        with open(file_path, 'r') as file:
            file_paths = [line.strip() for line in file if not line.startswith('#')]
        return file_paths
    except Exception as e:
        return str(e)


#For the extraction and saving of the new keywords in a list of fits files:
def extract_keyword(file_path, keyword):
    try:
        with fits.open(file_path) as hdul:
            header = hdul[0].header
            if keyword in header:
                return header[keyword]
            else:
                return f"Keyword '{keyword}' non trovata nel file."
    except Exception as e:
        return str(e)


def save_to_text_file(data, output_file):
    with open(output_file, 'w') as file:
        # writing the header
        file.write(f"#Spectrum {data[0]['keyword']}\n")

        for entry in data:
            file.write(f"{entry['file']} {entry['value']}\n")

#**************** FUNCTIONS FOR PLOTTING *********************
def linear_fit(x, m, b):
    return m * x + b

def get_column_names(file_path):
    try:
        data = pd.read_csv(file_path, sep=None, engine='python')
        return list(data.columns)
    except Exception as e:
        sg.popup_error(f'Error reading file: {str(e)}')
        return []

def plot_data(file_path, x_column, y_columns, x_label, y_label, marker_color, marker_size, plot_scale, x_label_size,
              y_label_size, x_tick_size, y_tick_size, legend, add_error_bars_x, add_error_bars_y, x_err, y_err, saveps,
              enable_linear_fit, x_log_scale, y_log_scale, x_range_min=None, x_range_max=None, y_range_min=None,
              y_range_max=None):

    try:
        #Load the data
        data = pd.read_csv(file_path, sep=None, engine='python')

        #plotting
        plt.figure(figsize=plot_scale)
        if x_log_scale:
            plt.xscale('log')
        if y_log_scale:
            plt.yscale('log')

        for column in y_columns:
            # Extract x and y values
            x = data[x_column].values
            y = data[column].values

            # Handlling the error bars
            if add_error_bars_y == True and add_error_bars_x == False :
                error_bar_data_y = data[y_err].values
                #I need 1D arrays
                x = np.squeeze(x)
                y = np.squeeze(y)
                error_bar_data_y = np.squeeze(error_bar_data_y)
                #Adding the error bars
                plt.scatter(x, y, label=column, color=marker_color, s=marker_size)
                plt.errorbar(x, y, yerr=error_bar_data_y, linestyle='None', ecolor = 'black', capsize=2)
            elif add_error_bars_x == True and add_error_bars_y == False:
                error_bar_data_x = data[x_err].values
                #I need 1D arrays
                x = np.squeeze(x)
                y = np.squeeze(y)
                error_bar_data_x = np.squeeze(error_bar_data_x)
                #Adding the error bars
                plt.scatter(x, y, label=column, color=marker_color, s=marker_size)
                plt.errorbar(x, y, xerr=error_bar_data_x, linestyle='None', ecolor = 'black', capsize=2)
            elif (add_error_bars_y == True and add_error_bars_x == True):
                error_bar_data_y = data[y_err].values
                error_bar_data_x = data[x_err].values
                x = np.squeeze(x)
                y = np.squeeze(y)
                error_bar_data_y = np.squeeze(error_bar_data_y)
                error_bar_data_x = np.squeeze(error_bar_data_x)
                #Adding the error bars
                plt.scatter(x, y, label=column, color=marker_color, s=marker_size)
                plt.errorbar(x, y, xerr=error_bar_data_x, yerr=error_bar_data_y, linestyle='None', ecolor = 'black', capsize=2)
            else:
                plt.scatter(x, y, label=column, color=marker_color, s=marker_size)

        if enable_linear_fit:
            popt, _ = curve_fit(linear_fit,  x.reshape(-1), y.reshape(-1))
            fit_x = np.linspace(min(x), max(x), 100)
            fit_y = linear_fit(fit_x, *popt)
            linear_regression = scipy.stats.linregress(x.reshape(-1), y.reshape(-1))
            pearson = linear_regression.rvalue
            pearson_str = ('R = ' + str(round(pearson,2)))
            x_data = x.reshape(-1)
            y_data = y.reshape(-1)

            #bootstrap for the error on R
            num_bootstrap_samples = 1000
            #removing the Nans
            nan_indices_xy = np.isnan(x_data) | np.isnan(y_data)
            x_data = x_data[~nan_indices_xy]
            y_data = y_data[~nan_indices_xy]
            bootstrap_results_xy = np.zeros(num_bootstrap_samples)
            for i in range(num_bootstrap_samples):
                #a) for the correlations with Mg2
                indices_xy = np.random.choice(len(x_data), len(y_data), replace=True)
                x_bootstrap = x_data[indices_xy]
                y_bootstrap = y_data[indices_xy]
                # Calcolo del coefficiente di correlazione per il campione bootstrap
                bootstrap_results_xy[i], _ = pearsonr(x_bootstrap, y_bootstrap)

            #mean_bootstrap_xy = np.mean(bootstrap_results_xy)
            std_bootstrap_xy = str(round(np.std(bootstrap_results_xy),2))
            plt.plot(fit_x, fit_y, linestyle='-', color='red', label=(pearson_str + '$\pm$'+std_bootstrap_xy))

        #Make the plot nice
        plt.xlabel(x_label, fontsize=x_label_size)
        plt.ylabel(y_label, fontsize=y_label_size)
        plt.xticks(fontsize=x_tick_size)
        plt.yticks(fontsize=y_tick_size)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        if x_range_min is not None and x_range_max is not None:
            plt.xlim(float(x_range_min), float(x_range_max))
        if y_range_min is not None and y_range_max is not None:
            plt.ylim(float(y_range_min), float(y_range_max))


        if legend:
            plt.legend()
        if saveps:
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            plt.savefig('plot_'+timestamp+'.eps', format='eps', dpi=300)
            sg.Popup ('eps file saved with success!')
        else:
            plt.show()

    except Exception as e:
        sg.popup_error(f'Error plotting data: {str(e)}')

    finally:
        plt.close()

#****************** FUNCTIONS FOR THE TEXT EDITOR ***************
def merge_files(file1_path, file2_path, common_column):
    try:
        #Load the files
        data1 = pd.read_csv(file1_path, sep=' ')
        data2 = pd.read_csv(file2_path, sep=' ')

        #Merge the files with the common row
        merged_data = pd.merge(data1, data2, on=common_column)

        #Saving the new file
        merged_file_path = sg.popup_get_file('Save the merged file', save_as=True, default_extension=".txt", file_types=(("Text Files", "*.txt"),))
        if merged_file_path:
            merged_data.to_csv(merged_file_path, sep=' ', index=False)
            sg.popup(f'The files have been merged and saved to {merged_file_path}.')
    except Exception as e:
        sg.popup_error(f'Error merging files: {str(e)}')



#**************** FUNCTIONS FOR THE 2D SPECTRA EXTRACTION **********

def open_fits(file_path):
    # Function to open FITS file and return 2D spectrum
    hdu = fits.open(file_path)
    spectrum = hdu[0].data
    header = hdu[0].header
    hdu.close()
    return spectrum, header

def find_and_fit_spectroscopic_trace(spectrum, y_range, poly_degree, first_iteration, with_plots):
    # Function to automatically find and fit spectroscopic trace
    x_axis = np.arange(len(spectrum[0]))
    y_axis = np.arange(len(spectrum))

    spectrum_subset = spectrum
    selected_y_range = slice(*y_range)
    spectrum_subset = spectrum[selected_y_range,:]

    #template = np.sum(spectrum[selected_y_range, :], axis=0)
    template = np.sum(spectrum[selected_y_range, :], axis=0)

    # Compute cross-correlation between each row and the template
    cross_correlation =  correlate2d(spectrum_subset, template[np.newaxis, :], mode='same')
    #cross_correlation = correlate2d(spectrum, template[np.newaxis, :], mode='same')

    # Find the row with maximum cross-correlation for each column
    trace_rows = np.argmax(cross_correlation, axis=0)

    # Fit a polynomial curve to the trace rows along the x-axis using numpy.polyfit
    coefficients = np.polyfit(x_axis, trace_rows, deg=poly_degree)
    trace_model = np.poly1d(coefficients)

    if with_plots == True:
        # Plot the cross-correlation and fitted trace for visualization
        plt.subplot(2, 1, 1)
        plt.imshow(cross_correlation, cmap='gray', aspect='auto')#, norm=LogNorm())
        plt.plot(trace_rows, color='r', linestyle='--', label='Trace Rows')
        plt.legend()
        plt.title("Cross-Correlation for Trace Detection")

        plt.subplot(2, 1, 2)
        plt.plot(x_axis, trace_rows, label="Original Trace Rows")
        plt.plot(x_axis, trace_model(x_axis), label="Trace Model")
        plt.legend()
        plt.title(f"Spectroscopic Trace Model (Degree {poly_degree})")

        plt.tight_layout()
        plt.show()
        plt.close()

    return trace_model

def get_wavelength_coordinates(header, x_axis):
    # Function to get wavelength coordinates from FITS header
    if 'CRVAL1' in header and 'CDELT1' in header:
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        wavelength_coordinates = crval1 + cdelt1 * x_axis
        return wavelength_coordinates
    else:
        return x_axis

def correct_distortion_slope(spectrum, trace_model, y_range):
    # Function to correct distortion and slope of 2D spectrum using the trace model
    x_axis = np.arange(len(spectrum[0]))
    y_axis = np.arange(len(spectrum))
    corrected_spectrum = np.zeros_like(spectrum)

    n_iter = 8
    h = 0

        # Compute correction factor for each column
    correction_factor = trace_model(x_axis) - np.mean(trace_model(x_axis))
    correction_delta = abs(np.min(correction_factor)-np.max(correction_factor))

    if correction_delta > 0.5:
        print ('Iterating the trace fit until the residuals are lower than 0.5 pix')
        print ('Residuals:')
        while correction_delta > 0.5 and h < n_iter:
            for j in range(len(x_axis)):
                # Use interpolation to shift the y-coordinates
                shift_amount = correction_factor[j]
                corrected_spectrum[:, j] = np.interp(y_axis + shift_amount, y_axis, spectrum[:, j], left=np.nan, right=np.nan)

            # Replace NaN values with zeros
            corrected_spectrum = np.nan_to_num(corrected_spectrum)

            #y_range_new = (
            trace_model = find_and_fit_spectroscopic_trace(corrected_spectrum, y_range, 1, False, False)
            spectrum = corrected_spectrum
            correction_factor = trace_model(x_axis) - np.mean(trace_model(x_axis))
            correction_delta = abs(np.min(correction_factor)-np.max(correction_factor))
            h= h+1

            print (correction_delta)

            if correction_delta < 0.5 and h<= n_iter:
                print ('Trace fitting convergence reached')
            if h == n_iter and correction_delta > 0.5:
                print ('Cannot reach the convergence of the fit after ', h, ' iterations. The spectrum could be distorted')

    if correction_delta < 0.5:
        for j in range(len(x_axis)):
            # Use interpolation to shift the y-coordinates
            shift_amount = correction_factor[j]
            corrected_spectrum[:, j] = np.interp(y_axis + shift_amount, y_axis, spectrum[:, j], left=np.nan, right=np.nan)

    # Plot the corrected 2D Spectrum
    plt.imshow(corrected_spectrum, cmap="gray", norm=LogNorm())
    plt.title("Corrected 2D Spectrum")
    plt.show()
    plt.close()

    return corrected_spectrum



def extract_1d_spectrum(corrected_spectrum, y_range, header, x_axis, output_fits_path=None):
    # Function to extract 1D spectrum along x-axis with user-defined y range

    selected_y_range = slice(*y_range)

    extracted_spectrum = np.sum(corrected_spectrum[selected_y_range, :], axis=0)

    # Get wavelength coordinates if available, otherwise use x coordinates
    x_coordinates = get_wavelength_coordinates(header, x_axis)

    # Plot the extracted 1D Spectrum
    plt.plot(x_coordinates, extracted_spectrum)
    plt.title("Extracted 1D Spectrum")
    plt.show()
    plt.close()

    if output_fits_path is not None:
        # Save the extracted 1D spectrum in a FITS file
        hdu = fits.PrimaryHDU(extracted_spectrum)
        hdu.header["CTYPE1"] = 'LINEAR'  # Linear wavelength spacing
        hdu.header["CRPIX1"] = 1  # Reference pixel is the first pixel
        hdu.header["CRVAL1"] = x_coordinates[0]  # Reference value is the first wavelength
        hdu.header["CDELT1"] = np.mean(np.diff(x_coordinates))  # Average wavelength interval

        fits.writeto(output_fits_path, hdu.data, hdu.header, overwrite=True)




def estimate_noise_level(corrected_spectrum, y_range):
    # Function to estimate noise level along the X-axis
    start, end = map(int, y_range)
    selected_y_range = slice(start, end)
    #noise_level_1 = abs(np.nanmean(corrected_spectrum[selected_y_range, :], axis=0))
    noise_level_1 = abs(np.nanmean(corrected_spectrum[selected_y_range, :], axis=0))
    #noise_level_1 = np.nanstd(corrected_spectrum[selected_y_range, :], axis=0)
    noise_level = np.mean(noise_level_1)
    return noise_level



def calculate_signal_to_noise(spectrum_row, noise_level):
    spectrum_row = abs(spectrum_row)
    signal_to_noise = (np.nanmean(spectrum_row)/noise_level)
    return signal_to_noise


# WITH TWO NOISE REGIONS TO SELECT. THE NOISE WILL BE A SIMPLE MEAN OF THE TWO REGIONS. COMMENT THIS AND UNCOMMENT THE FOLLOWING VERSION IF YOU PREFER TO SELECT ONLY ONE NOISE REGION.
def extract_and_save_snr_spectra(corrected_spectrum, trace_model, header, x_axis, snr_threshold, pixel_scale, file_path):
    # Function to extract and save 1D spectra based on the mean signal-to-noise threshold along the X-axis
    y_axis = np.arange(len(corrected_spectrum))
    n_rows = len(y_axis)
    spectra = []

    signal_profile = np.sum(corrected_spectrum, axis=1)
    print ('Please, select two regions containing noise. Two clicks for each region: one for the start and the other for the end')
    print ('WARNING: Do not close the plot window without selecting the noise regions. Otherwise the program will freeze')
    # Allow the user to click on the corrected spectrum to select the Y-values for noise estimation
    plt.plot(y_axis, signal_profile)
    plt.title("Corrected 2D Spectrum - Select Noise Region")
    plt.xlabel("Y-axis")
    plt.ylabel("Intensity")
    noise_region_points = plt.ginput(n=4, timeout=-1, show_clicks=True)
    noise_region_y_values_1 = [int(min(noise_region_points[0][0], noise_region_points[1][0])),
                                int(max(noise_region_points[0][0], noise_region_points[1][0]))]
    noise_region_y_values_2 = [int(min(noise_region_points[2][0], noise_region_points[3][0])),
                                int(max(noise_region_points[2][0], noise_region_points[3][0]))]

    plt.close()
    # Calculate mean noise level for the two regions
    noise_level_1 = estimate_noise_level(corrected_spectrum, noise_region_y_values_1)
    noise_level_2 = estimate_noise_level(corrected_spectrum, noise_region_y_values_2)
    noise_level = np.mean([noise_level_1, noise_level_2])

    y_positions = []
    y_positions_mean = []

    trace_mean_y_position = round(int(np.mean(trace_model(y_axis))))

    print ('')
    print ('Selected noise regions', noise_region_y_values_1, noise_region_y_values_2)
    print ('')
    print ('Mean noise level', noise_level)
    print ('')

    y_pos = 0
    i = 0
    while i < n_rows-1:
        if (min(noise_region_y_values_1) < np.min(y_axis) or max(noise_region_y_values_1) > np.max(y_axis) or min(noise_region_y_values_2) < np.min(y_axis) or max(noise_region_y_values_2) > np.max(y_axis)):
            sg.popup ('Noise region outside the spectrum!')
            break

        snr = calculate_signal_to_noise(corrected_spectrum[i, :], noise_level)
        if snr >= snr_threshold:
            # If the current row meets the threshold, add it to the spectra
            spectra.append(corrected_spectrum[i, :])
            y_pos = i
            y_positions.append(y_pos)
            y_positions_mean.append(y_pos)
            i += 1
        else:
            # If the current row does not meet the threshold, sum consecutive rows until the threshold is reached
            #The y position will be the snr weigthed mean position of the created bin.
            snr_for_mean = []
            y_for_mean = []
            summed_spectrum = np.copy(corrected_spectrum[i, :])
            n_bins = 1
            while i + 1 < n_rows and np.nanmean(snr) < snr_threshold:
                i += 1
                n_bins +=1 #IT IS CORRECT HERE? I AM LOOSING A ROW?
                snr_single = calculate_signal_to_noise(corrected_spectrum[i, :], noise_level)
                snr_for_mean.append(snr_single)
                y_for_mean.append(i)
                #i += 1
                summed_spectrum += corrected_spectrum[i, :]
                mean_spectrum = summed_spectrum/n_bins
                noise_level_new = noise_level/mt.sqrt(n_bins)
                snr = calculate_signal_to_noise(mean_spectrum, noise_level_new)

            y_pos_mean_var = 0
            for t in range (len(y_for_mean)):
                y_pos_mean_var += (y_for_mean[t]*snr_for_mean[t])

            y_pos_mean_var = y_pos_mean_var/sum(snr_for_mean)
            y_positions_mean.append(round(y_pos_mean_var,1))
            y_pos = i
            y_positions.append(y_pos)
            spectra.append(summed_spectrum)

    if (min(noise_region_y_values_1) < np.min(y_axis) or max(noise_region_y_values_1) > np.max(y_axis) or min(noise_region_y_values_2) < np.min(y_axis) or max(noise_region_y_values_2) > np.max(y_axis)):
        print ('No files saved')
    else:

        spectra = np.array(spectra)
        y_positions = np.array(y_positions)
        y_positions_mean = np.array(y_positions_mean)

        if pixel_scale != 0:
            arcsec_scale = y_positions*pixel_scale
            arcsec_scale_mean = y_positions_mean*pixel_scale

        print ('Spectral bins Y mean position: ', y_positions_mean)
        print ('')
        print ('Number of bins: ', len(y_positions_mean))

        # Get wavelength coordinates if available, otherwise use x coordinates
        x_coordinates = get_wavelength_coordinates(header, x_axis)

        # Save the 1D spectra in IRAF-style FITS format
        for i, spectrum_row in enumerate(spectra):
            hdu = fits.PrimaryHDU(spectrum_row)
            hdu.header["CTYPE1"] = 'LINEAR'  # Linear wavelength spacing
            hdu.header["CRPIX1"] = 1  # Reference pixel is the first pixel
            hdu.header["CRVAL1"] = x_coordinates[0]  # Reference value is the first wavelength
            hdu.header["CDELT1"] = np.mean(np.diff(x_coordinates))  # Average wavelength interval

            #adding the position of the 1d spectra with respect to the central trace to the fits header
            hdu.header.set("Y_POS", y_positions_mean[i] - trace_mean_y_position, "Pix position from the center")
            if pixel_scale != 0:
                hdu.header.set("R", arcsec_scale_mean[i] - (trace_mean_y_position*pixel_scale), "Arcsec position from the center")
            hdu.writeto(f"spectrum_snr_{i+1}.fits", overwrite=True)

        # Plot the extracted 1D Spectra
        plt.figure()
        for i, spectrum_row in enumerate(spectra):
            plt.plot(x_coordinates, spectrum_row, label=f"Spectrum {i+1}")

        plt.xlabel("Wavelength")
        plt.ylabel("Flux")
        plt.legend()
        plt.title("Extracted 1D SNR Spectra")
        plt.show()
        plt.close()
        sg.popup ('1D spectra saved in the working directory')



#UNCOMMENT THE FOLLOWING FUNCTION AND COMMENT THE PREVIOUS IF YOU WANT TO SELECT JUST ONE NOISE REGION INSTEAD OF TWO!
#def extract_and_save_snr_spectra(corrected_spectrum, trace_model, header, x_axis, snr_threshold, pixel_scale, file_path):
    ## Function to extract and save 1D spectra based on the mean signal-to-noise threshold along the X-axis
    #y_axis = np.arange(len(corrected_spectrum))
    #n_rows = len(y_axis)
    #spectra = []

    #signal_profile = np.sum(corrected_spectrum, axis=1)
    #print ('Please, select a region containing noise. Two clicks: one for the start and the other for the end')
    #print ('WARNING: Do not close the plot window without selecting the noise regions. Otherwise the program will freeze')
    ## Allow the user to click on the corrected spectrum to select the Y-values for noise estimation
    #plt.plot(y_axis, signal_profile)
    #plt.title("Corrected 2D Spectrum - Select Noise Region")
    #plt.xlabel("Y-axis")
    #plt.ylabel("Intensity")
    #noise_region_points = plt.ginput(n=2, timeout=-1, show_clicks=True)
    #noise_region_y_range = (int(min(noise_region_points[0][0], noise_region_points[1][0])),
                            #int(max(noise_region_points[0][0], noise_region_points[1][0])))
    #plt.close()

    #noise_level = estimate_noise_level(corrected_spectrum, noise_region_y_range)

    #y_positions = []
    #y_positions_mean = []

    #trace_mean_y_position = round(int(np.mean(trace_model(y_axis))))

    #print ('')
    #print ('Selected noise region', noise_region_y_range)
    #print ('')
    #print ('Mean noise level', noise_level)
    #print ('')
    ##print (len(corrected_spectrum))

    #y_pos = 0
    #i = 0
    #while i < n_rows-1:
        #if (min(noise_region_y_range) < np.min(y_axis) or max(noise_region_y_range) > np.max(y_axis)):
            #sg.popup ('Noise region outside the spectrum!')
            #break

        #snr = calculate_signal_to_noise(corrected_spectrum[i, :], noise_level)
        #if np.nanmean(snr) >= snr_threshold:
            ## If the current row meets the threshold, add it to the spectra
            #spectra.append(corrected_spectrum[i, :])
            #y_pos = i
            #y_positions.append(y_pos)
            #y_positions_mean.append(y_pos)
            #i += 1
        #else:
            ## If the current row does not meet the threshold, sum consecutive rows until the threshold is reached
            #snr_for_mean = []
            #y_for_mean = []
            #summed_spectrum = np.copy(corrected_spectrum[i, :])
            #n_bins = 0
            #while i + 1 < n_rows and np.nanmean(snr) < snr_threshold:
                #i += 1
                #n_bins +=1
                #snr_single = calculate_signal_to_noise(corrected_spectrum[i, :], noise_level)
                #snr_for_mean.append(snr_single)
                #y_for_mean.append(i)
                #summed_spectrum += corrected_spectrum[i, :]
                #mean_spectrum = summed_spectrum/n_bins
                #noise_level_new = noise_level/mt.sqrt(n_bins)
                #snr = calculate_signal_to_noise(mean_spectrum, noise_level_new)

            #y_pos_mean_var = 0
            #for t in range (len(y_for_mean)):
                #y_pos_mean_var += (y_for_mean[t]*snr_for_mean[t])

            #y_pos_mean_var = y_pos_mean_var/sum(snr_for_mean)
            #y_positions_mean.append(round(y_pos_mean_var,1))
            #y_pos = i
            #y_positions.append(y_pos)
            #spectra.append(summed_spectrum)

    #if (min(noise_region_y_range) < np.min(y_axis) or max(noise_region_y_range) > np.max(y_axis)):
        #print ('No files saved')
    #else:

        #spectra = np.array(spectra)
        #y_positions = np.array(y_positions)
        #y_positions_mean = np.array(y_positions_mean)

        #if pixel_scale != 0:
            #arcsec_scale = y_positions*pixel_scale
            #arcsec_scale_mean = y_positions_mean*pixel_scale

        ##print ('Spectral bins Y position: ', y_positions)
        #print ('Spectral bins Y mean position: ', y_positions_mean)
        #print ('')
        #print ('Number of bins: ', len(y_positions_mean))

        ## Get wavelength coordinates if available, otherwise use x coordinates
        #x_coordinates = get_wavelength_coordinates(header, x_axis)

        ## Save the 1D spectra in IRAF-style FITS format
        #for i, spectrum_row in enumerate(spectra):
            #hdu = fits.PrimaryHDU(spectrum_row)
            #hdu.header["CTYPE1"] = 'LINEAR'  # Linear wavelength spacing
            #hdu.header["CRPIX1"] = 1  # Reference pixel is the first pixel
            #hdu.header["CRVAL1"] = x_coordinates[0]  # Reference value is the first wavelength
            #hdu.header["CDELT1"] = np.mean(np.diff(x_coordinates))  # Average wavelength interval

            ##adding the position of the 1d spectra with respect to the central trace to the fits header
            #hdu.header.set("Y_POS", y_positions_mean[i] - trace_mean_y_position, "Pix position from the center")
            #if pixel_scale != 0:
                #hdu.header.set("R", arcsec_scale_mean[i] - (trace_mean_y_position*pixel_scale), "Arcsec position from the center")

            #base_filename = os.path.splitext(os.path.basename(file_path))[0]
            #hdu.writeto(f"{base_filename}_{i+1}.fits", overwrite=True)

        ## Plot the extracted 1D Spectra
        #plt.figure()
        #for i, spectrum_row in enumerate(spectra):
            #plt.plot(x_coordinates, spectrum_row, label=f"Spectrum {i+1}")

        #plt.xlabel("Wavelength")
        #plt.ylabel("Flux")
        #plt.legend()
        #plt.title("Extracted 1D SNR Spectra")
        #plt.show()
        #plt.close()
        #sg.popup ('1D spectra saved in the working directory')


# Correction EW for the sigma broadening of Lick/IDS indices
def corr_ew_lick(ew_values, err_values, ew_mag_values, coeff_file, sigma_value):

    #Converting the EW and EW_err to numpy
    ew_values_np = ew_values #ew_values.to_numpy(dtype = float)
    err_values_np = err_values #err_values.to_numpy(dtype = float)

    #Open and reading the file containing the coefficients
    data_corr = pd.read_csv(coeff_file, header=None, sep = ' ')
    data_number = len(ew_values)

    corr_ew_values = data_corr.iloc[1:, :].values
    corr_ew_values_np = corr_ew_values.astype(float)

    #creating the array of the corrected ews and errors
    ew_correction = np.zeros(data_number)
    new_ew_err = np.zeros(data_number)
    new_ew = np.zeros(data_number)


    # CORRECTING THE EWS AND ERRORS, CONSIDERING BOTH IN ANGSTROM AND MAGNITUDES

    for i in range(data_number):
        #correcting ews
        ew_correction[i] = (corr_ew_values_np[3,i] + corr_ew_values_np[2,i]*sigma_value+corr_ew_values_np[1,i]*sigma_value**2+corr_ew_values_np[0,i]*sigma_value**3)
        new_ew[i] = ew_values_np[i]*ew_correction[i]

    new_ew_err = err_values_np*ew_correction
    new_ew_mag = ew_mag_values*ew_correction
    new_ew_err_mag = np.zeros(len(new_ew_err))
    for i in range(data_number):
        if new_ew[i] == 0:
            new_ew[i] = 0
        else:
            new_ew_err_mag[i]= 0.434*abs(new_ew_err[i]/new_ew[i])
        #new_ew_err_mag[i] = 0.

    return new_ew, new_ew_err, new_ew_mag, new_ew_err_mag


# TEXT EDITOR FUNCTIONS
def save_file(filename, text):
    with open(filename, 'w') as file:
        file.write(text)

def undo(text_history, text_index, window_editor):
    if text_index > 0:
        text_index -= 1
        window_editor['-TEXT-'].update(text_history[text_index])
    return text_index

def find_replace(text, find, replace, replace_all):
    if replace_all:
        return text.replace(find, replace)
    else:
        return text.replace(find, replace, 1)

def create_new_column(df, new_column_name, col1_name, col2_name, expression):
    try:
        df[new_column_name] = pd.eval(expression, engine='python',
                                    local_dict={col1_name: df[col1_name], col2_name: df[col2_name]})
        return df
    except Exception as e:
        sg.popup_error(f'Error creating the new column: {str(e)}')
        return None


# functions for spectra list

def get_files_in_folder(folder_path):
    file_list = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            absolute_path = os.path.join(root, file)
            file_list.append(absolute_path)
    return sorted(file_list, key=str.lower)

def save_to_text_file(file_list, output_file):
    with open(output_file, 'w') as f:
        f.write("#Spectrum\n")
        for absolute_path in file_list:
            f.write(f"{absolute_path}\n")


###############################################################################

#simple cropping function
def crop_spec(wavelength, flux, wave_interval):
    wave1 = np.min(wave_interval)
    wave2 = np.max(wave_interval)
    cropped_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]
    cropped_flux = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    return cropped_wave, cropped_flux

#wavelets denoising
def wavelet_cleaning(wavelength, flux, sigma, wavelets_layers):
    #normalizing the spectrum
    epsilon_norm = 2
    #norm_flux = norm_spec(wavelength, flux, np.mean(wavelength), epsilon_norm, flux)
    #flux = norm_flux
    denoised_flux = denoise_wavelet(flux, sigma=sigma, wavelet='sym5', mode='soft', wavelet_levels=wavelets_layers)

    return denoised_flux


def degradeRtoFWHM(wavelength, flux, R_resolution, FWHM_resolution):

    fwhm_to_sigma = 2.3548
    #original_resolution_lambda_nm = original_resolution/10. #converting to nm!
    FWHM_resolution = FWHM_resolution/10. #converting to nm!

    R_resolution_to_FWHM_nm = wavelength/R_resolution #array contenente le risoluzioni in FWHM
    final_resolution_FWHM_nm = np.full_like(R_resolution_to_FWHM_nm, FWHM_resolution, dtype=float)

    real_value_to_broad = np.zeros_like(wavelength)
    degraded_flux = np.zeros_like(flux)

    try:
        for i in range (len(wavelength)):

            real_value_to_broad[i] = mt.sqrt(final_resolution_FWHM_nm[i]**2-R_resolution_to_FWHM_nm[i]**2)

        gauss_sigma = real_value_to_broad/fwhm_to_sigma

        #using the varsmooth function of ppxf.util for convolution with variable sigma. Works great!
        degraded_flux = util.varsmooth(wavelength, flux, gauss_sigma, xout=None, oversample=1)
        return wavelength, degraded_flux
    except Exception:
        print ('You want to improve the resolution? That''s impossible! Skypping...')
        return wavelength, flux

#def degradeFWHMtoR(wavelength, flux, FWHM_resolution, R_resolution):

    #fwhm_to_sigma = 2.3548
    ##original_resolution_lambda_nm = original_resolution/10. #converting to nm!
    #FWHM_resolution = FWHM_resolution/10. #converting to nm!

    #FWHM_resolution_to_R_nm = wavelength/FWHM_resolution
    #final_resolution_R_nm = np.full_like(FWHM_resolution_to_R_nm, R_resolution, dtype=float)

def mask_spectrum(wavelength, flux, mask_ranges):
    mask = np.ones_like(wavelength, dtype=bool)

    for mask_range in mask_ranges:
        mask = mask & ((wavelength < mask_range[0]) | (wavelength > mask_range[1]))

    return mask

def continuum(wavelength, flux, want_to_maks, mask_ranges, poly_degree, math_operation, with_plots):

    if want_to_maks  == True:
        mask = mask_spectrum(wavelength, flux, mask_ranges)
        wavelength_masked = wavelength[mask]
        flux_masked = flux[mask]
    else:
        wavelength_masked = wavelength
        flux_masked = flux

    #poly_degree = 5

    coefficients = np.polyfit(wavelength_masked, flux_masked, poly_degree)
    continuum_model = np.polyval(coefficients, wavelength)

    if math_operation == 'subtract':
        new_flux = flux - continuum_model
    elif math_operation == 'divide':
        new_flux = flux/continuum_model

    if with_plots == True:
        # Visualizza i risultati con le regioni mascherate
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot dello spettro originale
        ax.plot(wavelength, flux, label='Original spectrum')

        # Plot del modello del continuo
        ax.plot(wavelength, continuum_model, label='Continuum model')

        # Plot dello spettro corretto
        ax.plot(wavelength, new_flux, label='Corrected spectrum')

        # Evidenzia le regioni mascherate utilizzando axvspan
        if want_to_maks == True:
            for mask_range in mask_ranges:
                ax.axvspan(mask_range[0], mask_range[1], color='gray', alpha=0.5)

        ax.legend()
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        plt.title('Continuum fitting')
        plt.show()
        plt.close()
    return new_flux, continuum_model

#lowpass
def lowpass(wavelength, flux, cut_off, order):
    try:
        b, a = butter(order, cut_off, btype='lowpass', analog=False)
        denoised_flux = filtfilt(b, a, flux)
        return denoised_flux
    except:
        print ('Error applying the lowpass filter. Skipping...')
        return flux

def mov_avg_gauss(wavelength, flux, sigma):
    try:
        denoised_flux = gaussian_filter1d(flux, sigma=sigma)
        return denoised_flux
    except:
        print ('Error applying the gaussian moving average. Skipping...')
        return flux

#bandpass
def bandpass(wavelength, flux, lower_cut_off, upper_cut_off, order):
    try:
        b, a = butter(order, [lower_cut_off, upper_cut_off], btype='band', analog=False)
        denoised_flux = filtfilt(b, a, flux)
        return denoised_flux
    except:
        print ('Error applying the bandpass filter. Skipping...')
        return flux

