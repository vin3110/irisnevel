#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:39:57 2024

@author: heikamp


start up script for astronomy image analysis with python

"""
#import packages

from astropy.io import fits #to read in FITS files
import os # os.path to manipulate file paths 
import glob # finding pathnames (to search for certain fits files in folders and subfolders)
import numpy as np # math applied to arrays (important, no need to read pixel for pixel!)
from matplotlib import pyplot as plt #plot library 
from astropy.visualization import ZScaleInterval #create minimum and maximum Z values for plotting 
from photutils.aperture import aperture_photometry, CircularAperture
from photutils.aperture import CircularAnnulus, ApertureStats
import math
from skimage.registration import phase_cross_correlation
from scipy.ndimage import fourier_shift

def brightness_i(data_im):

    hdu = fits.open(data_im)[0]
    header = hdu.header 
    data = hdu.data

    data = data[:4046,:4046]

    #print(header)

    z = ZScaleInterval()
    z1,z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')


    #import list of fits files:
    # set path to calibration folder (path_in) and path to save files (path_out)
    # path_in = " /Users/Documents/DATA/"
    # path_out=" /Users/Documents/DATA/python/"

    # file_list = glob.glob(path_in+'/**/'+str("*.fit*"),recursive=True) #/**/ means to search subfolders if recursive = True
    # sorted_list = sorted(file_list, key=os.path.getmtime) 

    # hdu_list = fits.open(sorted_list[10])[0] # 10 = 10th item in list
    # header = hdu_list.header 
    # data = hdu_list.data


    #plate solve image: upload image to nova.astrometry.net/upload


    #simple manipulation and save fits file
    # data = data.T 
    # data_sub = data_frame1 - data_frame2
    # data_median = np.median(data)


    # hdu = fits.PrimaryHDU(data_sub, header=header)
    # hdu.writeto(path_out + "filename.fits", overwrite=True)

    def instrumental_mag(x_pixel, y_pixel, radius, r_in, r_out):
        position = (x_pixel, y_pixel) # Kan meerdere tegelijk met lijst
        aperture = CircularAperture(position, r= radius)

        phot_table = aperture_photometry(data, aperture)
        # aperature_sum column is the sum of the counts

        annulus_aperature = CircularAnnulus(position, r_in, r_out) #r_in moet groter dan radius zijn
        aperstats = ApertureStats(data, annulus_aperature)
        bkg_mean = aperstats.mean
        total_bkg = bkg_mean * aperture.area

        phot_bkgsub = phot_table['aperture_sum'] - total_bkg

        # print(phot_bkgsub)
        # print(total_bkg)

        im = -2.5*math.log10(phot_bkgsub) # ook deze funcite voor im van nevel
        err_zeropoint = 100 * (np.sqrt((phot_table['aperture_sum'] - total_bkg) + aperture.area * (bkg_mean + 4.2 + 100))) / total_bkg

        return im, bkg_mean, err_zeropoint


    # UCAC4 791-034810
    abs_mag_i_uca = 12.768
    x_pixel_uca = 2782
    y_pixel_uca = 1688
    r_uca = 18
    r_in_uca = 22
    r_out_uca = 30


    # IRAS 21023+6754
    abs_mag_i_iras = 15.02
    x_pixel_iras = 2545
    y_pixel_iras = 3183
    r_iras = 12
    r_in_iras = 16
    r_out_iras = 24
    # abs_mag_I_iras = 14.46
    # abs_mag_R_iras = 15.69
    # abs_mag_G_iras = 15.619860

    # EM* LkHA 428B
    abs_mag_I_EM = 13.95
    abs_mag_R_EM = 15.86
    abs_mag_G_EM = 15.502530

    # Cl* NGC 7023 RS 5
    abs_mag_I_Cl = 14.35
    abs_mag_R_Cl = 16.28

    # def ratios(magI, magR, magG):

    UCA_im, achtergrond_uca, err_zpt_UCA = instrumental_mag(x_pixel_uca, y_pixel_uca, r_uca, r_in_uca, r_out_uca)
    IRAS_im, achtergrond_iras, err_zpt_IRAS = instrumental_mag(x_pixel_iras, y_pixel_iras, r_iras, r_in_iras, r_out_iras)
    print(UCA_im)
    print(IRAS_im)

    delta_m_uca = abs_mag_i_uca - UCA_im
    delta_m_iras = abs_mag_i_iras - IRAS_im
    delta_m = (delta_m_iras+delta_m_uca) / 2

    achtergrond = (achtergrond_iras + achtergrond_uca) / 2
    # print(delta_m_iras)
    print(delta_m)

    # in g:
    # 1708,1716 bij 2565,2252
    # maskeren: 2028,2042 bij 2097,2100
    # g naar i: -15, +58

    # Afgekaderd gebied
    # x_start, x_end = 2139, 2230
    # y_start, y_end = 2005, 2182

    x_start, x_end = 1692, 2550
    y_start, y_end = 1774, 2310
    sub_data = data[y_start:y_end, x_start:x_end]

    # Define the mask region within the sub-region
    x_mask_start, x_mask_end = 2013 - x_start, 2082 - x_start
    y_mask_start, y_mask_end = 2100 - y_start, 2158 - y_start

    # # Apply the mask (setting the specified region to zero)
    sub_data[y_mask_start:y_mask_end, x_mask_start:x_mask_end] = 0

    # data_kader = data[y_start:y_end, x_start:x_end]

    z = ZScaleInterval()
    z1,z2 = z.get_limits(sub_data)
    plt.imshow(sub_data, vmin=z1, vmax=z2, cmap='gray')
    plt.colorbar()
    #plt.show()

    nevel_zonder_achtergrond = sub_data - achtergrond
    im_nevel = -2.5 * np.log10(nevel_zonder_achtergrond) 
    m_nevel = im_nevel + delta_m 


    # Calculate the total area of the sub-region
    total_area = (x_end - x_start)*0.4310157882516833 * (y_end - y_start)*0.4310157882516833

    # Calculate the area of the mask region
    mask_area = (x_mask_end - x_mask_start)*0.4310157882516833 * (y_mask_end - y_mask_start)*0.4310157882516833

    # Calculate the area of the unmasked region
    unmasked_area = (total_area - mask_area)

    surface_brightness = m_nevel + 2.5 * np.log10(0.4310157882516833**2)

    # z = ZScaleInterval()
    # z1,z2 = z.get_limits(surface_brightness)
    plt.imshow(surface_brightness, cmap='gray')
    plt.colorbar()
    #plt.show()

    print(err_zpt_UCA)
    print(err_zpt_IRAS)
    err_zpt = (err_zpt_UCA + err_zpt_IRAS) / 2
    print(err_zpt)
    lijn = nevel_zonder_achtergrond[214:403, 513]
    err_mag = np.sqrt(((-2.5 / np.log(10)) * (1 / np.sqrt(lijn))) ** 2 + err_zpt**2)

    print(err_mag)

    # Save the masked data to a new FITS file
    # masked_hdu = fits.PrimaryHDU(sub_data, header=header)
    # masked_hdu.writeto('masked_data.fits', overwrite=True)

    return surface_brightness, sub_data, im_nevel, err_mag



def brightness_g(data_im):

    hdu = fits.open(data_im)[0]
    header = hdu.header 
    data = hdu.data

    data = data[:4046,:4046]

    #print(header)

    z = ZScaleInterval()
    z1,z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')


    #import list of fits files:
    # set path to calibration folder (path_in) and path to save files (path_out)
    # path_in = " /Users/Documents/DATA/"
    # path_out=" /Users/Documents/DATA/python/"

    # file_list = glob.glob(path_in+'/**/'+str("*.fit*"),recursive=True) #/**/ means to search subfolders if recursive = True
    # sorted_list = sorted(file_list, key=os.path.getmtime) 

    # hdu_list = fits.open(sorted_list[10])[0] # 10 = 10th item in list
    # header = hdu_list.header 
    # data = hdu_list.data


    #plate solve image: upload image to nova.astrometry.net/upload


    #simple manipulation and save fits file
    # data = data.T 
    # data_sub = data_frame1 - data_frame2
    # data_median = np.median(data)


    # hdu = fits.PrimaryHDU(data_sub, header=header)
    # hdu.writeto(path_out + "filename.fits", overwrite=True)

    def instrumental_mag(x_pixel, y_pixel, radius, r_in, r_out):
        position = (x_pixel, y_pixel) # Kan meerdere tegelijk met lijst
        aperture = CircularAperture(position, r= radius)

        phot_table = aperture_photometry(data, aperture)
        # aperature_sum column is the sum of the counts

        annulus_aperature = CircularAnnulus(position, r_in, r_out) #r_in moet groter dan radius zijn
        aperstats = ApertureStats(data, annulus_aperature)
        bkg_mean = aperstats.mean
        total_bkg = bkg_mean * aperture.area

        phot_bkgsub = phot_table['aperture_sum'] - total_bkg

        # print(phot_bkgsub)
        # print(total_bkg)

        im = -2.5*math.log10(phot_bkgsub) # ook deze funcite voor im van nevel

        err_zeropoint = 100 * (np.sqrt((phot_table['aperture_sum'] - total_bkg) + aperture.area * (bkg_mean + 4.2 + 100))) / total_bkg

        return im, bkg_mean, err_zeropoint


    # UCAC4 791-034810
    abs_mag_g_uca = 15.128
    # x_pixel_uca = 2797 # Van voor alignment
    # y_pixel_uca = 1631
    x_pixel_uca = 2782
    y_pixel_uca = 1688
    r_uca = 18
    r_in_uca = 22
    r_out_uca = 30


    # IRAS 21023+6754
    abs_mag_g_iras = 16.3196 # berekend met formule
    # x_pixel_iras = 2561
    # y_pixel_iras = 3125
    x_pixel_iras = 2545
    y_pixel_iras = 3183
    r_iras = 12
    r_in_iras = 16
    r_out_iras = 24
    # abs_mag_I_iras = 14.46
    # abs_mag_R_iras = 15.69
    # abs_mag_G_iras = 15.619860

    # EM* LkHA 428B
    abs_mag_I_EM = 13.95
    abs_mag_R_EM = 15.86
    abs_mag_G_EM = 15.502530

    # Cl* NGC 7023 RS 5
    abs_mag_I_Cl = 14.35
    abs_mag_R_Cl = 16.28

    # def ratios(magI, magR, magG):

    UCA_im, achtergrond_uca, err_zpt_UCA = instrumental_mag(x_pixel_uca, y_pixel_uca, r_uca, r_in_uca, r_out_uca)
    IRAS_im, achtergrond_iras, err_zpt_IRAS = instrumental_mag(x_pixel_iras, y_pixel_iras, r_iras, r_in_iras, r_out_iras)
    print(UCA_im)
    print(IRAS_im)

    delta_m_uca = abs_mag_g_uca - UCA_im
    delta_m_iras = abs_mag_g_iras - IRAS_im
    delta_m = (delta_m_iras+delta_m_uca) / 2

    achtergrond = (achtergrond_iras + achtergrond_uca) / 2
    # print(delta_m_iras)
    print(delta_m)

    # in g:
    # 1708,1716 bij 2565,2252
    # maskeren: 2028,2042 bij 2097,2100
    # g naar i: -15, +58

    # Afgekaderd gebied
    # x_start, x_end = 2139, 2230
    # y_start, y_end = 2005, 2182

    # x_start, x_end = 1708, 2565
    # y_start, y_end = 1716, 2252
    # sub_data = data[y_start:y_end, x_start:x_end]

    # # Define the mask region within the sub-region
    # x_mask_start, x_mask_end = 2028 - x_start, 2097 - x_start
    # y_mask_start, y_mask_end = 2042 - y_start, 2100 - y_start

    x_start, x_end = 1692, 2550
    y_start, y_end = 1774, 2310
    sub_data = data[y_start:y_end, x_start:x_end]

    # Define the mask region within the sub-region
    x_mask_start, x_mask_end = 2013 - x_start, 2082 - x_start
    y_mask_start, y_mask_end = 2100 - y_start, 2158 - y_start

    # # Apply the mask (setting the specified region to zero)
    sub_data[y_mask_start:y_mask_end, x_mask_start:x_mask_end] = 0

    # data_kader = data[y_start:y_end, x_start:x_end]

    z = ZScaleInterval()
    z1,z2 = z.get_limits(sub_data)
    plt.imshow(sub_data, vmin=z1, vmax=z2, cmap='gray')
    plt.colorbar()
    #plt.show()

    nevel_zonder_achtergrond = sub_data - achtergrond
    im_nevel = -2.5 * np.log10(nevel_zonder_achtergrond) 
    m_nevel = im_nevel + delta_m 


    # Calculate the total area of the sub-region
    total_area = (x_end - x_start)*0.4310157882516833 * (y_end - y_start)*0.4310157882516833

    # Calculate the area of the mask region
    mask_area = (x_mask_end - x_mask_start)*0.4310157882516833 * (y_mask_end - y_mask_start)*0.4310157882516833

    # Calculate the area of the unmasked region
    unmasked_area = (total_area - mask_area)

    surface_brightness = m_nevel + 2.5 * np.log10(0.4310157882516833**2)

    # z = ZScaleInterval()
    # z1,z2 = z.get_limits(surface_brightness)
    plt.imshow(surface_brightness, cmap='gray')
    plt.colorbar()
    #plt.show()

    err_zpt = (err_zpt_UCA + err_zpt_IRAS) / 2
    lijn = nevel_zonder_achtergrond[214:403, 513]
    err_mag = np.sqrt(((-2.5 / np.log(10)) * (1 / np.sqrt(lijn))) ** 2 + err_zpt**2)

    # Save the masked data to a new FITS file
    # masked_hdu = fits.PrimaryHDU(sub_data, header=header)
    # masked_hdu.writeto('masked_data.fits', overwrite=True)

    return surface_brightness, sub_data, im_nevel, err_mag

def brightness_r(data_im):

    hdu = fits.open(data_im)[0]
    header = hdu.header 
    data = hdu.data

    data = data[:4046,:4046]


    z = ZScaleInterval()
    z1,z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')


    #import list of fits files:
    # set path to calibration folder (path_in) and path to save files (path_out)
    # path_in = " /Users/Documents/DATA/"
    # path_out=" /Users/Documents/DATA/python/"

    # file_list = glob.glob(path_in+'/**/'+str("*.fit*"),recursive=True) #/**/ means to search subfolders if recursive = True
    # sorted_list = sorted(file_list, key=os.path.getmtime) 

    # hdu_list = fits.open(sorted_list[10])[0] # 10 = 10th item in list
    # header = hdu_list.header 
    # data = hdu_list.data


    #plate solve image: upload image to nova.astrometry.net/upload


    #simple manipulation and save fits file
    # data = data.T 
    # data_sub = data_frame1 - data_frame2
    # data_median = np.median(data)


    # hdu = fits.PrimaryHDU(data_sub, header=header)
    # hdu.writeto(path_out + "filename.fits", overwrite=True)

    def instrumental_mag(x_pixel, y_pixel, radius, r_in, r_out):
        position = (x_pixel, y_pixel) # Kan meerdere tegelijk met lijst
        aperture = CircularAperture(position, r= radius)

        phot_table = aperture_photometry(data, aperture)
        # aperature_sum column is the sum of the counts

        annulus_aperature = CircularAnnulus(position, r_in, r_out) #r_in moet groter dan radius zijn
        aperstats = ApertureStats(data, annulus_aperature)
        bkg_mean = aperstats.mean
        total_bkg = bkg_mean * aperture.area

        phot_bkgsub = phot_table['aperture_sum'] - total_bkg

        # print(phot_bkgsub)
        # print(total_bkg)

        im = -2.5*math.log10(phot_bkgsub) # ook deze funcite voor im van nevel
        

        return im, bkg_mean,


    # UCAC4 791-034810
    abs_mag_r_uca = 13.700
    x_pixel_uca = 2782
    y_pixel_uca = 1688
    r_uca = 18
    r_in_uca = 22
    r_out_uca = 30


    # IRAS 21023+6754
    abs_mag_r_iras = 15.5806 # Berekend met formule
    x_pixel_iras = 2545
    y_pixel_iras = 3183
    r_iras = 12
    r_in_iras = 16
    r_out_iras = 24
    # abs_mag_I_iras = 14.46
    # abs_mag_R_iras = 15.69
    # abs_mag_G_iras = 15.619860

    # EM* LkHA 428B
    abs_mag_I_EM = 13.95
    abs_mag_R_EM = 15.86
    abs_mag_G_EM = 15.502530

    # Cl* NGC 7023 RS 5
    abs_mag_I_Cl = 14.35
    abs_mag_R_Cl = 16.28

    # def ratios(magI, magR, magG):

    UCA_im, achtergrond_uca = instrumental_mag(x_pixel_uca, y_pixel_uca, r_uca, r_in_uca, r_out_uca)
    IRAS_im, achtergrond_iras = instrumental_mag(x_pixel_iras, y_pixel_iras, r_iras, r_in_iras, r_out_iras)
    print(UCA_im)
    print(IRAS_im)

    delta_m_uca = abs_mag_r_uca - UCA_im
    delta_m_iras = abs_mag_r_iras - IRAS_im
    delta_m = (delta_m_iras+delta_m_uca) / 2

    achtergrond = (achtergrond_iras + achtergrond_uca) / 2
    # print(delta_m_iras)
    print(delta_m)

    # in g:
    # 1708,1716 bij 2565,2252
    # maskeren: 2028,2042 bij 2097,2100
    # g naar i: -15, +58

    # Afgekaderd gebied
    # x_start, x_end = 2139, 2230
    # y_start, y_end = 2005, 2182

    x_start, x_end = 1692, 2550
    y_start, y_end = 1774, 2310
    sub_data = data[y_start:y_end, x_start:x_end]

    # Define the mask region within the sub-region
    x_mask_start, x_mask_end = 2013 - x_start, 2082 - x_start
    y_mask_start, y_mask_end = 2100 - y_start, 2158 - y_start

    # # Apply the mask (setting the specified region to zero)
    sub_data[y_mask_start:y_mask_end, x_mask_start:x_mask_end] = 0

    # data_kader = data[y_start:y_end, x_start:x_end]

    z = ZScaleInterval()
    z1,z2 = z.get_limits(sub_data)
    plt.imshow(sub_data, vmin=z1, vmax=z2, cmap='gray')
    plt.colorbar()
    #plt.show()

    nevel_zonder_achtergrond = sub_data - achtergrond
    im_nevel = -2.5 * np.log10(nevel_zonder_achtergrond) 
    m_nevel = im_nevel + delta_m 


    # Calculate the total area of the sub-region
    total_area = (x_end - x_start)*0.4310157882516833 * (y_end - y_start)*0.4310157882516833

    # Calculate the area of the mask region
    mask_area = (x_mask_end - x_mask_start)*0.4310157882516833 * (y_mask_end - y_mask_start)*0.4310157882516833

    # Calculate the area of the unmasked region
    unmasked_area = (total_area - mask_area)

    surface_brightness = m_nevel + 2.5 * np.log10(0.4310157882516833**2)

    # z = ZScaleInterval()
    # z1,z2 = z.get_limits(surface_brightness)
    plt.imshow(surface_brightness, cmap='gray')
    plt.colorbar()
    #plt.show()


    # Save the masked data to a new FITS file
    # masked_hdu = fits.PrimaryHDU(sub_data, header=header)
    # masked_hdu.writeto('masked_data.fits', overwrite=True)

    return surface_brightness, sub_data

#import fits file
sb_i, sub_data_i, im_nevel_i, err_mag_i = brightness_i("/Users/ilias/OneDrive/Documenten/GitHub/irisnevel/data/processed/ngc7023 i filter.fit")
# sb_r, sub_data_r = brightness_r("/Users/ilias/OneDrive/Documenten/GitHub/irisnevel/data/processed/ngc7023 r filter aligned.fit")
sb_g, sub_data_g, im_nevel_g, err_mag_g = brightness_g('/Users/ilias/OneDrive/Documenten/GitHub/irisnevel/data/processed/g60s aligned.fit')

print(err_mag_i)
print(err_mag_g)


sb_i_lijn = sb_i[214:403, 513]
sb_g_lijn = sb_g[214:403, 513]

err_mag = np.sqrt((err_mag_g)**2 + (err_mag_i)**2)
# err_mag = np.abs(sb_g_lijn/sb_i_lijn) * np.sqrt((err_mag_g / sb_g_lijn)**2 + (err_mag_i / sb_i_lijn)**2)
print(err_mag)

# print(sb_i_lijn)
# print(sb_g_lijn)

# Crop images to the same shape
min_shape = np.minimum(sub_data_i.shape, sub_data_g.shape)
sub_data_i_cropped = sub_data_i[:min_shape[0], :min_shape[1]]
sub_data_g_cropped = sub_data_g[:min_shape[0], :min_shape[1]]

# Align the images using cross-correlation
shift, error, diffphase = phase_cross_correlation(sub_data_i_cropped, sub_data_g_cropped, upsample_factor=100)

# Apply the shift to the g band image
shifted_g = np.fft.ifftn(fourier_shift(np.fft.fftn(sub_data_g_cropped), shift)).real

# Compute the ratio and plot
# verhouding_gi = shifted_g - sub_data_i_cropped

verhouding_gi = sb_g - sb_i
# verhouding_gr = sb_g - sb_r
verhouding_im_gi = im_nevel_g - im_nevel_i


figure1 = plt.figure('g-i')
z = ZScaleInterval()
z1,z2 = z.get_limits(verhouding_gi)
plt.imshow(verhouding_gi, vmin=-1, vmax=1, cmap='gray_r')
plt.colorbar()

hdu = fits.PrimaryHDU(verhouding_gi)
hdul = fits.HDUList([hdu])
hdul.writeto('verhouding_gi.fits', overwrite=True)

# figure2 = plt.figure('g-r')
# z = ZScaleInterval()
# z1,z2 = z.get_limits(verhouding_gr)
# plt.imshow(verhouding_gr, vmin=-1, vmax=1, cmap='gray_r')
# plt.colorbar()

# hdu = fits.PrimaryHDU(verhouding_gr)
# hdul = fits.HDUList([hdu])
# hdul.writeto('verhouding_gr.fits', overwrite=True)

#plt.imshow(verhouding_im_gi, vmin=z1, vmax=z2, cmap='gray')


# plt.show()

# Extract values from the specified line
x = 513
y_start, y_end = 214, 403

# Extract values along the specified vertical line
line_values = verhouding_gi[y_start:y_end, x]

# Create a new figure for the plot
plt.figure(figsize=(10, 6))

# # Plot the values along the vertical line
# plt.plot(range(y_start, y_end), line_values, color = 'w', label='g/i ratio in the nebula', linestyle = 'o')
plt.errorbar(range(y_start, y_end), line_values, yerr=err_mag, ecolor='lightgrey', fmt='o', color='w', label='g/i ratio in the nebula')
plt.axhline(y=-0.24, color='r', linestyle='-', label='g/i ratio in HD200775')
plt.xlabel('Y Pixel Coordinate', size = 20, color = 'w')
plt.ylabel("g'- i'", size = 20, color = 'w')
plt.ylim(-3,2)
plt.grid(True)

# # Adjust the ticks on the axes for the image
plt.xticks(fontsize=20, color='white')  # Adjust the font size for x-axis ticks
plt.yticks(fontsize=20, color='white')  # Adjust the font size for y-axis ticks

plt.gca().spines['bottom'].set_color('white')
plt.gca().spines['top'].set_color('white') 
plt.gca().spines['right'].set_color('white')
plt.gca().spines['left'].set_color('white')

plt.gca().tick_params(axis='x', colors='white')
plt.gca().tick_params(axis='y', colors='white')
plt.savefig('gi verhouding2', transparent=True, bbox_inches='tight')
plt.show()

# # z = ZScaleInterval()
# # z1,z2 = z.get_limits(sub_data_i)
# # plt.imshow(sub_data_i, vmin=z1, vmax=z2, cmap='gray')

# plt.imshow(verhouding_gi, vmin=-2, vmax=2, cmap='gray')
# cbar = plt.colorbar()

# # Set the colorbar ticks and label color to white
# cbar.ax.yaxis.set_tick_params(color='white', labelsize=15)
# cbar.ax.yaxis.set_tick_params(labelcolor='white')
# cbar.outline.set_edgecolor('white')

# # Mark the selected line on the image
# plt.plot([x, x], [y_start, y_end], color='red', linestyle='-', linewidth=1.5)
# plt.xlabel('X Pixel', size = 15, color='white')
# plt.ylabel('Y Pixel', size = 15, color='white')

# # Adjust the ticks on the axes for the image
# plt.xticks(fontsize=15, color='white')  # Adjust the font size for x-axis ticks
# plt.yticks(fontsize=15, color='white')  # Adjust the font size for y-axis ticks

# plt.gca().spines['bottom'].set_color('white')
# plt.gca().spines['top'].set_color('white') 
# plt.gca().spines['right'].set_color('white')
# plt.gca().spines['left'].set_color('white')

# plt.gca().tick_params(axis='x', colors='white')
# plt.gca().tick_params(axis='y', colors='white')
# plt.savefig('lijn op plaatje2', transparent=True, bbox_inches='tight')
# # plt.show()

# # print(verhouding_gi, verhouding_gr)

# # g naar r +2 -1/+1









