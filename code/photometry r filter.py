from photutils.aperture import aperture_photometry, CircularAperture
from photutils.aperture import CircularAnnulus, ApertureStats
import numpy as np
import glob
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib.colors import LogNorm
import math as mt

# from astropy.stats import sigma_clipped_stats
# from photutils.detection import DAOStarFinder
# from photutils.aperture import CircularAperture



# data_in ="/Users/vinve/OneDrive/Documenten/GitHub/irisnevel/data/processed/ngc7023 r filter goed.fit"

data_g ="/NGC7023/g60s aligned.fit"
data_r ="/NGC7023/ngc7023 r filter aligned.fit"
data_i ="/NGC7023/ngc7023 i filter.fit"


hdu = fits.open(data_r)[0]
header = hdu.header 
data_r = hdu.data

hdu = fits.open(data_g)[0]
header = hdu.header 
data_g = hdu.data


hdu = fits.open(data_i)[0]
header = hdu.header 
data_i = hdu.data

# data = -2.5 * np.log10(data)

def run(data):

    def subregion(data):
        x_start, x_end = 1692, 2550
        y_start, y_end = 1774, 2310
        sub_data = data[y_start:y_end, x_start:x_end]

        return sub_data

    # Coordinaten van de ijk sterren: IRAS...6758 (2546, 3185) - UCA...4810 (2783, 1689)


    def dM(positions, radius, data, bkgrnd):

        aperture = CircularAperture(positions, r = radius)

        phot_table = aperture_photometry(data, aperture)

        annulus_aperture = CircularAnnulus(positions, r_in = 13, r_out = 22)
        aperstats = ApertureStats(data, annulus_aperture)
        if bkgrnd == 0:
            bkg_mean = aperstats.mean
        else:
            bkg_mean = bkgrnd
        total_bkg = bkg_mean * aperture.area


        star_data = aperture_photometry(data, aperture)
        star_data['total_bkg'] = total_bkg

        for col in star_data.colnames:
            star_data[col].info.format = '%8g'

        star_data.pprint()

        phot_bkgsub = phot_table['aperture_sum'] - total_bkg
        # print(total_bkg)
        # print(phot_table['aperture_sum'])

        px_in_aperture = aperture.area

        im = -2.5 * np.log10(phot_bkgsub)
        print(im)

        if np.all(data == data_r):
            V = 15.88 # IRAS
            B = 16.77 # IRAS

            am_UCA = 13.700
            am_IRAS = V - 0.46*(B - V) + 0.11 

        elif np.all(data == data_i):
            am_UCA = 12.768
            am_IRAS = 15.02

        elif np.all(data == data_g):
            V = 15.88 # IRAS
            B = 16.77 # IRAS
            
            am_UCA =  15.128
            am_IRAS = V + 0.64*(B-V) - 0.13

        dM_UCA = am_UCA - im[0]
        dM_IRAS = am_IRAS - im[1]

        print(dM_UCA, dM_IRAS)
        dM = (dM_UCA + dM_IRAS) / 2
        return dM, bkg_mean, px_in_aperture, im


    positions_ijk = [(2783, 1689), (2546, 3185)]
    radius_ijk = 10

    dM, bkg_mean, px_in_aperture, im = dM(positions_ijk, radius_ijk, data, bkgrnd=0)

    bkg_mean = (bkg_mean[0] + bkg_mean[1]) / 2
    print(dM, bkg_mean)

    data = data - bkg_mean
    data = -2.5 * np.log10(data) + dM

    sub_data = subregion(data)

    # 320 390 390 320
    # Define the mask region within the sub-region
    x_mask_start, x_mask_end = 320 , 390 #2013+17 - x_start, 2082-58 - x_start
    y_mask_start, y_mask_end = 320 , 390 #2100+17 - y_start, 2158-58 - y_start

    # # Apply the mask
    sub_data[y_mask_start:y_mask_end, x_mask_start:x_mask_end] = None

    return sub_data


sub_data_g = run(data_g)
sub_data_r = run(data_i)

sub_data = sub_data_g - sub_data_r

z = ZScaleInterval()
z1,z2 = z.get_limits(sub_data)
figure = plt.figure()
# aperture.plot(color='red', lw=1.5, alpha=0.5)

plt.imshow(sub_data, cmap='viridis', vmin=z1, vmax=z2)
plt.colorbar()
plt.show()





# position_filament = (2215, 2025)
# radius_filament = 25

# dM_neh, bkg_mean, px_in_aperture_, im = dM(position_filament, radius_filament, data_r, bkgrnd=bkg_mean)

# print()

# print(dM, bkg_mean, px_in_aperture)






