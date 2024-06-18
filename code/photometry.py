from photutils.aperture import aperture_photometry, CircularAperture
from photutils.aperture import CircularAnnulus, ApertureStats
import numpy as np
# import glob
# import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib.colors import LogNorm
import math as mt

# from astropy.stats import sigma_clipped_stats
# from photutils.detection import DAOStarFinder
# from photutils.aperture import CircularAperture



data_in ="/Users/ilias/OneDrive/Documenten/GitHub/irisnevel/data/processed/g60s.fit"

hdu = fits.open(data_in)[0]
header = hdu.header 
data = hdu.data

z = ZScaleInterval()
z1,z2 = z.get_limits(data)




positions = [(2798, 1630), (2561, 3125)]
radius = 10
aperture = CircularAperture(positions, r = radius)

phot_table = aperture_photometry(data, aperture)

annulus_aperture = CircularAnnulus(positions, r_in = 13, r_out = 22)
aperstats = ApertureStats(data, annulus_aperture)
bkg_mean = aperstats.mean
total_bkg = bkg_mean * aperture.area


star_data = aperture_photometry(data, aperture)
star_data['total_bkg'] = total_bkg

figure = plt.figure()
aperture.plot(color='red', lw=1.5, alpha=0.5)
plt.imshow(data, cmap='Greys', vmin=z1, vmax=z2)
plt.show()


for col in star_data.colnames:
    star_data[col].info.format = '%8g'

star_data.pprint()

phot_bkgsub = phot_table['aperture_sum'] - total_bkg
# print(total_bkg)
# print(phot_table['aperture_sum'])


im = -2.5 * np.log10(phot_bkgsub)
print(im)