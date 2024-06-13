import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib.colors import LogNorm

from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture

data_in ='/NGC7023/L60s.fit'

hdu = fits.open(data_in)[0]
header = hdu.header 
data = hdu.data

z = ZScaleInterval()
z1,z2 = z.get_limits(data)

# plt.imshow(data, cmap='Greys', vmin=z1, vmax=z2)

def positions(section):

    mean, median, std = sigma_clipped_stats(section, sigma=3.0)
    
    daofind = DAOStarFinder(fwhm=4.5, threshold=5.*std)

    sources = daofind(section - median)

    # for col in sources.colnames:
    #     if col not in ('id', 'npix'):
    #         sources[col].info.format = '%.2f'
    # sources.pprint(max_width=76)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    return positions


def sectioning_and_starfind(data, sqrt_of_parts):
    dimention_x = data.shape[0]
    dimention_y = data.shape[1]
    x = int(dimention_x/sqrt_of_parts)
    y = int(dimention_y/sqrt_of_parts)
    positions_list = []

    for i in range(sqrt_of_parts):
        for j in range(sqrt_of_parts):
            section = data[i*x:(i+1)*x, j*y:(j+1)*y]
            # print(i*x,(i+1)*x, j*y,(j+1)*y)
            positions_list.append(positions(section))

            ##---- Uncomment this section to see the sections and the stars found in them ----##
            # figure = plt.figure()
            # apatures = CircularAperture(positions_list[-1], r=4.0)
            # apatures.plot(color='red', lw=1.5, alpha=0.5)
            # plt.imshow(section, cmap='Greys', vmin=z1, vmax=z2)

    return positions_list

positions_list = sectioning_and_starfind(data, 4)

# for i in range(len(positions_list)):
#     apatures = CircularAperture(positions_list[i], r=4.0)
#     apatures.plot(color='red', lw=1.5, alpha=0.5)

# apatures = CircularAperture(positions_list[5], r=4.0)
# apatures.plot(color='red', lw=1.5, alpha=0.5)

plt.show()