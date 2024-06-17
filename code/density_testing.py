import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib.colors import LogNorm

from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture

# data_in ="/processed/L60s.fit"
data_in ="/NGC7023/L60s.fit"


hdu = fits.open(data_in)[0]
header = hdu.header 
data = hdu.data

data = data[:4046,:4046]

z = ZScaleInterval()
z1,z2 = z.get_limits(data)

fig, ax = plt.subplots()

fig = plt.imshow(data, cmap='Greys', vmin=z1, vmax=z2)

def positions(section):

    mean, median, std = sigma_clipped_stats(section, sigma=3.0)
    
    daofind = DAOStarFinder(fwhm=4.5, threshold=5.*std)

    sources = daofind(section - median)

    if sources is not None:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    else:
        positions = []
    
    return positions


def sectioning_and_starfind(data, sqrt_of_parts):
    dimention_x = data.shape[0]
    dimention_y = data.shape[1]
    x = int(dimention_x/sqrt_of_parts)
    y = int(dimention_y/sqrt_of_parts)
    positions_list = []
    densities = []
    

    for i in range(sqrt_of_parts):
        for j in range(sqrt_of_parts):
            x_begin = i*x
            x_end = (i+1)*x
            y_begin = j*y
            y_end = (j+1)*y
            section = data[x_begin:x_end, y_begin:y_end]
            # print(i*x,(i+1)*x, j*y,(j+1)*y)
            position = positions(section)



            densities.append(density(x, y, position, x_begin, x_end, y_begin, y_end))


            if position is not None:
                for pos in position:
                    pos[0] += x_begin
                    pos[1] += y_begin

            positions_list.append(position)



            ##---- Uncomment this section to see the sections and the stars found in them ----##
            # figure = plt.figure()
            # apatures = CircularAperture(positions_list[-1], r=4.0)
            # apatures.plot(color='red', lw=1.5, alpha=0.5)
            # plt.imshow(section, cmap='Greys', vmin=z1, vmax=z2)

    return positions_list, densities



def density(xlen, ylen, position, x_begin, x_end, y_begin, y_end):
    area = xlen * ylen

    # gives density per arcsec^2
    density = (len(position) / area) * 0.4310157882516833

    threshold = 1.3574445334205194e-05

    if density > threshold:
        ax.fill_between([x_begin, x_end], y_begin, facecolor='green', alpha=.2)
    else: 
        ax.fill_between([x_begin, x_end], y_begin, facecolor='red', alpha=.2)

    # if density > threshold:
    #     plt.axhspan(x_begin,x_end, facecolor='g', alpha=0.2)

    #     ##--onderstaande lijn is miss dubbel--##
    #     plt.axvspan(y_begin, y_end, facecolor='g', alpha=0.2)
    # else:
    #     plt.axhspan(x_begin,x_end, facecolor='r', alpha=0.2)

    #     ##--onderstaande lijn is miss dubbel--##
    #     plt.axvspan(y_begin, y_end, facecolor='r', alpha=0.2)


    return density

positions_list, densities = sectioning_and_starfind(data, 16)
# print(densities)
# print(len(densities))
# print(densities)

for positions in positions_list:
    for position in positions:
        apatures = CircularAperture(position, r=4.0)
        apatures.plot(color='red', lw=1.5, alpha=0.5)

# plt.plot(positions_list)
# for i in range(len(positions_list)):
#     apatures = CircularAperture(positions_list[i], r=4.0)
#     apatures.plot(color='red', lw=1.5, alpha=0.5)

# apatures = CircularAperture(positions_list[5], r=4.0)
# apatures.plot(color='red', lw=1.5, alpha=0.5)

plt.show()


