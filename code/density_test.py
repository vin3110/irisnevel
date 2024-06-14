import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib.colors import LogNorm

from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture

data_in ="ngc7023 i filter.fit"

hdu = fits.open(data_in)[0]
header = hdu.header 
data = hdu.data

data = data[:4046,:4046]

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

    # if sources is not None:
        # Probeer de xcentroid te krijgen
        # print(len(sources['xcentroid']))
    # else:
        # Als de xcentroid niet bestaat, print een bericht en return None
        # print('no stars found')

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
    
    map = np.empty((64,64))

    for i in range(sqrt_of_parts):
        for j in range(sqrt_of_parts):
            map[i,j] = len(positions_list)
            section = data[i*x:(i+1)*x, j*y:(j+1)*y]
            # print(i*x,(i+1)*x, j*y,(j+1)*y)
            position = positions(section)
            positions_list.append(position)
            densities.append(density(x, y, position))

            # ##---- Uncomment this section to see the sections and the stars found in them ----##
            # figure = plt.figure()
            # apatures = CircularAperture(positions_list[-1], r=4.0)
            # apatures.plot(color='red', lw=1.5, alpha=0.5)
            # plt.imshow(section, cmap='Greys', vmin=z1, vmax=z2)

    return positions_list, densities




def density(xlen, ylen, position):
    area = xlen * ylen
    density = (len(position) / area) * 0.4310157882516833

    return density

threshold = 1.5844615292067567e-05



positions_list, densities = sectioning_and_starfind(data, 16)
#print(densities)

sqrt_of_parts = 16
dimension_x = data.shape[0]
dimension_y = data.shape[1]
x = int(dimension_x / sqrt_of_parts)
y = int(dimension_y / sqrt_of_parts)

fig, ax = plt.subplots(sqrt_of_parts, sqrt_of_parts, figsize=(16, 16))

# for i in range(sqrt_of_parts):
#     for j in range(sqrt_of_parts):
#         section = data[i * x:(i + 1) * x, j * y:(j + 1) * y]
#         ax[i, j].imshow(section, cmap='Greys', vmin=z1, vmax=z2)
        
#         color = 'green' if densities[i * sqrt_of_parts + j] > threshold else 'red'
#         apatures = CircularAperture(positions_list[i * sqrt_of_parts + j], r=4.0)
#         apatures.plot(ax=ax[i, j], color=color, lw=1.5, alpha=0.5)
        
#         ax[i, j].axis('off')

plt.tight_layout()
plt.show()