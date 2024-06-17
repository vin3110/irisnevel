import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LogNorm

#TODO Zoek extent uit
#TODO 2 color maps
#TODO gaussian blur, if yes welke sigma etc
#TODO Data normaliseren ja of nee

# import fits file
data_in ='/NGC7023/L60s.fit'
data_in2 ='/NGC7023/L60s_GausBlur2.fit'

def run(data1,sigma,naam):
    hdu = fits.open(data1)[0]
    header = hdu.header 
    data1 = hdu.data
    # data1 = data1[:2500,:]
    # print(data.shape)

    figure1 = plt.figure()

    # Apply Gaussian blur
    sigma = sigma  # Standard deviation for Gaussian kernel. Adjust this value to your needs.
    data_blurred = gaussian_filter(data1, sigma)

    z = ZScaleInterval()
    z1,z2 = z.get_limits(data_blurred)

    data_normalized = (data_blurred - z1)/(z2 - z1)
    # data_normalized = data_blurred

    # data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))

    # z = ZScaleInterval()
    # z1,z2 = z.get_limits(data_normalized)

    quantiles = np.quantile(data_normalized,[0.1,0.95])
    plt.imshow(data_normalized, cmap='Greys', vmin=quantiles[0], vmax=quantiles[1], norm='log')


    # levels = []
    # for i in np.arange(0, 1, 0.1):
    #     levels.append(i)

    s1 = 0.44
    s2 = 0.62
    num_levels = 10

    levels = np.geomspace(s1, s2, num=num_levels)

    levels_center = np.geomspace(s2, z2, num=num_levels)

    plt.contour(data_normalized, levels=levels)
    plt.contour(data_normalized, levels=levels_center)
    plt.title(naam)

    # Plot de afbeelding
    
    plt.colorbar()
    

    # print(levels)
    # print(np.median(data_normalized),np.mean(data_normalized))
    # print(data)
    # print(data_normalized)
    # print(data_normalized.shape)
    # print(z1,z2)

    # print(np.max(data_normalized), np.min(data_normalized))

run(data_in, 0, "blur0")
run(data_in, 2, "blur2")
plt.show()