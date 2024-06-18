import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits


data_g ="/NGC7023/20240307/NGC7023_reduced_g'_7s.fit"
data_r ="/NGC7023/L60s.fit"
data_i ="/NGC7023/ngc7023 i filter.fit"



hdu = fits.open(data_g)[0]
header = hdu.header
data_g = hdu.data

hdu = fits.open(data_r)[0]
header = hdu.header
data_r = hdu.data


hdu = fits.open(data_i)[0]
header = hdu.header
data_i = hdu.data


data1 = data_g - data_i

figure1=plt.figure()
z = ZScaleInterval()
z1,z2 = z.get_limits(data1)
plt.imshow(data1, vmin=z1, vmax=z2, cmap='gray')

figure2=plt.figure()
data2 = data_g - data_r

z = ZScaleInterval()
z1,z2 = z.get_limits(data2)
plt.imshow(data2, vmin=z1, vmax=z2, cmap='gray')

plt.show()