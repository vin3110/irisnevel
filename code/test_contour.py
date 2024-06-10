import glob
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits

from matplotlib.colors import LogNorm

# import fits file
# data_in ='/NGC7023/20240308/Output/iris_nevel_reduced_Ha_60s_minSNRof50.fit'
data_in ='/NGC7023/L60s_GausBlur2.fit'
hdu = fits.open(data_in)[0]
header = hdu.header 
data = hdu.data

# print(header)
# min=int(data.min())
# max=int(data.max())

z = ZScaleInterval()
z1,z2 = z.get_limits(data)


data_normalized = (data - z1)/(z2 - z1) # [(i-min)/(max-min) for i in data]


# data_normalized = (data - np.min(data))/(np.max(data) -np.min(data))

# z = ZScaleInterval()
# z1,z2 = z.get_limits(data_normalized)

section = data_normalized[1000:3200,1000:3200]
plt.imshow(section, vmin=0, vmax=1, cmap='Greys')

# levels = []
# for i in np.arange(0, 1, 0.1):
#     levels.append(i)

levels = np.linspace(0.45, 0.6, num=10)

plt.contourf(section, levels)

plt.title('NGC7023 L-band')

# Plot de afbeelding

plt.colorbar()
plt.show()

# print(levels)
# print(np.median(data_normalized),np.mean(data_normalized))
# print(data)
# print(data_normalized)
# print(data_normalized.shape)
# print(z1,z2)
