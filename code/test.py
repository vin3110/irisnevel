import numpy as np

levels = np.linspace(0.44, 0.62, num=10)

z1 = 0.44
z2 = 0.62
num_levels = 10

levels = np.geomspace(z1, z2, num=num_levels)

print(levels)
print(np.diff(levels))

