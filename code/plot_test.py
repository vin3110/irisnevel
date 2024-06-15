import matplotlib 
import matplotlib.pyplot as plt
import numpy as np

# plt.figure()
# plt.xlim(0, 5)
# plt.ylim(0, 5)

# for i in range(0, 5):
#     plt.axhspan(i, i+.2, facecolor='0.2', alpha=0.5)
#     plt.axvspan(i, i+.5, facecolor='b', alpha=0.5)




t = np.arange(0.0, 2, 0.01)
s = np.sin(2*np.pi*t)

fig, ax = plt.subplots()

ax.plot(t, s, color='black')
ax.axhline(0, color='black')

ax.fill_between([0,1], 2, facecolor='green', alpha=.5)
ax.fill_between(t, -1,  facecolor='red', alpha=.5)


# plt.xlim(0,10)
# plt.ylim(-10,0)

plt.show()