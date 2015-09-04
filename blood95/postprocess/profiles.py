import matplotlib.pyplot as plt
import numpy as np
import glob
import convex_hull

rvol = '95'
rmax = '90'

fname = './profiles_bem/shape_conf' + rmax + '_v1eb10.dat'
print 'Reading ' + fname + ' ...'

# load xyz coordinates
data = np.loadtxt(fname, unpack=True)

hull = convex_hull.convex_hull(data)
x = hull[:,0]
y = hull[:,1]

plt.plot(x, y, '.')
plt.show()

