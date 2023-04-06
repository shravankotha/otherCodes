import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

a = np.random.randint(0,255,[200,200,200])  # random integer matrix of 200X200X200 size containing integers between 0 and 255(excluding)
print(a)
fig = plt.figure()
ax = plt.subplot(projection='3d')
cmap = plt.get_cmap("viridis")
#print('cmap.colors:',cmap.colors)
norm= plt.Normalize(a.min(), a.max())   # normalizing the a matrix with (max-min). This is done element-wise
ax.voxels(np.ones_like(a), facecolors=cmap(norm(a)), edgecolor="black") # norm(a) has the same size as 'a'. cmap(norm(a)) has a size of 200X200X200X4. i.e. each voxel gets a color value represented by 4 numbers. cmap(norm(a)) performs the mapping of colors from the space of cmap.colors in the cmaps object above to norm(a) space
print('cmap(norm(a)):',cmap(norm(a)).shape)
#plt.savefig('dummy_prediction_voxels.png')
plt.show()