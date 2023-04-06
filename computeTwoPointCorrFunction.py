import sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import porespy as ps    # Calculates the two-point correlation function using Fourier transforms. see https://porespy.org/examples/metrics/reference/two_point_correlation.html
import inspect
import random
inspect.signature(ps.metrics.two_point_correlation)

#np.set_printoptions(threshold=sys.maxsize)

np.random.seed(10)
#im = ps.generators.blobs(shape=[100,100,100])
#print('im: ',im)
#fig, ax = plt.subplots(1, 1, figsize=[6, 6])
#ax.imshow(im[:,:,9], origin='lower', interpolation='none')
#ax.axis(False);
#plt.show()
#data = ps.metrics.two_point_correlation(im)
#print('data.distance:',data.distance)
#print('data.probability:',data.probability)
#fig, ax = plt.subplots(1, 1, figsize=[6, 6])
#ax.plot(data.distance, data.probability, 'r.')
#ax.set_xlabel("distance")
#ax.set_ylabel("two point correlation function");
#plt.show()



#data = np.array([[[1,1,0], [0,1,0], [0,1,1]], [[0,0,0], [0,1,0], [1,0,0]], [[0,0,1], [1, 1, 1], [0, 1, 0]]], np.int32)
#data = np.array([[[255,255,0], [0,255,0], [0,255,255]], [[0,0,0], [0,255,0], [255,0,0]], [[0,0,255], [255, 255, 255], [0, 255, 0]]], np.int32)
#data = np.random.randint(2, size=(80, 80, 80))
#print('data: ', data)
#fig, ax = plt.subplots(1, 1, figsize=[6, 6])
#ax.imshow(data[:,:,2], origin='lower', interpolation='none')
#ax.axis(False);
#plt.show()
#dataCorr = ps.metrics.two_point_correlation(data)
#print('data.distance:',dataCorr.distance)
#print('data.probability:',dataCorr.probability)
#fig, ax = plt.subplots(1, 1, figsize=[6, 6])
#ax.plot(dataCorr.distance, dataCorr.probability, 'r.')
#ax.set_xlabel("distance")
#ax.set_ylabel("two point correlation function");
#plt.show()



# Read the synthetic image and compute the 2-point corr function in 2D or for different slices
performSliceWiseComputation = 'false' # slices are taken along one of the directions
dataSynthetic = np.load('WC-720-51-final-imageOutput.npy')

fig, ax = plt.subplots(1, 1, figsize=[6, 6])
if performSliceWiseComputation == 'true':
    for iSlice in range(0,dataSynthetic.shape[0]):
        dataCorrSynthetic = ps.metrics.two_point_correlation(dataSynthetic[:,:,iSlice],voxel_size=1, bins=100)    
        ax.plot(dataCorrSynthetic.distance, dataCorrSynthetic.probability, 'r-')
else:        
    dataCorrSynthetic = ps.metrics.two_point_correlation((1/255)*dataSynthetic,voxel_size=1, bins=100)    
    ax.plot(dataCorrSynthetic.distance, dataCorrSynthetic.probability, 'r-', label='Synthetic Image - 3D')
    
# Read the reference image and compute the 2-point corr function
dataReferenceImg = mpimg.imread('WC-720-51-final.png') # shape of dataReferenceImg = (pixelsY,pixelsX,nChannels), nChannels=3 for RGB image and nChannels=4 for RGBA image
if dataReferenceImg.shape[2] >= 3:
    dataReferenceImg_blackWhite = (1/3)*dataReferenceImg[:,:,0] + (1/3)*dataReferenceImg[:,:,1] + (1/3)*dataReferenceImg[:,:,2]

dataCorrReal = ps.metrics.two_point_correlation(dataReferenceImg_blackWhite,voxel_size=1, bins=100)

ax.plot(dataCorrReal.distance, dataCorrReal.probability, 'r.', label='Reference Image - 2D')
print(len(dataCorrReal.distance))
print(dataCorrSynthetic.distance[len(dataCorrSynthetic.distance)-1])
ax.set_xlim(0, 1.1*np.amin([dataCorrReal.distance[len(dataCorrReal.distance)-1],dataCorrSynthetic.distance[len(dataCorrSynthetic.distance)-1]]))
#ax.set_xlim(0,50)
ax.set_xlabel("Distance (r)")
ax.set_ylabel("$S_2(r)$");
ax.legend()
plt.show()