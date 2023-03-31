import sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import porespy as ps
import inspect
import random
inspect.signature(ps.metrics.two_point_correlation)

#np.set_printoptions(threshold=sys.maxsize)

#np.random.seed(10)
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





dataSynthetic = (1/255)*np.load('WC-720-51-final-imageOutput.npy')
sizeSynthImage = dataSynthetic.shape[0]
print('dataSyntheticImage: ',dataSynthetic.shape) 
sliceNumber = random.randrange(sizeSynthImage)
fig, ax = plt.subplots(1, 1, figsize=[6, 6])
ax.imshow(dataSynthetic[:,:,sliceNumber], origin='lower', interpolation='none')
ax.axis(False);
plt.show()
dataCorrSynthetic = ps.metrics.two_point_correlation(dataSynthetic)
fig, ax = plt.subplots(1, 1, figsize=[6, 6])
ax.plot(dataCorrSynthetic.distance, dataCorrSynthetic.probability, 'r.')
ax.set_xlabel("distance")
ax.set_ylabel("two point correlation function");
plt.show()

## Read Images
#dataReferenceImg = mpimg.imread('WC-720-51-final.png')
#print('dataReferenceImg: ',dataReferenceImg.shape) 
## Output Images
#plt.imshow(dataReferenceImg)
#plt.show()
#dataCorrReal = ps.metrics.two_point_correlation(dataReferenceImg)
#fig, ax = plt.subplots(1, 1, figsize=[6, 6])
#ax.plot(dataCorrReal.distance, dataCorrReal.probability, 'r.')
#ax.set_xlabel("distance")
#ax.set_ylabel("two point correlation function");
#plt.show()