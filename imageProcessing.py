### https://www.geeksforgeeks.org/opencv-python-tutorial/
import cv2
import numpy as np  
  
# Read the image.
img = cv2.imread('TC-as-fabricated-cropped.png')
  
  
# Apply bilateral filter with d = 15, 
# sigmaColor = sigmaSpace = 75.
bilateral = cv2.bilateralFilter(img, 15, 75, 75)  
# Save the output.
cv2.imwrite('img_bilateral.png', bilateral)


# Gray-scaling of an image (https://www.geeksforgeeks.org/python-grayscaling-of-images-using-opencv/)
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
cv2.imwrite('img_grayScale.png', img1)

# 
ret, thresh1 = cv2.threshold(img1, 120, 255, cv2.THRESH_BINARY)
cv2.imwrite('img_BinaryThreshold.png', thresh1)

# 
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
cv2.imwrite('img_segmentThreshold.png', thresh)

# edge-detection (https://www.geeksforgeeks.org/image-processing-in-python-scaling-rotating-shifting-and-edge-detection/)
edges = cv2.Canny(img, 100, 200)
cv2.imwrite('img_edges.png', edges)

# adaptive thresholding (https://www.geeksforgeeks.org/python-thresholding-techniques-using-opencv-set-2-adaptive-thresholding/)
thresh1 = cv2.adaptiveThreshold(img1, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 199, 5)  
thresh2 = cv2.adaptiveThreshold(img1, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 199, 5)
cv2.imwrite('img_AdaptiveMean.png', thresh1)
cv2.imwrite('img_AdaptiveGaussian.png', thresh2)

# Erosion (https://www.geeksforgeeks.org/python-opencv-cv2-erode-method/)
kernel = np.ones((5, 5), np.uint8)
# Using cv2.erode() method 
image = cv2.erode(img, kernel)   
# Displaying the image 
cv2.imwrite('img_erode.png', image) 

# Dilation (https://www.geeksforgeeks.org/erosion-dilation-images-using-opencv-python/)
img_dilation = cv2.dilate(img, kernel, iterations=1)
cv2.imwrite('img_dilate.png', img_dilation)   

# Blurring (https://www.geeksforgeeks.org/python-image-blurring-using-opencv/)
# Gaussian Blur
Gaussian = cv2.GaussianBlur(img, (7, 7), 0)
cv2.imwrite('img_GaussianBlurring.png', Gaussian)
# Median Blur
median = cv2.medianBlur(img, 5)
cv2.imwrite('img_MedianBlurring.png', median)
# Bilateral Blur
bilateral = cv2.bilateralFilter(img, 9, 75, 75)
cv2.imwrite('img_BilateralBlurring.png', bilateral)


