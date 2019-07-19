import skimage.measure as msr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as style
import pandas as pd
from PIL import Image, ImageCms
import time

import glob
import os
from pathlib import Path
import imageio
import cv2
from wmark import WaterMark
import multiprocessing

# Source Directory
srcFolder = 'TestSetFull/'
# Source Path
srcPth = Path(srcFolder).resolve()

def procImg(img):
    # Read original image
    imgOriginal = WaterMark.imread(img)
    print(img)

    randomNum = 5
    counter = 0
    wObject = WaterMark(randomNum)
    
    # Find impact factor within a PSNR range
    ImpactFactor = wObject.findImpactFactor(imgOriginal, rangePSNR = (35,40))
    
    # Mark image
    imgMarked = wObject.embedMark(imgOriginal, factor = ImpactFactor)
        
    # GCR Masking
    imgGCR = WaterMark.gcrMasking(imgOriginal, imgMarked)
    
    #CROPPING
    #CROPPING
    #CROPPING
    
    # Width
    lx = imgGCR.shape[0]
    # Height
    ly = imgGCR.shape[1]
    
    # Generating cropped image
    # 1/8 from each side (result is by 1/4 smaller than original)
    imgCropped = imgGCR[lx // 8: - lx // 8, ly // 8: - ly // 8]
    
    #SCALING
    #SCALING
    #SCALING
    
    # Wanted value of scaled image (in pixels)
    height = 768
    width = 768
    
    # Generating scaled image
    imgScaled = cv2.resize(imgGCR, (width, height))

    #ROTATION
    #ROTATION
    #ROTATION
    
    # Rotation angle
    angle = 90
    
    # Defines center of image
    rows, cols = imgGCR.shape[:2]
    M = cv2.getRotationMatrix2D((cols / 2, rows / 2), angle, 1)
    
    # Generating rotated image
    imgRotated = cv2.warpAffine(imgGCR, M, (cols, rows))    

    # AFFINE TRANSFORMATION
    # AFFINE TRANSFORMATION
    # AFFINE TRANSFORMATION
    
    # Defining points in original image
    point1a = [50, 50]
    point1b = [200, 50]
    point1c = [50, 200]
    # Defining points in transformed image
    point2a = [10, 100]
    point2b = [200, 50]
    point2c = [100, 250]
    
    rows, cols, ch = imgGCR.shape
    # Compiling original points into single array
    pts1 = np.float32([point1a, point1b, point1c])
    # Compiling transformed points into single array
    pts2 = np.float32([point2a, point2b, point2c])
    M = cv2.getAffineTransform(pts1, pts2)
    
    #Generating affined image
    imgAffined = cv2.warpAffine(imgGCR, M, (cols, rows))
        
    # 2D CONVOLUTION
    # 2D CONVOLUTION
    # 2D CONVOLUTION
    
    # Horizontal deviation (intensity)
    kernelX = 3
    # Vertical deviation (intensity)
    kernelY = 3
    
    kernelMultiplication = kernelX * kernelY
    kernel = np.ones((kernelX, kernelY), np.float32)/kernelMultiplication
    
    #Generating convoluted image
    imgConvoluted = cv2.filter2D(imgGCR, -1, kernel)  
    
    # GAUSSIAN BLUR
    # GAUSSIAN BLUR
    # GAUSSIAN BLUR
    
    # Horizontal deviation (blur intensity)
    kernelX = 3
    # Vertical deviation (blur intensity)
    kernelY = 3
    
    # Generating blurred image
    imgBlurred = cv2.GaussianBlur(imgGCR, (kernelX, kernelY), cv2.BORDER_DEFAULT)
    
    # GAUSSIAN NOISE
    # GAUSSIAN NOISE
    # GAUSSIAN NOISE
    
    # Wanted value of rotation angle (in degrees)
    noiseAmount = 0.5
    
    # Generate Gaussian noise
    gauss = np.random.normal(0, noiseAmount, imgGCR.size)
    gauss = gauss.reshape(
        imgGCR.shape[0], imgGCR.shape[1], imgGCR.shape[2]).astype('uint8')
    
    # Generating Gaussian noised image
    imgNoised = cv2.add(imgGCR, gauss)
    
    metricOriginal = wObject.decodeMark(imgOriginal, 'CORR')
    results[counter, 0] = metricOriginal
    
    metricCropped = wObject.decodeMark(imgCropped, 'CORR')
    results[counter, 1] = metricCropped
    
    metricScaled = wObject.decodeMark(imgScaled, 'CORR')
    results[counter, 2] = metricScaled
    
    metricRotated = wObject.decodeMark(imgRotated, 'CORR')
    results[counter, 3] = metricRotated
    
    metricAffined = wObject.decodeMark(imgAffined, 'CORR')
    results[counter, 4] = metricAffined
    
    metricConvoluted = wObject.decodeMark(imgConvoluted, 'CORR')
    results[counter, 5] = metricConvoluted
    
    metricBlurred = wObject.decodeMark(imgBlurred, 'CORR')
    results[counter, 6] = metricBlurred
    
    metricNoised = wObject.decodeMark(imgNoised, 'CORR')
    results[counter, 7] = metricNoised
    
    metricGCR = wObject.decodeMark(imgGCR, 'CORR')
    results[counter, 8] = metricGCR
    
    # counter += 1
    

    return results

metricDataframe = pd.DataFrame()

# for i in range(1, 2):

    
# Generating random seed as integer up to 10000
# randomNums = np.random.randint(10000, size=(1))
# randomNum = randomNums[0]
# print(randomNum)

# wObject = WaterMark(randomNum)

# All TIF files in the src_path are now imgs
listOfImgs = glob.glob(srcFolder+'*.tif')

# Counting the number of images in src_folder
# listofimgs = os.listdir(srcPth)
# -1 because there os one hidden file in directory
numImages = len(listOfImgs)

# Creating zeros Numpy Array for results
results = np.zeros([ 1,9])

pool = multiprocessing.Pool() # no parameter provided in constructor, use all threads
results = pool.map(procImg, listOfImgs) # return data as list of arrays
pool.close()
pool.join()
        
metricValues = ["Original", "Cropped", "Scaled", "Rotated", "Affined", "2D Convoluted", "Blurred", "Noised", "GCR / Marked"]
# randomNum = 5
# metricSeeds = [f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}", f"Marked {randomNum}"]

# stack all rows fo the array to create dataframe, pass titles for columns
metricDataframe = pd.DataFrame(np.row_stack(results),columns=metricValues)   
print(metricDataframe.tail())
metricDataframe.to_pickle('metric_data_frame_GCR.pkl')
metricDataframe.to_csv('metric_data_frame_GCR.csv')
# metricDataframe = pd.concat([metricDataframe, finalMetric], axis = 1)