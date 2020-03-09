#%%
from PIL import Image, ImageCms
import skimage.metrics as msr
from matplotlib import pyplot as plt
import numpy as np
from wmark import WaterMark
#%%
# Create wmark object with seed 5
wmarkObj = WaterMark(5)

# Read and image
img = wmarkObj.imread('TestSet/cmykExample.tif')

# Find appropriate implementation factor to get quality (PSNR) between 20, 30 dB
factorQuality = wmarkObj.findImpactFactor(img,(20,30))
print('for impact factor: {:.0f}, PSNR quality is: {:.2f}'.format(factorQuality[0],factorQuality[1]))

# Embed mark in an image
img_marked = wmarkObj.embedMark(img, factor=factorQuality[0])

# Mask wateremark using gcr masking
img_masked = wmarkObj.gcrMasking(img, img_marked, 'ISOcoated_v2_eci.icc')

# Show images
WaterMark.compareImages(img, img_marked, img_masked)

# Calculate and print psrn and ssim in cmyk 
psnr_marked_cmyk = msr.peak_signal_noise_ratio(img,img_marked)
psnr_masked_cmyk = msr.peak_signal_noise_ratio(img,img_masked)

print(f'Marked Image CMYK >>  PSNR = {psnr_marked_cmyk:.2f}')
print(f'Masked Image CMYK >>  PSNR = {psnr_marked_cmyk:.2f}')


# Calculate and print psrn in lab
psnr_marked_lab = wmarkObj.labPSNR(img, img_marked, profileName='ISOcoated_v2_eci.icc')
psnr_masked_lab = wmarkObj.labPSNR(img, img_masked, profileName='ISOcoated_v2_eci.icc')

print(f'Marked Image lab >> PSNR = {psnr_marked_lab:.2f}')
print(f'Masked Image lab >> PSNR = {psnr_masked_lab:.2f}')

# Determine if the image has watermark
maxCorr=wmarkObj.decodeMark(img_masked,'CORR')
isMarked=wmarkObj.detectOutlier(img_masked,'CORR',alpha=0.0001)
print('Watermark is embeded: {}\nMax correlation value is {:.3f}'.format(isMarked,maxCorr))