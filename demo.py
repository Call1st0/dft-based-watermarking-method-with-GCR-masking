#%%
from PIL import Image, ImageCms
import skimage.measure as msr
from matplotlib import pyplot as plt
import numpy as np
from wmark import WaterMark
#%%
# Create wmark object with seed 5
wmarkObj = WaterMark(5)

# Read and image
img = wmarkObj.imread('TestSet/cmykExample.tif')

# Find appropriate implementation factor to get quality (PSNR) between 20, 30 dB
factor = wmarkObj.findImpactFactor(img,(20,30))
print(factor)

# Embed mark in an image
img_marked = wmarkObj.embedMark(img, factor=factor)

# Mask wateremark using gcr masking
img_masked = wmarkObj.gcrMasking(img, img_marked, 'ISOcoated_v2_eci.icc')

# Show images
WaterMark.compareImages(img, img_marked, img_masked)

# Calculate and print psrn and ssim in cmyk 
ssim_marked_cmyk = msr.compare_ssim(img,img_marked,multichannel=True)
psnr_marked_cmyk = msr.compare_psnr(img,img_marked)
ssim_masked_cmyk = msr.compare_ssim(img,img_masked,multichannel=True)
psnr_masked_cmyk = msr.compare_psnr(img,img_masked)

print(f'Marked Image CMYK >> SSIM = {ssim_marked_cmyk:.2f}, PSNR = {psnr_marked_cmyk:.2f}')
print(f'Masked Image CMYK >> SSIM = {ssim_marked_cmyk:.2f}, PSNR = {psnr_marked_cmyk:.2f}')

# Convert to sRGB
img_srgb = WaterMark.profileCmyk2srgb(img, profileCmyk='profiles/ISOcoated_v2_eci.icc')
img_marked_srgb = WaterMark.profileCmyk2srgb(img_marked, profileCmyk='profiles/ISOcoated_v2_eci.icc')
img_masked_srgb = WaterMark.profileCmyk2srgb(img_masked, profileCmyk='profiles/ISOcoated_v2_eci.icc')

# Calculate and print psrn and ssim in srgb 
ssim_marked_srgb = msr.compare_ssim(img_srgb, img_marked_srgb,multichannel=True)
psnr_marked_srgb = msr.compare_psnr(img_srgb, img_marked_srgb)
ssim_masked_srgb = msr.compare_ssim(img_srgb, img_masked_srgb,multichannel=True)
psnr_masked_srgb = msr.compare_psnr(img_srgb, img_masked_srgb)

print(f'Marked Image sRGB >> SSIM = {ssim_marked_cmyk:.2f}, PSNR = {psnr_marked_srgb:.2f}')
print(f'Masked Image sRGB >> SSIM = {ssim_masked_cmyk:.2f}, PSNR = {psnr_masked_srgb:.2f}')

