#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2019 Ante Poljicak
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.


import skimage
import skimage.measure as msr
import skimage.color as color
import skimage.transform as trs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import imageio
from PIL import Image
from scipy import fftpack
import math
from enum import Enum


class WaterMark:

    # instance atributes

    # constructor with seed
    def __init__(self, seed):
        self.seed = seed
        pass

    @staticmethod
    def compareImages(originalImg, modifiedImg):
        """ Helper method for displaying comparison of images """
        fig, axes = plt.subplots(nrows=1, ncols=2, sharex='all', sharey='all')
        # ax = axes.ravel()

        psnr_orig = msr.compare_psnr(originalImg, originalImg)
        ssim_orig = msr.compare_ssim(originalImg, originalImg)

        psnr_mod = msr.compare_psnr(originalImg, modifiedImg)
        ssim_mod = msr.compare_ssim(originalImg, modifiedImg)

        label = 'PSNR: {:.2f}, SSIM: {:.2f}'

        axes[0].imshow(originalImg, cmap=plt.cm.gray)
        axes[0].set_xlabel(label.format(psnr_orig, ssim_orig))
        axes[0].set_title('Original image')

        axes[1].imshow(modifiedImg, cmap=plt.cm.gray)
        axes[1].set_xlabel(label.format(psnr_mod, ssim_mod))
        axes[1].set_title('Modified image')

        plt.show()

    def imread(self, imgName):
        """ read image and return np array """
        return np.array(Image.open(imgName))

    def inputProc(self, img):
        """ Takes image and returns its magnitude and phase in fourier domain"""
        # convert to float64 to keep precision
        skimage.img_as_float64(img)
        fft2 = fftpack.fft2(img)
        # shift lowest frequency to the center
        magnitude = fftpack.fftshift(np.absolute(fft2))
        phase = np.angle(fft2)
        return magnitude, phase

    def pseudoGen(self, length):
        """ Takes lenght and returns binary pseudorandom vector of that length """
        np.random.seed(self.seed)
        return np.random.randint(2, size=(length, 1)).astype(np.float32)

    def outputProc(self, magnitude, phase):
        """ Takes magnitude and phase of the image and returns img in spatial domain"""
        img = fftpack.ifft2(np.multiply(
            fftpack.ifftshift(magnitude), np.exp(1j * phase)))
        img = np.real(img)  # Image still is complex so take only real part
        return img

    @staticmethod
    def getMarkChannel(img):
        """getter for channel used for embeding watermark

        Arguments:
            img {[type]} -- [description]

        Raises:
            TypeError: [description]
            TypeError: [description]

        Returns:
            [type] -- [description]
        """
        # check if the image is grayscale, RGB or CMYK image
        if len(img.shape) == 2:
            image_type = 'GRAY'
            img_y = img
        elif len(img.shape) == 3:
            if img.shape[2] == 4:  # CMYK image
                image_type = 'CMYK'
                img_y = img[:, :, 3]
            elif img.shape[2] == 3:  # RBG image
                image_type = 'RGB'
                img_ycbcr = color.rgb2ycbcr(img)
                img_y = img_ycbcr[:, :, 0]
            else:
                raise TypeError(
                    "Image has wrong number of channels. Only 3 or 4 channels allowed")
        else:
            raise TypeError(
                "Image isn't of correct type. Only grayscale, RGB and CMYK allowed")

        return img_y, image_type

    def embedMark(self, img, length=200, frequencies='MEDIUM', factor=5):
        """ Method for embeding watermark in magnitude spectrum of the image

        Arguments:
            img {ndarray} -- host image

        Keyword Arguments:
            length {int} -- length of the watermark (default: {200})
            frequencies {str} -- frequencies in which watermark will be embedded (default: {'MEDIUM'})
            factor {int} -- implementation strength (default: {5})

        Returns:
            {ndarray} -- [description]
        """

        img_y, image_type = WaterMark.getMarkChannel(img)

        radius = WaterMark.frequency2Radius(img, frequencies)

        # Generate prnd sequence
        mark = self.pseudoGen(length)

        # Transformation in Fourier domain
        magnitude, phase = self.inputProc(img_y)

        # Embed the vector into the image
        watermark_mask = np.zeros(img_y.shape)

        # 2 masks are needed to ensure simetric modification of magnitude spectrum
        mask_shape = (3, 3)
        mask_up = np.zeros(mask_shape)
        mask_down = np.zeros(mask_shape)

        for ind in range(length):
            x1 = int((img_y.shape[0]/2) +
                     np.around(radius*math.cos(ind*math.pi/length)))
            y1 = int((img_y.shape[0]/2) +
                     np.around(radius*math.sin(ind*math.pi/length)))
            x2 = int((img_y.shape[0]/2) + np.around(radius *
                                                    math.cos(ind*math.pi/length+math.pi)))
            y2 = int((img_y.shape[0]/2) + np.around(radius *
                                                    math.sin(ind*math.pi/length+math.pi)))
            for ind_m_x in range(3):
                for ind_m_y in range(3):
                    mask_up[ind_m_x, ind_m_y] = img_y[(
                        x1-1+ind_m_x), (y1-1+ind_m_y)]
                    mask_down[ind_m_x, ind_m_y] = img_y[(
                        x2-1+ind_m_x), (y2-1+ind_m_y)]

            # build watermark mask by multiplying mark bits
            # with mean value of the surrounding coefficients
            watermark_mask[x1, y1] = mark[ind]*np.mean(mask_up)
            watermark_mask[x2, y2] = mark[ind]*np.mean(mask_down)

        # add watermark mask multiplied by factor to the original magnitude
        magnitude_m = magnitude+factor*watermark_mask

        # Transformation to Spatial domain
        img_y_marked = self.outputProc(magnitude_m, phase)

        if image_type == 'RGB':
            img_ycbcr = color.rgb2ycbcr(img)
            img_ycbcr_m = np.dstack((img_y_marked, img_ycbcr[:, :, 1:]))
            img_m = color.ycbcr2rgb(img_ycbcr_m)
        elif image_type == 'CMYK':
            img_m = np.dstack((img[:, :, 0:3], img_y_marked))
        else:
            img_m = img_y_marked

        # TODO this processing will change the data of unprocessed channels too. Figure out a way to fix that
        # check if image is in the -1 , 1 interval
        if np.amax(img_m) > 1:
            img_m = img_m/np.amax(img_m)

        # return image as uint8
        return skimage.img_as_ubyte(img_m)

    def decodeMark(self, img, length=200, frequencies='MEDIUM'):
        """Method for extracting ana decoding hidden watermark
        
        Arguments:
            img {ndarray} -- host image
        
        Keyword Arguments:
            length {int} -- length of embeded (watermark) 1D vector  (default: {200})
            frequencies {str} -- frequencies at which to search for vector ('LOW', 'MEDIUM', 'HIGH') (default: {'MEDIUM'})
        
        Returns:
            float -- correlation value of the extracted vector and generated watermark
        """
        img_y, image_type = WaterMark.getMarkChannel(img)
        magnitude, phase = self.inputProc(img_y)

        mark = self.pseudoGen(length)
        radius = WaterMark.frequency2Radius(img, frequencies)

        corr = np.zeros(32)
        counter = 0
        for ind in range(int(radius-16), int(radius+16)):
            vec = WaterMark.extractMark(magnitude, ind)
            reshapedMark = WaterMark.generateReshapedMark(mark, magnitude, ind)
            corr[counter] = WaterMark.covarMark(reshapedMark, vec)
            counter += 1

        return np.amax(corr)

    @staticmethod
    def extractMark(img, radius):
        """ Method for extracting vector from the ndarray given the radius

        Arguments:
            img {ndarray} -- host ndarray
            radius {int} -- radius at which to extract vector

        Returns:
            {ndarray} -- Extracted vector
        """
        # step between coefficients used for embedding.
        step = math.pi/(2*math.asin(1/(2*radius)))
        vec = np.zeros((math.ceil(step), 1))
        mask = np.zeros((3, 3))

        for ind in range(math.ceil(step)):
            x1 = int((img.shape[0]/2) +
                     np.around(radius*math.cos(ind*math.pi/step)))
            y1 = int((img.shape[0]/2) +
                     np.around(radius*math.sin(ind*math.pi/step)))

            for ind_m_x in range(3):
                for ind_m_y in range(3):
                    mask[ind_m_x, ind_m_y] = img[(
                        x1-1+ind_m_x), (y1-1+ind_m_y)]
            vec[ind, 0] = np.amax(mask)
        return vec

    def generateReshapedMark(mark, img, radius=128):
        """Generator for watermark mask
        
        Arguments:
            mark {ndarray} -- 1D PRND vector generated using pseudoGen method
            img {ndarray} -- host image
        
        Keyword Arguments:
            radius {int} -- radius at which to embed PRND vector (default: {128})
        
        Returns:
            [type] -- [description]
        """
        watermark_mask = np.zeros(img.shape)
        length = len(mark)

        for ind in range(length):
            x1 = int((img.shape[0]/2) +
                     np.around(radius*math.cos(ind*math.pi/length)))
            y1 = int((img.shape[0]/2) +
                     np.around(radius*math.sin(ind*math.pi/length)))

            watermark_mask[x1, y1] = mark[ind]

        return WaterMark.extractMark(watermark_mask, radius)

    def frequency2Radius(img, frequencies):
        """convert frequency to radius

        Arguments:
            img {ndarray} -- image
            frequency {str} -- 'LOW', 'MEDIUM' or 'HIGH'

        Raises:
            AttributeError: Raises if wrong frequency zone is passed

        Returns:
            float -- radius representing provided frequency zone
        """
        # define radius for wanted zone
        if frequencies == 'LOW':
            radius = 1/8*(img.shape[0])
        elif frequencies == 'MEDIUM':
            radius = 1/4*(img.shape[0])
        elif frequencies == 'HIGH':
            radius = 3/8*(img.shape[0])
        else:
            raise AttributeError(
                "Unknown frequency zone. Use 'LOW', 'MEDIUM', or 'HIGH' ")

        return radius

    @staticmethod
    def covarMark(mark, vector):
        """Cross correlation of two 1D sequences that are first made zero-mean

        Arguments:
            mark {ndarray} -- array like input sequences
            vector {ndarray} -- array like input sequences

        Returns:
            ndarray -- 1D vector of full cross correlation of two args
        """
        # normalize to zero mean
        mark = mark - np.mean(mark)
        vector = vector - np.mean(vector)
        return (np.cov(mark[:, 0], vector[:, 0])[0][1])  

    def isReplaceable(img, repCMYmin, repKmax=100):
        """Checks if pixels in an image can be GCR-d
        
        Arguments:
            im {ndarray} -- image
            repCMYmin {float} -- minimal value for CMY to be considered replacable
            repKmax {float} -- maximum value of K to be considered replacable 
        
        Returns:
            ndarray -- array of bools with the same shape as img inidicating replacebility
        """

        img = np.ndarray.astype(img, dtype=float)
        imsh0 = img.shape[0]
        imsh1 = img.shape[1]
        img = img.reshape(np.int((imsh0*imsh1)), img.shape[2])
        rpl = np.logical_and(np.logical_and(np.logical_and(img[:, 0]/2.55 > repCMYmin, img[:, 1]/2.55 > repCMYmin),
                                            img[:, 2]/2.55 > repCMYmin), img[:, 3]/2.55 < repKmax)
        return rpl.reshape((imsh0, imsh1))
