# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:42:49 2018

@author: Davor
"""

import numpy as np
import scipy as sp
from PIL import Image
import time
import timeit
import ctypes
from ctypes import cdll
# from ctypes import windll
import matplotlib as mpl
import matplotlib.pyplot as plt
from sys import platform

# Define auxiliary functions
PRINTPROGRESS = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_int)
def PrintProgress(current, total):
    print(str(current + 1) + ' / ' + str(total))


# Load correct dynamic library depending on the system
if platform == "linux" or platform == "linux2":
    dynamicLibrary = 'libGcr.so'
elif platform == "darwin":
    dynamicLibrary = 'libGcr.dylib'
elif platform == "win32":
    dynamicLibrary = 'libGcr.dll' #TODO fix this to make usuable on win platform

# Setup LittleCMS wrapper
<<<<<<< HEAD
lcmswrap = cdll.LoadLibrary('lib/libGcr.dylib')
||||||| merged common ancestors
lcmswrap = cdll.LoadLibrary('../lib/libGcr.dylib')
=======
lcmswrap = cdll.LoadLibrary('lib/'+ dynamicLibrary)
>>>>>>> a8907e89f4c657cafb0e1859cc10568a2903c263
lcmswrap.makecform.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lcmswrap.makecform.restype = ctypes.c_void_p
lcmswrap.applycform.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lcmswrap.applycform.restype = None
lcmswrap.deletecform.argtypes = [ctypes.c_void_p]
lcmswrap.deletecform.restype = None
lcmswrap.makeabscform.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_uint]
lcmswrap.makeabscform.restype = ctypes.c_void_p
lcmswrap.applyabscform.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
lcmswrap.applyabscform.restype = None
lcmswrap.deleteabscform.argtypes = [ctypes.c_void_p]
lcmswrap.deleteabscform.restype = None

#Setup RevInt wrapper
<<<<<<< HEAD
revint = cdll.LoadLibrary('lib/libGcr.dylib')
||||||| merged common ancestors
revint = cdll.LoadLibrary('../lib/libGcr.dylib')
=======
revint = cdll.LoadLibrary('lib/'+ dynamicLibrary)
>>>>>>> a8907e89f4c657cafb0e1859cc10568a2903c263
revint.makeRevInt.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
revint.makeRevInt.restype = ctypes.c_void_p
revint.deleteRevInt.argtypes=[ctypes.c_void_p]
revint.deleteRevInt.restype = None
revint.makegammapForm.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_uint]
revint.makegammapForm.restype = ctypes.c_void_p
revint.makelab2kForm.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_uint]
revint.makelabk2cmykForm.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t]
revint.makelabk2cmykForm.restype = ctypes.c_void_p
revint.applygammapForm.argtypes= [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_size_t]
revint.applygammapForm.restype = None
revint.applylab2kForm.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_uint]
revint.applylab2kForm.restype = None
revint.applylabk2cmykForm.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
revint.applylabk2cmykForm.restype = None
revint.deletegammapForm.argtypes = [ctypes.c_void_p]
revint.deletegammapForm.restype = None
revint.deletelab2kForm.argtypes = [ctypes.c_void_p]
revint.deletelab2kForm.restype = None
revint.deletelabk2cmykForm.argtypes = [ctypes.c_void_p]
revint.deletelabk2cmykForm.restype = None
revint.getLabExtmRng.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
revint.getLabExtmRng.restype = None
# Setup Gamut mapping algorithm from RevInt wrapper
revint.makegma.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32), ctypes.c_size_t]
revint.makegma.restype = ctypes.c_void_p
revint.map.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint16), PRINTPROGRESS]
revint.map.restype = None
revint.deletegma.argtypes = [ctypes.c_void_p]
revint.deletegma.restype = None


class cform:
    
    def __init__(self, profpath, luttype):
        profpath_p = profpath.encode('utf-8')
        luttype_p = luttype.encode('utf-8')
        self.address = ctypes.c_void_p(lcmswrap.makecform(profpath_p, luttype_p))
        
        if luttype.lower() == 'atob1':
            self.no_in_ch = 4
            self.no_out_ch = 3
        elif luttype.lower() == 'btoa1':
            self.no_in_ch = 3
            self.no_out_ch = 4
        
    def apply(self, inmat):
        if len(inmat.shape) == 1: #Expects rows, correct if one column vector is passed in
            inmat = inmat.reshape((1,inmat.shape[0]))
        outmat = np.ndarray(shape=(1, np.int(np.size(inmat, 0)*self.no_out_ch)), dtype=ctypes.c_double, order='F')
        inmat = inmat.reshape((1, np.int(np.size(inmat, 0)*self.no_in_ch)))
        outmat.setflags(write=1)
        inmat_p = inmat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        outmat_p = outmat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lcmswrap.applycform(self.address, inmat_p, outmat_p, ctypes.c_size_t(np.int(np.size(inmat, 1)/self.no_in_ch)))
        inmat = inmat.reshape((np.int(np.size(inmat, 1)/self.no_in_ch), self.no_in_ch))
        outmat = outmat.reshape((np.int(np.size(outmat, 1)/self.no_out_ch), self.no_out_ch))
        return outmat
        
    def delete(self):
        lcmswrap.deletecform(self.address)
        
        
class labk2cmykform:
    
    def __init__(self, filename, LabExtm):
        self.LabExtm = LabExtm
        CMYKlut = np.fromfile(filename, dtype=ctypes.c_double, count=-1, sep="")
        CMYKlut = CMYKlut.astype(ctypes.c_float) / 100
        CMYKlut_p = CMYKlut.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        self.address = ctypes.c_void_p(revint.makelabk2cmykForm(CMYKlut_p, ctypes.c_size_t(45)))
        
    def apply(self, labk):
        if len(labk.shape) == 1: #Expects rows, correct if one column vector is passed in
            labk = labk.reshape((1,4))
        labk_scl = labk.copy()
        labk_scl = labk_scl.astype(ctypes.c_double)
        labk_scl[:, 0:3] = labk_scl[:, 0:3] - self.LabExtm[0, :]
        labk_scl[:, 0:3] = labk_scl[:, 0:3] / (self.LabExtm[1, :] - self.LabExtm[0, :]) * 100
        labk_scl = labk_scl.reshape((1, np.int(np.size(labk_scl, 0)*4)))
        labk_scl_p = labk_scl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        cmyk = np.ndarray(shape=(1, np.size(labk_scl, 1)), dtype=ctypes.c_double, order='F')
        cmyk.setflags(write=1)
        cmyk_p = cmyk.ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
        revint.applylabk2cmykForm(self.address, labk_scl_p, cmyk_p, ctypes.c_size_t(np.int(np.size(labk_scl, 1) / 4)))
        cmyk = cmyk.reshape((np.int(np.size(cmyk, 1)/4), 4))
        return cmyk
    
    def delete(self):
        revint.deletelabk2cmykForm(self.address)
        

class gammapform:
    
    def __init__(self, filename):
        Lab_grid_LUT_gm_f = np.fromfile('bin/gm_lab.bin', dtype=ctypes.c_float, count=-1, sep='')
        Lab_grid_LUT_gm_f_p = Lab_grid_LUT_gm_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        self.address = ctypes.c_void_p(revint.makegammapForm(Lab_grid_LUT_gm_f_p, ctypes.c_uint(45)))
        
    def apply(self, labin):
        lab = labin.copy()
        lab[:, 0] = lab[:, 0] / 100 #Get in [0, 1] clut range
        lab[:, 1] = (lab[:, 1] + 128) / 255
        lab[:, 2] = (lab[:, 2] + 128) / 255
    
        lab = lab.reshape((1,np.int(np.size(lab,0)*3)))
        lab = lab.astype(ctypes.c_float)
        lab_gm = np.ndarray(shape=(1,np.size(lab,1)), dtype=ctypes.c_float, order='F')
        lab_p = lab.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        lab_gm_p = lab_gm.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        revint.applygammapForm(self.address, lab_p, lab_gm_p, ctypes.c_size_t(np.int(np.size(lab,1)/3)))
        lab = lab.reshape((np.int(np.size(lab,1)/3),3))
        
        lab[:, 0] = lab[:, 0] * 100 #Restore normal L*a*b* range
        lab[:, 1] = lab[:, 1] * 255 - 128
        lab[:, 2] = lab[:, 2] * 255 - 128
        
        lab_gm = lab_gm.reshape((np.int(np.size(lab_gm,1)/3),3))
        return lab_gm
        
    def delete(self):
        revint.deletegammapForm(self.address)