#Transforms only replaceable colors

import numpy as np
import scipy as sp
from PIL import Image
import time
import timeit
import ctypes
from ctypes import cdll
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import gcrpywrap as gw

# Define auxiliary functions
PRINTPROGRESS = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_int)
def PrintProgress(current, total):
    print(str(current + 1) + ' / ' + str(total))

def lutindex(l, k, j, nopts):
		return nopts**2*l + nopts*k + j

class wmgcr:
  def __init__(self, prof_name):
      
    #Define min CMY and maxK for replacement
    #E.g. don't replace colors where CMY<10
    self.repCMYmin=0
    self.repKmax=100
    #Set total ink limit
    self.TIL = 400
    
    # Setup input and output folders
    self.prof_name = prof_name

	# Setup CMYK->Lab transform
    self.profpath = Path('profiles/' + self.prof_name)
    if not self.profpath.exists():
        raise FileNotFoundError(self.profpath.name + ' does not exist in the profiles/ folder')
        pass
    else:
        self.cfAtoB = gw.cform(self.profpath.resolve().as_posix(), 'AtoB1')

	# Permutations of CMYK values
    import itertools
    x = np.linspace(0, 100, 16)
    CMYK_grid = np.array([p for p in itertools.product(x, repeat=4)], dtype=ctypes.c_double, order='F')
    Lab_grid = self.cfAtoB.apply(CMYK_grid)

	#Get extreme Lab values (LUT endpoints)
    self.LabExtm = np.ndarray(shape=(2, 3), dtype=ctypes.c_double, order='F')
    self.LabExtm[0, :] = np.min(Lab_grid, 0)
    self.LabExtm[1, :] = np.max(Lab_grid, 0)

	#Setup gamut mapping form
    self.gmform = gw.gammapform('bin/gm_lab.bin')

	# Setup LabK->CMYK transform
    self.lk2cform = gw.labk2cmykform('bin/cmyklut.bin', self.LabExtm)

	#Setup Lab to Kmin Kmax transform
    Klut = np.fromfile('bin/kminmax.bin', dtype=ctypes.c_double, count=-1, sep="")
    Klut = Klut.reshape(np.int(np.size(Klut,0)/2), 2)
    from scipy.interpolate import RegularGridInterpolator
    Lnode = np.linspace(self.LabExtm[0, 0], self.LabExtm[1, 0], 45)
    anode = np.linspace(self.LabExtm[0, 1], self.LabExtm[1, 1], 45)
    bnode = np.linspace(self.LabExtm[0, 2], self.LabExtm[1, 2], 45)
    Vmin = np.zeros((45, 45, 45))
    Vmax = np.zeros((45, 45, 45))
    for l in range(45):
    	for k in range(45):
    		for j in range(45):
    			Vmin[l, k, j] = Klut[np.int(lutindex(l, k, j, 45)) ,0]
    			Vmax[l, k, j] = Klut[np.int(lutindex(l, k, j, 45)) ,1]

    self.Kminfn = RegularGridInterpolator((Lnode, anode, bnode), Vmin)
    self.Kmaxfn = RegularGridInterpolator((Lnode, anode, bnode), Vmax)

  # Define image transform function
  def transformImage(self, im_orig, im_wm, rpl):
    #Get image data
    pix = im_orig
    pixwm = im_wm
    rpl = rpl.reshape((np.int(rpl.shape[0]*rpl.shape[1])))
    pix_dbl = np.ndarray.astype(pix, dtype=ctypes.c_double)
    pix_dbl = pix_dbl.reshape((np.size(pix, 0)*np.size(pix, 1), np.size(pix, 2)))
    I_rpl = rpl
    pix_rpl = pix_dbl[I_rpl,:]
    pix_unq, Iunq, Iunq_inv = np.unique(pix_rpl, return_inverse=True, return_index=True, axis=0)
    pix_unq = pix_unq / 2.55
    #Tansform CMYK->Lab
    pix_Lab = self.cfAtoB.apply(pix_unq)
  
    #Gamut mapping using form
    pix_Lab_gm = self.gmform.apply(pix_Lab)
    
    #Get Kmin and Kmax
    Kminmax = np.ndarray(shape=(np.size(pix_Lab_gm, 0), 2), dtype=ctypes.c_double, order='F')
    Kminmax[:,0] = self.Kminfn(pix_Lab_gm)
    Kminmax[:,1] = self.Kmaxfn(pix_Lab_gm)
    
    
    #Try function estimated Kmin and iteratively find real Kmin
    LabK = np.ndarray(shape=(np.size(pix_Lab_gm, 0), 4), dtype=ctypes.c_double, order='F')
    LabK[:, 0:3] = pix_Lab_gm
    LabK[:, 3] = Kminmax[:,0]
    CMYK_BG = self.lk2cform.apply(LabK)
    acov =  np.round(np.sum(CMYK_BG, 1), 2)  #Area coverage
    Indsfix = np.logical_or(CMYK_BG[:,0] > 100, np.logical_or(CMYK_BG[:,1] > 100, CMYK_BG[:,2] > 100))
    Indicesfix = np.where(Indsfix==True)  
    Kprev = Kminmax[Indsfix,0]
    Knext = Kminmax[Indsfix,1]
    
    loopcount=0
    while (np.size(np.where(Indsfix==True)) > 0 and loopcount<100):
        print(np.size(np.where(Indsfix==True)))

        Kcur = (Knext + Kprev) / 2
        LabKcur = np.ndarray(shape=(np.size(Kprev), 4), dtype=ctypes.c_double, order='F')
        LabKcur[:, 0:3] = LabK[Indsfix, 0:3]
        LabKcur[:,3] = Kcur
        CMYKcur = self.lk2cform.apply(LabKcur)
        #Round to second decimal
        CMYKcurrnd = CMYKcur.copy()
        acov =  np.round(np.sum(CMYKcurrnd, 1), 2)  #Area coverage
        CMYKcurrnd = np.round(CMYKcurrnd, 2)
        #Condition 1 - C,M,Y > 100; Increase K
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] > 100, np.logical_or(CMYKcurrnd[:,1] > 100, CMYKcurrnd[:,2] > 100))
        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 2 - C,M,Y > 0 AND acov > TIL; Increase K (acov - area coverage; TIL - total ink limit)
#        Indsfixcur = np.logical_and(np.logical_and(CMYKcur[:,0] > 0, np.logical_and(CMYKcur[:,1] > 0, CMYKcur[:,2] > 0)), acov > TIL)
#        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
#        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 3 - C,M,Y < 100 & acov < TIL
        Indsfixcur = np.logical_and(CMYKcurrnd[:,0] < 100, np.logical_and(CMYKcurrnd[:,1] < 100, CMYKcurrnd[:,2] < 100))
        Knext[Indsfixcur]=np.min(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Conditon 4 - C,M,Y == 100 OR acov = TIL; Real Kmin found
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] == 100, np.logical_or(CMYKcurrnd[:,1] == 100, CMYKcurrnd[:,2] == 100))
        Kminmax[Indicesfix[0][Indsfixcur], 0] = Kcur[Indsfixcur]
        Indsfix[Indicesfix[0][Indsfixcur]] = False
        Indicesfix = np.where(Indsfix==True)
        Kprev=np.delete(Kprev, np.where(Indsfixcur==True)[0])
        Knext=np.delete(Knext, np.where(Indsfixcur==True)[0])
        loopcount+=1
    
    
    #Try function estimated Kmax and iteratively find real Kmax
    LabK[:, 3] = Kminmax[:,1]
    CMYK_BG = self.lk2cform.apply(LabK)
    Indsfix = np.logical_or(CMYK_BG[:,0] < 0, np.logical_or(CMYK_BG[:,1] < 0, CMYK_BG[:,2] < 0))
    Indicesfix = np.where(Indsfix==True)
    Kprev = Kminmax[Indsfix,1]
    Knext = Kminmax[Indsfix,0] - 10 #Kminfn is inaccurate so take lower values to ensure convergence
                                    #Taking too large Kmin from Kminfn will prevent convergence
    loopcount=0
    while (np.size(np.where(Indsfix==True)) > 0 and loopcount<100):
        print(np.size(np.where(Indsfix==True)))

        Kcur = (Knext + Kprev) / 2
        LabKcur = np.ndarray(shape=(np.size(np.where(Indsfix==True)), 4), dtype=ctypes.c_double, order='F')
        LabKcur[:, 0:3] = LabK[Indsfix, 0:3]
        LabKcur[:,3] = Kcur
        CMYKcur = self.lk2cform.apply(LabKcur)
        #Round to second decimal
        CMYKcurrnd = CMYKcur.copy()
        CMYKcurrnd = np.round(CMYKcurrnd, 2)
        #Condition 1 - C,M,Y < 0; Decrease K
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] < 0, np.logical_or(CMYKcurrnd[:,1] < 0, CMYKcurrnd[:,2] < 0))
        Knext[Indsfixcur]=np.min(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 2 - C,M,Y > 0; Increase K
        Indsfixcur = np.logical_and(CMYKcurrnd[:,0] > 0, np.logical_and(CMYKcurrnd[:,1] > 0, CMYKcurrnd[:,2] > 0))
        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 3 - C,M,Y == 0
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] == 0, np.logical_or(CMYKcurrnd[:,1] == 0, CMYKcurrnd[:,2] == 0))
        Kminmax[Indicesfix[0][Indsfixcur], 1] = Kcur[Indsfixcur]
        Indsfix[Indicesfix[0][Indsfixcur]] = False
        Indicesfix = np.where(Indsfix==True)
        Kprev=np.delete(Kprev, np.where(Indsfixcur==True)[0])
        Knext=np.delete(Knext, np.where(Indsfixcur==True)[0])
        loopcount+=1
    
    #Black generation
    pixwm_dbl = np.ndarray.astype(pixwm, dtype=ctypes.c_double)
    pixwm_dbl = pixwm_dbl.reshape((np.size(pix, 0)*np.size(pix, 1), np.size(pix, 2)))
    #K = pixwm_dbl[Iunq, 3]/2.55
    K = pixwm_dbl[I_rpl,:][Iunq, 3]/2.55
    
    #Transform LabK->CMYK
    LabK = np.ndarray(shape=(np.size(pix_Lab_gm, 0), 4), dtype=ctypes.c_double, order='F')
    LabK[:, 0:3] = pix_Lab_gm
    LabK[:, 3] = K
    CMYK_BG = self.lk2cform.apply(LabK)    
    CMYK_BG[CMYK_BG < 0] = 0
    CMYK_BG[CMYK_BG > 100] = 100
    pix_new = pix_dbl.copy()
    pix_new[I_rpl,:] = CMYK_BG[Iunq_inv,:] * 2.55
    pix_new_ui8 = (pix_new).astype(np.uint8)
    im_new = Image.fromarray(pix_new_ui8.reshape((np.size(pix, 0), np.size(pix, 1), 4)))
    im_new.mode = 'CMYK'    
    rpl = rpl.reshape((np.size(pix, 0), np.size(pix, 1)))
    return im_new

  def isReplaceable(self, im, Krng, repCMYmin, repKmax=100):
    im = np.ndarray.astype(im, dtype=ctypes.c_double)
    pix=im
    pix_dbl = np.ndarray.astype(pix, dtype=ctypes.c_double)
    pix_dbl = pix_dbl.reshape((np.size(pix, 0)*np.size(pix, 1), np.size(pix, 2)))
    I_rpl = np.logical_and(np.logical_and(np.logical_and(pix_dbl[:,0]/2.55>self.repCMYmin, pix_dbl[:,1]/2.55>self.repCMYmin), pix_dbl[:,2]/2.55>self.repCMYmin), pix_dbl[:,3]/2.55<self.repKmax)
    pix_rpl = pix_dbl[I_rpl,:]
    pix_unq, Iunq, Iunq_inv = np.unique(pix_rpl, return_inverse=True, return_index=True, axis=0)
    pix_unq = pix_unq / 2.55
    #Tansform CMYK->Lab
    pix_Lab = self.cfAtoB.apply(pix_unq)
  
    #Gamut mapping using form
    pix_Lab_gm = self.gmform.apply(pix_Lab)
    
    #Get Kmin and Kmax
    Kminmax = np.ndarray(shape=(np.size(pix_Lab_gm, 0), 2), dtype=ctypes.c_double, order='F')
    Kminmax[:,0] = self.Kminfn(pix_Lab_gm)
    Kminmax[:,1] = self.Kmaxfn(pix_Lab_gm)
    
    
    #Try function estimated Kmin and iteratively find real Kmin
    LabK = np.ndarray(shape=(np.size(pix_Lab_gm, 0), 4), dtype=ctypes.c_double, order='F')
    LabK[:, 0:3] = pix_Lab_gm
    LabK[:, 3] = Kminmax[:,0]
    CMYK_BG = self.lk2cform.apply(LabK)
    acov =  np.round(np.sum(CMYK_BG, 1), 2)  #Area coverage
    Indsfix = np.logical_or(CMYK_BG[:,0] > 100, np.logical_or(CMYK_BG[:,1] > 100, CMYK_BG[:,2] > 100))
    Indicesfix = np.where(Indsfix==True)  
    Kprev = Kminmax[Indsfix,0]
    Knext = Kminmax[Indsfix,1]
    
    loopcount=0
    while (np.size(np.where(Indsfix==True)) > 0 and loopcount<100):
        print(np.size(np.where(Indsfix==True)))

        Kcur = (Knext + Kprev) / 2
        LabKcur = np.ndarray(shape=(np.size(Kprev), 4), dtype=ctypes.c_double, order='F')
        LabKcur[:, 0:3] = LabK[Indsfix, 0:3]
        LabKcur[:,3] = Kcur
        CMYKcur = self.lk2cform.apply(LabKcur)
        #Round to second decimal
        CMYKcurrnd = CMYKcur.copy()
        acov =  np.round(np.sum(CMYKcurrnd, 1), 2)  #Area coverage
        CMYKcurrnd = np.round(CMYKcurrnd, 2)
        #Condition 1 - C,M,Y > 100; Increase K
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] > 100, np.logical_or(CMYKcurrnd[:,1] > 100, CMYKcurrnd[:,2] > 100))
        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 2 - C,M,Y > 0 AND acov > TIL; Increase K (acov - area coverage; TIL - total ink limit)
#        Indsfixcur = np.logical_and(np.logical_and(CMYKcur[:,0] > 0, np.logical_and(CMYKcur[:,1] > 0, CMYKcur[:,2] > 0)), acov > TIL)
#        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
#        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 3 - C,M,Y < 100 & acov < TIL
        Indsfixcur = np.logical_and(CMYKcurrnd[:,0] < 100, np.logical_and(CMYKcurrnd[:,1] < 100, CMYKcurrnd[:,2] < 100))
        Knext[Indsfixcur]=np.min(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Conditon 4 - C,M,Y == 100 OR acov = TIL; Real Kmin found
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] == 100, np.logical_or(CMYKcurrnd[:,1] == 100, CMYKcurrnd[:,2] == 100))
        Kminmax[Indicesfix[0][Indsfixcur], 0] = Kcur[Indsfixcur]
        Indsfix[Indicesfix[0][Indsfixcur]] = False
        Indicesfix = np.where(Indsfix==True)
        Kprev=np.delete(Kprev, np.where(Indsfixcur==True)[0])
        Knext=np.delete(Knext, np.where(Indsfixcur==True)[0])
        loopcount+=1
    
    
    #Try function estimated Kmax and iteratively find real Kmax
    LabK[:, 3] = Kminmax[:,1]
    CMYK_BG = self.lk2cform.apply(LabK)
    Indsfix = np.logical_or(CMYK_BG[:,0] < 0, np.logical_or(CMYK_BG[:,1] < 0, CMYK_BG[:,2] < 0))
    Indicesfix = np.where(Indsfix==True)
    Kprev = Kminmax[Indsfix,1]
    Knext = Kminmax[Indsfix,0] - 10 #Kminfn is inaccurate so take lower values to ensure convergence
                                    #Taking too large Kmin from Kminfn will prevent convergence
    loopcount=0
    while (np.size(np.where(Indsfix==True)) > 0 and loopcount<100):
        print(np.size(np.where(Indsfix==True)))

        Kcur = (Knext + Kprev) / 2
        LabKcur = np.ndarray(shape=(np.size(np.where(Indsfix==True)), 4), dtype=ctypes.c_double, order='F')
        LabKcur[:, 0:3] = LabK[Indsfix, 0:3]
        LabKcur[:,3] = Kcur
        CMYKcur = self.lk2cform.apply(LabKcur)
        #Round to second decimal
        CMYKcurrnd = CMYKcur.copy()
        CMYKcurrnd = np.round(CMYKcurrnd, 2)
        #Condition 1 - C,M,Y < 0; Decrease K
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] < 0, np.logical_or(CMYKcurrnd[:,1] < 0, CMYKcurrnd[:,2] < 0))
        Knext[Indsfixcur]=np.min(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 2 - C,M,Y > 0; Increase K
        Indsfixcur = np.logical_and(CMYKcurrnd[:,0] > 0, np.logical_and(CMYKcurrnd[:,1] > 0, CMYKcurrnd[:,2] > 0))
        Knext[Indsfixcur]=np.max(np.array([Kprev[Indsfixcur], Knext[Indsfixcur]]).transpose(), 1)
        Kprev[Indsfixcur]=Kcur[Indsfixcur]
        #Condition 3 - C,M,Y == 0
        Indsfixcur = np.logical_or(CMYKcurrnd[:,0] == 0, np.logical_or(CMYKcurrnd[:,1] == 0, CMYKcurrnd[:,2] == 0))
        Kminmax[Indicesfix[0][Indsfixcur], 1] = Kcur[Indsfixcur]
        Indsfix[Indicesfix[0][Indsfixcur]] = False
        Indicesfix = np.where(Indsfix==True)
        Kprev=np.delete(Kprev, np.where(Indsfixcur==True)[0])
        Knext=np.delete(Knext, np.where(Indsfixcur==True)[0])
        loopcount+=1
    rpl = np.ndarray(shape=(pix.shape[0]*pix.shape[1], 1), dtype=bool, order='F')
    rpl.fill(False)
    I_rng = np.abs(Kminmax[:,0]-Kminmax[:,1]) >= Krng
    rpl[I_rpl,0] = I_rng[Iunq_inv]
    rpl = rpl.reshape((pix.shape[0], pix.shape[1]))
    return rpl

  # Destructor 
  def delete(self): 
    self.cfAtoB.delete()
    self.gmform.delete()
    self.lk2cform.delete()
