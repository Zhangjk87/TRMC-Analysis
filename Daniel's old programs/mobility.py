# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Daniel's data analysis script for TRMC

#IMPORTANT: STORE PARAMETERS THAT ARE FOR ALL DATA POINTS TO BE ANALYZED IN A FILE CALLED sampleparams.py IN THE SAME FOLDER AS THIS SCRIPT

import numpy as np


import math
from sampleparams import *



def mobility(dP, laserPower):
    #define absolute constants that won't change between data sets
    beta = 0.9/0.4#--width/height of cavity
    L = 0.0913638#meters--cavity length
    pi=math.pi
    e0=8.85418782e-12#farads/meter
    q=1.602176565e-19
    spotdiameter=0.9#cm
    photonenergy=3.734e-16#mJ/photon, this is for 532nm right now
    #For now what this will do is just put dP in here, and it will calculate a mobility. Next step is build this into the fitting program for one click data analysis
    
    if negativesignal is True:
        sign = 1
    else:
        sign = -1
    
    
    #dP=0.0014#volts
    #laserPower=0.1#mJ/pulse on power meter
    
    spotarea=pi*(spotdiameter/2)**2#cm^2
    #laserpower read in in microjoules
    I0=laserPower*1e-3/photonenergy/spotarea
    print(I0)
    #print(format(I0, "e"))
    #careful with signs
    print(dP)
    K=sign*Q*(1+1/np.sqrt(R0))/(beta*e0*er*f0*pi*L)
    print(K)
    dG=dP/P0/K
    print(dG)
    
    phimu=0.619703*dG/beta/I0/Fa/q
    #print(phimu)
    return(phimu)    
    #ADD IN CORRECTION FOR ILLUMINATION FRACTION
    
