# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Daniel's data analysis script for TRMC

#IMPORTANT: STORE PARAMETERS THAT ARE FOR ALL DATA POINTS TO BE ANALYZED IN A FILE CALLED sampleparams.py IN THE SAME FOLDER AS THIS SCRIPT

import numpy as np


import math
from tempdata.sampleparams import *
from tempdata.resonanceparams import *
from cavityparams import *
from scipy.constants import *


def mobility(dP, laserPower):
    #define absolute constants that won't change between data sets

    pi=math.pi
    e0=epsilon_0#farads/meter
    q=e
    freq = c/(wavelength*1e-9) #Hz
    photonenergy = h*freq*1000
    print(photonenergy)#in mJ
    #mJ/photon, this is for 532nm right now
    #For now what this will do is just put dP in here, and it will calculate a mobility. Next step is build this into the fitting program for one click data analysis
    
    if dP<0:
        sign = -1
    else:
        sign = 1
    
    
    #dP=0.0014#volts
    #laserPower=0.1#mJ/pulse on power meter
    
    spotarea=pi*(spotdiameter/2)**2#cm^2
    #laserpower read in in microjoules
    print(laserPower)
    I0=laserPower*1000/photonenergy/spotarea
    print('I0='+str(I0))
    #print(format(I0, "e"))
    #careful with signs
    print('dP = ' + str(dP))
    K=-sign*Q*(1+sign*1/np.sqrt(R0))/(beta*e0*er*f0*pi*L)#see absorption vs emission; I'm assume p0 is negative since that's what we get out of the detector, and a positive dP is "less negative" power so corresponds to an absorption
    print('K=' + str(K))
    dG=dP/P0/K
    #print(dG)
    
    phimu=illuminationFactor*dG/beta/I0/Fa/q
    print('mu = ' + str(phimu))
    return(phimu)        
