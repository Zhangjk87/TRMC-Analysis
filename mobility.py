# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Daniel's data analysis script for TRMC

#IMPORTANT: STORE PARAMETERS THAT ARE FOR ALL DATA POINTS TO BE ANALYZED IN A FILE CALLED sampleparams.py IN THE SAME FOLDER AS THIS SCRIPT

import numpy as np


import math

from cavityparams import *
from scipy.constants import *


def mobility(dP, laserPower, folder, P0):
    import sys
    sys.path.insert(0, folder)
    import sampleparams
    import resonanceparams
    #define absolute constants that won't change between data sets
    
    #correction for dP/P=ndV/V    
    dPCorrection = 1.42
    dP=dPCorrection*dP    
    
    pi=math.pi
    e0=epsilon_0#farads/meter
    q=e
    freq = c/(sampleparams.wavelength) #Hz
    photonenergy = h*freq
    #print(photonenergy)#in J
    #J/photon, this is for 532nm right now
    #For now what this will do is just put dP in here, and it will calculate a mobility. Next step is build this into the fitting program for one click data analysis
    
    if dP<0:
        sign = 1
    else:
        sign = -1
        
    print('sign = ', sign)
    
    #dP=0.0014#volts
    #laserPower=    #J/pulse on power meter
    
    spotarea=pi*(sampleparams.spotdiameter/2)**2#m^2
    #laserpower read in in joules
    #print(laserPower)
    I0=laserPower/photonenergy/spotarea
    print('I0 = '+str(I0))
    #print(format(I0, "e"))
    #careful with signs
    print('dP = ' + str(dP))
    K=sign*resonanceparams.Q*(1/np.sqrt(resonanceparams.R0)+sign*1)/(beta*e0*sampleparams.er*resonanceparams.f0*pi*L)#see absorption vs emission; I'm assume p0 is negative since that's what we get out of the detector, and a negative dP corresponds to an absorption with HP 462A amplifier
    print('K = ' + str(K))
    dG=sampleparams.illuminationFactor*dP/P0/K
    #print(dG)
    
    phimu=dG/beta/I0/sampleparams.Fa/q/1e-4
    print('mu = ' + str(phimu))
    return(phimu,I0, dPCorrection*sampleparams.illuminationFactor/P0/K/I0/sampleparams.Fa/q/beta/1e-4)        
