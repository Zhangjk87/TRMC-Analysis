# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:02:46 2015

@author: Daniel
"""

from scipy.constants import *
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
import matplotlib.pyplot as plt
import math

def singleexp(x, a, tau, offset):
    return a*np.exp(-(x)/tau) + offset
    
def doubleexp(x, a, t1, b, t2, offset):
    return a*np.exp(-(x)/t1) + b*np.exp(-(x)/t2) + offset

def readdata(f, numberofarrays, filenameconditions):
    if filenameconditions is True:
        #print(f)
        print(f[f.index('_')+1:f.index('J')], 'J') #current laser power
        pulseenergy = format(float(f[f.index('_')+1:f.index('J')]))
        TRMCdata = np.genfromtxt(f, delimiter=',', skip_header=4, usecols=(0,2))
        for i in range(2,numberofarrays+1):
            #print(f[:f.index('_1.csv')]+'_' + str(i) +'.csv')
            tempTRMCdata = np.genfromtxt(f[:f.index('_1.csv')]+'_' + str(i) +'.csv', delimiter=',', skip_header=4, usecols=(0,2))
            TRMCdata[:,1] = np.add(TRMCdata[:,1], tempTRMCdata[:,1])
        TRMCdata[:,1] = np.divide(TRMCdata[:,1], numberofarrays)
        return {'TRMCdata':TRMCdata, 'pulseenergy':pulseenergy}
    
    def findmaxormin(array, maxormin):
        if (abs(np.amax(array[:,1])) > abs(np.amin(array[:,1]))):
            return np.argmax(array[:,1])
        elif (abs(np.amax(array[:,1])) < abs(np.amin(array[:,1]))):
            return np.argin(array[:,1])
        else:
            print("Error in findmaxormin(). Max not less than or greater than min array value.")
            raise RuntimeError('findmaxormin')



        
        
        #datamean = np.mean(TRMCdata[0:5001,1])
        #TRMCdata[:,1] = np.subtract(TRMCdata[:,1], datamean)
        ## np.savetxt('corrected.txt', TRMCdata, delimiter=',')

        #decay = TRMCdata[maxindex:]
        
