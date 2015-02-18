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

def singleExp(x, a, tau, offset):
    return a*np.exp(-(x)/tau) + offset
    
def doubleExp(x, a, t1, b, t2, offset):
    return a*np.exp(-(x)/t1) + b*np.exp(-(x)/t2) + offset

def readdata(f, numberofarrays):
    #print(f)
    print(f[f.index('_')+1:f.index('J')], 'J') #current laser power
    pulseenergy = format(float(f[f.index('_')+1:f.index('J')]))
    TRMCdata = np.genfromtxt(f, delimiter=',', skip_header=4, usecols=(0,2))
    for i in range(2,numberofarrays+1):
        #print(f[:f.index('_1.csv')]+'_' + str(i) +'.csv')
        tempTRMCdata = np.genfromtxt(f[:f.index('_1.csv')]+'_' + str(i) +'.csv', delimiter=',', skip_header=4, usecols=(0,2))
        TRMCdata[:,1] = np.add(TRMCdata[:,1], tempTRMCdata[:,1])
    TRMCdata[:,1] = np.divide(TRMCdata[:,1], numberofarrays)
    return (TRMCdata, pulseenergy)
    
def findmaxormin(array):
    if (abs(np.amax(array[:,1])) > abs(np.amin(array[:,1]))):
        return np.argmax(array[:,1])
    elif (abs(np.amax(array[:,1])) < abs(np.amin(array[:,1]))):
        return np.argmin(array[:,1])
    else:
        print("Error in findmaxormin(). Max not less than or greater than min array value.")
        raise RuntimeError('findmaxormin')
        #return None

def subtractOffset(array):
     datamean = np.mean(array[0:5001,1])
     array[:,1] = np.subtract(array[:,1], datamean)
     return array
       
def saveArray(f, array):
    np.savetxt(f, array, delimiter=',')

def fitSingle(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(singleExp, xdata.T, ydata.T, p0=guess)


        print('a=' + str(popt[0]))
        print('t1='+ str(popt[1]))
        print('offset='+str(popt[2]))
 
        #np.savetxt(f+'_fit.csv', popt, delimiter = ',')
        #np.savetxt(f+'_fit_cov.csv', popt, delimiter = ',')
        
        # from http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

        alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
    
        n = len(ydata)    # number of data points
        p = len(popt) # number of parameters
        
        dof = max(0, n - p) # number of degrees of freedom
        
        # student-t value for the dof and confidence level
        tval = t.ppf(1.0-alpha/2., dof)
        params = []
        for i, p,var in zip(range(n), popt, np.diag(pcov)):
            sigma = var**0.5
            params.append('p{0}: {1} [{2}  {3}]'.format(i, p,p - sigma*tval,p + sigma*tval))
        #print(params)
        #with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        #    newfile.write('\n'.join(params))
        return(popt, pcov, params)
    except RuntimeError:
        print("Error - curve_fit failed")

def fitDouble(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(doubleExp, xdata.T, ydata.T, p0=guess)

        print('a=' + str(popt[0]))
        print('t1='+ str(popt[1]))
        print('b=' + str(popt[2]))
        print('t2='+ str(popt[3]))
        print('offset='+str(popt[4]))
 
        #np.savetxt(f+'_fit.csv', popt, delimiter = ',')
        #np.savetxt(f+'_fit_cov.csv', popt, delimiter = ',')
        
        # from http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

        alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
    
        n = len(ydata)    # number of data points
        p = len(popt) # number of parameters
        
        dof = max(0, n - p) # number of degrees of freedom
        
        # student-t value for the dof and confidence level
        tval = t.ppf(1.0-alpha/2., dof)
        params = []
        for i, p,var in zip(range(n), popt, np.diag(pcov)):
            sigma = var**0.5
            params.append('p{0}: {1} [{2}  {3}]'.format(i, p,p - sigma*tval,p + sigma*tval))
        #print(params)
        #with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        #    newfile.write('\n'.join(params))
        return(popt, pcov, params)
    except RuntimeError:
        print("Error - curve_fit failed")
        
def generateFitData(singleExponential, xdata, *popt):
    if singleExponential==True:
        return (np.array([xdata, singleExp(xdata, *popt)]).T)
    elif singleExponential==False:
        return (np.array([xdata, doubleExp(xdata, *popt)]).T)
    else:
        raise RuntimeError('generateFitData')

def saveFitParams(f, params):
    with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        newfile.write('\n'.join(params))
        
def deconvolve(array, cavityLifetime):
    cavitydecayfft=np.fft.fft(np.exp(-array[:,0]/cavityLifetime)/sum(np.exp(-array[:,0]/cavityLifetime)))
    datafft=np.fft.fft(array[:,1])
    fftconv=np.divide(datafft,cavitydecayfft)
    deconv=np.fft.ifft(fftconv)
    deconvolvedData=np.array([array[:,0], deconv]).T
    return(deconvolvedData)

def binData(array, width):
    #Try to use averaging to make the data look nicer--http://stackoverflow.com/questions/21921178/binning-a-numpy-array
    ydata=array[:,1]
    ydataavg = ydata[:(ydata.size // width) * width].reshape(-1, width).mean(axis=1)
    xdata=array[:,0]
    xdataavg = xdata[:(xdata.size // width) * width].reshape(-1, width).mean(axis=1)
    return(np.array([xdataavg, ydataavg]).T)