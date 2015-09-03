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
#from scipy import signal
import scipy.signal

def trim(mydata, minwavenumber, maxwavenumber):
    for counter, row in enumerate(mydata):
            if row[0] > minwavenumber:
                mydata = np.delete(mydata, np.s_[:counter], axis=0)
                break
    for counter, row in enumerate(mydata):
            if row[0] > maxwavenumber:
                mydata = np.delete(mydata, np.s_[counter:], axis=0)
                break
    return mydata

def singleExp(x, a, tau, offset):
    return a*np.exp(-(x)/tau) + offset
    
def doubleExp(x, a, t1, b, t2, offset):
    return a*np.exp(-(x)/t1) + b*np.exp(-(x)/t2) + offset

def tripleExp(x, a, t1, b, t2, c, t3, offset):
    return a*np.exp(-(x)/t1) + b*np.exp(-(x)/t2) + c*np.exp(-(x)/t3) + offset 

def readdata(f):
    print(f)
    print(f[f.index('_')+1:f.index('J')], 'J') #current laser power
    pulseenergy = float(f[f.index('nm_')+3:f.index('J')])
    p0=f[f.index('p0_')+3:f.index('V')] #current laser power
    print('P0 = ' + p0 + 'V')
    p0val=float(p0)
    TRMCdata = np.genfromtxt(f, delimiter=',', skip_header=4, usecols=(0,1))
    #tempTRMCdata = np.genfromtxt(f[:f.index('_p0')]+ '_' ++'_p0_' + p0 + 'V.csv', delimiter=',', skip_header=4, usecols=(0,2))
    #TRMCdata[:,1] = np.add(TRMCdata[:,1], tempTRMCdata[:,1])
    #TRMCdata[:,1] = np.divide(TRMCdata[:,1], numberofarrays)
    #print('numberofarrays='+str(numberofarrays))
    return(TRMCdata, pulseenergy, p0val)
    
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
     zerotime=np.argwhere(abs(array[:,0])<1e-17)
     datamean = np.mean(array[0:zerotime,1])
     array[:,1] = np.subtract(array[:,1], datamean)
     return array
       
def saveArray(f, array):
    np.savetxt(f, array, delimiter=',')

def fitSingle(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(singleExp, xdata.T, ydata.T, p0=guess)


        print('a = ' + str(popt[0]))
        print('t1 = '+ str(popt[1]))
        print('offset = '+str(popt[2]))
 
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
        return(0,0,0)

def fitDouble(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(doubleExp, xdata.T, ydata.T, p0=guess)
        print('a = ' + str(popt[0]))
        print('t1 = '+ str(popt[1]))
        print('b = ' + str(popt[2]))
        print('t2 = '+ str(popt[3]))
        print('offset = '+str(popt[4]))
 
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
        
def fitTriple(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(tripleExp, xdata.T, ydata.T, p0=guess)
        print('a = ' + str(popt[0]))
        print('t1 = '+ str(popt[1]))
        print('b = ' + str(popt[2]))
        print('t2 = '+ str(popt[3]))
        print('c = ' + str(popt[4]))
        print('t3 = '+ str(popt[5]))
        print('offset = '+str(popt[6]))
 
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
    if singleExponential==1:
        return (np.array([xdata, singleExp(xdata, *popt)]).T)
    elif singleExponential==2:
        return (np.array([xdata, doubleExp(xdata, *popt)]).T)
    elif singleExponential==3:
        return (np.array([xdata, tripleExp(xdata, *popt)]).T)
    else:
        raise RuntimeError('generateFitData')

def saveFitParams(f, params):
    with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        newfile.write('\n'.join(params))
        
def deconvolve(array, cavityLifetime):
    cavitydecayfft=np.fft.fft(np.exp(-array[:,0]/cavityLifetime)/sum(np.exp(-array[:,0]/cavityLifetime)))
    #datafft=np.fft.fft(signal.wiener(array[:,1],mysize=13))
    datafft=np.fft.fft(array[:,1])
    fftconv=np.divide(datafft,cavitydecayfft)
    deconv=np.fft.ifft(fftconv)
    deconvolvedData=np.array([array[:,0], deconv]).T
    return(deconvolvedData.real)

def binData(array, width):
    #Try to use averaging to make the data look nicer--http://stackoverflow.com/questions/21921178/binning-a-numpy-array
    ydata=array[:,1]
    ydataavg = ydata[:(ydata.size // width) * width].reshape(-1, width).mean(axis=1)
    xdata=array[:,0]
    xdataavg = xdata[:(xdata.size // width) * width].reshape(-1, width).mean(axis=1)
    return(np.array([xdataavg, ydataavg]).T)
    
def chargePerQD(I0, Fa, radius, packingFraction, thickness):
    #charges per QD = (Io*Fa*Volume_of_one_nanocrystal) / (packing_fraction*thickness)
    #packing_fraction = 0.6
    charge = (I0*Fa*(4/3*math.pi*radius**3)/(packingFraction*thickness))
    print('chargePerQD='+str(charge))
    return(charge)

#see http://stackoverflow.com/questions/25191620/creating-lowpass-filter-in-scipy-understanding-methods-and-units
def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y