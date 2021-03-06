# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 09:39:21 2015

@author: Daniel
"""

import numpy as np

import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
#from scipy.stats.distributions import  t
import math

def func(f, R0, f0, dW, m, b):
    return( (R0 + np.power((2*(f-f0))/dW, 2))/(1+np.power((2*(f-f0))/dW,2)) +m*f +b)

#res_peak_folder = 'tempdata'
res_peak_folder = input('path containing res peak data: ')

os.chdir(res_peak_folder)

f='res peak.txt'

respeak = np.genfromtxt(f, delimiter=',', skip_header=7)
#respeak2=trim(respeak, 8.89e9, 8.925e9 )
#print(respeak)
#print(np.amax(respeak[:,1]))
respeak[:,1] = np.power(respeak[:,1], 2)/50
respeak[:,1] = np.divide(respeak[:,1], np.amax(respeak[:,1])) #normalize resonance peak so maximum at 1
xdata = respeak[:,0]
ydata = respeak[:,1]

guess = [.1, 9.015e9, 1e7, 1e-8, 0]
plt.plot(xdata, ydata)


#plt.show()
try:
    popt, pcov = curve_fit(func, xdata.T, ydata.T, p0=guess)
    print(popt)
    #print(pcov)
    fit=func(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    np.savetxt('respeakfit.txt', fit)
    plt.plot(xdata, fit)
    
    R0=np.amin(respeak[:,1])
    Q=popt[1]/popt[2]
    #print('Q = ' + str(Q))
    responseTime=(Q/math.pi / popt[1])
    #print('response time = ' + str(responseTime))
    if Q<0:
        print('\n\nWARNING: Q<0. DO NOT USE THESE PARAMETERS FOR FITTING. TRY AGAIN. resonanceparams.py not written.\n\n')
    else:
        with open('resonanceparams.py', 'w') as newfile:
            newfile.write('Q = '+str(Q)+'\nR0 = '+str(R0)+'\nR0 = '+str(R0)+'\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    print('Q = '+str(Q)+'\nR0 = '+str(R0)+'\nR0 = '+str(R0)+'\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    #np.savetxt(res_peak_folder + '\\respeakfit.txt',fit, delimiter=',')
    #np.savetxt(res_peak_folder + '\\ydata.txt',ydata, delimiter=',')

    
except RuntimeError:
    print("Error - curve_fit failed")


plt.savefig('res peak.png')
plt.show()



