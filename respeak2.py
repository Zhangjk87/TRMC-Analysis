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
import sys

def func(f, R0, f0, dW, m, b):
    return( (R0 + np.power((2*(f-f0))/dW, 2))/(1+np.power((2*(f-f0))/dW,2)) + m*f + b)

def back_line(f, m, b):
    return(m*f + b)
    
#res_peak_folder = 'tempdata'
    
try:
    res_peak_folder=sys.argv[1]
    print('res_peak_folder = ', res_peak_folder, '\n')
except:
    print('no folder command line argument\n')
    res_peak_folder = input('path containing TRMC data: ')

os.chdir(res_peak_folder)
file_name='respeak.txt'

respeak = np.genfromtxt(file_name, delimiter=',', skip_header=12)

# now finds min frequency automatically
guess = [.1, 0, 1e7, 1e-8, 0]

respeak[:,1] = np.power(respeak[:,1], 2)/50
respeak[:,1] = np.divide(respeak[:,1], np.amax(respeak[:,1]))  # normalize resonance peak so maximum at 1
xdata = respeak[:,0]
ydata = respeak[:,1]

guess[1] = xdata[np.argmin(ydata)]

backgroundx = np.array(respeak[0:50,0])  # + respeak[-1:-9, 0]])
backgroundy = np.array(respeak[0:50,1])  # + respeak[-1:-9, 1]])
backgroundx = np.append(backgroundx, respeak[-50:, 0])
backgroundy = np.append(backgroundy, respeak[-50:, 1])

line_guess = [-2e9, 0]

try:
    popt, pcov = curve_fit(back_line, np.array(backgroundx), np.array(backgroundy), p0=line_guess)
    print(popt)
    #print(pcov)
    linearfit=back_line(xdata, popt[0], popt[1])
    np.savetxt('linear_background.txt', linearfit)
    plt.plot(xdata, linearfit)
    
except RuntimeError:
    print("Error - curve_fit failed")


plt.plot(xdata, ydata)

#plt.show()
try:
    popt, pcov = curve_fit(func, xdata.T, ydata.T, p0=guess)
    print(popt)
    #print(pcov)
    fit=func(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    np.savetxt('respeakfit.txt', fit)
    plt.plot(xdata, fit)
    
    print(linearfit[np.argmin(respeak[:,1])])
    print('uncorrected R0 = ' + str(np.amin(respeak[:,1])))
    R0=np.amin(respeak[:,1])/linearfit[np.argmin(respeak[:,1])]
    #R0=np.amin(respeak[:,1])
    Q=popt[1]/popt[2]
    #print('Q = ' + str(Q))
    responseTime=(Q/math.pi / popt[1])
    #print('response time = ' + str(responseTime))
    if Q<0:
        print('\n\nWARNING: Q<0. DO NOT USE THESE PARAMETERS FOR FITTING. TRY AGAIN. resonanceparams.py not written.\n\n')
    else:
        with open('resonanceparams.py', 'w') as newfile:
            newfile.write('Q = '+str(Q)+ '\nR0 = ' +str(R0)+'\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    print('Q = '+str(Q)+'\nR0 = '+str(R0) + '\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    #np.savetxt(res_peak_folder + '\\respeakfit.txt',fit, delimiter=',')
    #np.savetxt(res_peak_folder + '\\ydata.txt',ydata, delimiter=',')

    
except RuntimeError:
    print("Error - curve_fit failed")


plt.savefig('res peak.png')
plt.show()



