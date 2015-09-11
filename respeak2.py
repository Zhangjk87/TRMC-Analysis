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
    return( (R0 + np.power((2*(f-f0))/dW, 2))/(1+np.power((2*(f-f0))/dW,2)) + m*f + b)

def back_line(f, m, b):
    return(m*f + b)
    
#res_peak_folder = 'tempdata'
res_peak_folder = input('path containing res peak data: ')

os.chdir(res_peak_folder)
f='res peak.txt'

respeak = np.genfromtxt(f, delimiter=',', skip_header=7)
#respeak2=trim(respeak, 8.89e9, 8.925e9 )
#print(respeak)
#print(np.amax(respeak[:,1]))

#backgroundx = []
#for i in range(10):
#    j = np.argmin(respeak[:,0]) + i
#    print(j)
#    backgroundx.append(respeak[j,0])
#    
#for i in np.arange(10):
#    j = np.argmax(respeak[:,0:-1]) - i
#    print(j)
#    backgroundx.append(respeak[j,0])
#    
#backgroundy = []
#for i in range(10):
#    j = respeak[:,0] + i
#    backgroundy.append(respeak[j,1])
#    
#for i in np.arange(10):
#    j = np.argmax(respeak[:,0:-1]) - i
#    print(j)
#    backgroundy.append(respeak[j,1])
#    

#now finds min frequency automatically
guess = [.1, 0, 1e7, 1e-8, 0]

respeak[:,1] = np.power(respeak[:,1], 2)/50
respeak[:,1] = np.divide(respeak[:,1], np.amax(respeak[:,1])) #normalize resonance peak so maximum at 1
xdata = respeak[:,0]
ydata = respeak[:,1]

guess[1] = xdata[np.argmin(ydata)]

backgroundx = np.array(respeak[0:50,0])# + respeak[-1:-9, 0]])
backgroundy = np.array(respeak[0:50,1])# + respeak[-1:-9, 1]])
backgroundx = np.append(backgroundx, respeak[-50:, 0])
backgroundy = np.append(backgroundy, respeak[-50:, 1]) 
  
#start = np.argmin(respeak[:,0])
#end = np.argmax(respeak[:,0])
#print('start = ',start)
#print(respeak[start,0])
#print(respeak[start,1])
#print(respeak[end,0])
#print(respeak[end,1])

#print(backgroundx)
#print(backgroundy)

#background = []
#background = np.c_[backgroundx, backgroundy]
#print(background)
#background = background.T
#print(background)

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
            newfile.write('Q = '+str(Q)+str(R0)+'\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    print('Q = '+str(Q)+'\nR0 = '+str(R0) + '\nresponseTime = '+str(responseTime)+'\nf0 = '+str(popt[1]))
    #np.savetxt(res_peak_folder + '\\respeakfit.txt',fit, delimiter=',')
    #np.savetxt(res_peak_folder + '\\ydata.txt',ydata, delimiter=',')

    
except RuntimeError:
    print("Error - curve_fit failed")


plt.savefig('res peak.png')
plt.show()



