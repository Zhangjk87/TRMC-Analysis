# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 09:39:21 2015

@author: Daniel
"""

import numpy as np

import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
import math

def trim(mydata, minwavenumber, maxwavenumber):
    for counter, row in enumerate(mydata):
            if row[0] > maxwavenumber:
                mydata = np.delete(mydata, np.s_[counter:], axis=0)
                break
    for counter, row in enumerate(mydata):
        if row[0] > minwavenumber:
            mydata = np.delete(mydata, np.s_[:counter], axis=0)
            break

    return mydata

def func(f, R0, f0, dW,a):
    return( a*(R0 + np.power((2*(f-f0))/dW, 2))/(1+np.power((2*(f-f0))/dW,2)))

f='res peak.txt'

respeak = np.genfromtxt(f, delimiter=',', skip_header=7)
#respeak2=trim(respeak, 8.89e9, 8.925e9 )
#print(respeak)
xdata = respeak[:,0]
ydata = respeak[:,1]
guess = [.05, 8.95e9, 1e7,.5]

#plt.show()
popt, pcov = curve_fit(func, xdata.T, ydata.T, p0=guess)
print(popt)
#print(pcov)
plt.plot(respeak[:,0], respeak[:,1])
fit=func(xdata, popt[0], popt[1], popt[2],popt[3])
plt.plot(xdata, fit)

plt.savefig(f+'.png')

Q=popt[1]/popt[2]
R2=(func(popt[2], popt[0], popt[1], popt[2], popt[3]))

print('response time =' + str(Q/math.pi / popt[1]))

#print(R2)
#print(popt[3]*popt[1]+popt[4])
print('Q = ' + str(Q))
with open(f+'_params.txt', 'w') as newfile:
    newfile.write(''.join(str(popt.tolist())) + '\n\n' + str(Q))  
plt.show()