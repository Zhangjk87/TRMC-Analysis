# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 10:46:28 2014

@author: Daniel
"""

# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t

import matplotlib.pyplot as plt
import os
import math

import csv

#put mobility.py file in the same folder
from mobility import *
#from sampleparams import *

# some code from http://python4esac.github.io/fitting/examples1d.html
#def func(x, a, x0, sigma, offset):
#    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + offset

def func(x, a, tau, offset):
    return a*np.exp(-(x)/tau) + offset

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
j=0
lifetimelist=[]
powerlist = []
mobilitylist = []
chargelist = []
for f in os.listdir("."):
    if "p0" not in f and "res" not in f and '.png' not in f and 'fit' not in f and ".csv" in f and 'J_' in f and '_1.csv' in f and not 'averaged' in f:
        print(f)
        datatofit = np.genfromtxt(f, delimiter=',', skip_header=4, usecols=(0,2))
        numberofarrays=2
        for i in range(2,numberofarrays+1):
            print(f[:f.index('_1.csv')]+'_' + str(i) +'.csv')
            tempdatatofit = np.genfromtxt(f[:f.index('_1.csv')]+'_' + str(i) +'.csv', delimiter=',', skip_header=4, usecols=(0,2))
            datatofit[:,1]= np.add(datatofit[:,1], tempdatatofit[:,1])
        datatofit[:,1]=np.divide(datatofit[:,1], numberofarrays)
        #datatofit = np.genfromtxt(f, delimiter=',', skip_header=4, usecols=(0,2))
        baseline =  trim(datatofit, -9e9, 0)
        averagepower=np.average(datatofit[:,1])
        datatofit2 = trim(datatofit, 0, 9e9)
        datatofit3=(np.copy(datatofit2))
        datatofit3[:,0]=datatofit3[:,0]/1e-9
        datatofit3[:,1]=-(datatofit3[:,1]-averagepower)/P0
        np.savetxt(f+'averaged.csv', datatofit3, delimiter=',')
        maxindex=np.argmax(datatofit2[:,1])
        maxval=datatofit2[maxindex,1]
        print(maxindex)
        print(maxval)
        datatofit2 = trim(datatofit2, datatofit2[maxindex, 0], 9e9)
        #print(datatofit2)
        xdata = datatofit2[:,0]
        #print(xdata)
        ydata = datatofit2[:,1]
        guess = [2e-4, 10e-8, 0]
        #plt.plot(xdata, ydata)
        #plt.show()
        try:
            popt, pcov = curve_fit(func, xdata.T, ydata.T, p0=guess)
        #print(popt)
        #peak = popt[3]+popt[0]
        #print(i)
        #i=i//2 #for outputting datafromfit
        
        #datafromfit[i,0] = peak
        #datafromfit[i,1] = popt[1]
        #print(i)
        #print(datafromfit[i,0])
        #print(datafromfit[i,1])
            print(popt[0])
            print(popt[1])
            print(popt[2])
 
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
            print(params)
            with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
                newfile.write('\n'.join(params))                    
            fit = func(xdata, popt[0], popt[1], popt[2])
            plt.figure(0)            
            plt.plot(xdata, fit, 'g')
            plt.figure(1)            
            plt.plot(xdata, fit, 'g')
            laserpowernum=float(f[f.index(' _')+2:f.index('J')])*0.125167/1e-6#in microJ
            laserpower=format(laserpowernum, '.0f')#in microJ
            powerlist.append(laserpower)
            lifetimelist.append(popt[1])
            
            
            spotdiameter=0.9#cm
            photonenergy=3.734e-16#mJ/photon, this is for 532nm right now
            # fit mobility  
            
            spotarea=math.pi*(spotdiameter/2)**2#cm^2
    #laserpower read in in microjoules
            I0=laserpowernum*1e-3/photonenergy/spotarea
            packing_fraction=0.6
            thickness=40e-7
            Diameter=3.66e-7
            Volume = (4/3)*math.pi*(Diameter/2)**3            
            
            charge_per_QD = (I0*Fa*Volume)/(packing_fraction*thickness)
            chargelist.append(charge_per_QD)

            mobval = mobility(maxval-averagepower,laserpowernum)
            mobilitylist.append(mobval)            
             
            print(mobval)
            j=j+1
            print('\n\n')
            
        except RuntimeError:
            print("Error - curve_fit failed")
            #lifetimelist.append(0)


        plt.figure(1)        
        plt.plot(datatofit[:,0], datatofit[:,1], label=str(laserpower))
        plt.figure(0)        
        plt.plot(datatofit[:,0], datatofit[:,1], 'b')
        plt.savefig(f + '.png')
        #plt.show()
        plt.close()
        
        

        
plt.figure(1)
plt.xlim(0,5e-7)
plt.legend()
plt.savefig('all.png')
        #plt.show()
plt.close()
plt.figure(2)
plt.plot(powerlist,lifetimelist, 'bo')
plt.xlabel('laser power (uJ)')
plt.ylabel('lifetime (s)')
#plt.ylim(0,1e-7)
plt.savefig('lifetimes.png')
with open('lifetimes.txt', 'w') as newfile: #http://stackoverflow.com/questions/19302612/how-to-write-data-from-two-lists-into-columns-in-a-csv
    writer = csv.writer(newfile)
    writer.writerows(zip(powerlist, lifetimelist, chargelist))  
plt.figure(3)
plt.plot(powerlist,mobilitylist, 'bo')
plt.xlabel('laser power (uJ)')
plt.ylabel('mobility($cm^2V^{-1}s^{-1}$)')
plt.savefig('mobilities.png')
with open('mobilities.txt', 'w') as newfile:
    writer = csv.writer(newfile)
    writer.writerows(zip(powerlist, mobilitylist, chargelist))    
#with open('filename' + 'fitdata.csv', 'w') as outputfile:
#np.savetxt('overnight_open_fitdata.csv', datafromfit, delimiter=',')
    