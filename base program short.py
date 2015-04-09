# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:38:40 2015

@author: Daniel
"""

from cavityparams import *

from functions_short import *
from mobility import *
from scipy.constants import *
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
import matplotlib.pyplot as plt
import math

#to make fonts from plots look normal
import matplotlib as mpl
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
#matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.rc('font', **{'sans-serif' : 'Arial', 'family' : 'sans-serif'})
#matplotlib.rcparams()
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 26}

mpl.rc('font', **font)
mpl.rcParams.update({'figure.autolayout': True})
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally
mpl.rcParams['lines.markersize'] = 12  


folder = input('path containing TRMC data: ')

import sys
sys.path.insert(0, folder)
from sampleparams import *
from resonanceparams import *

exponential=1

#create lists that are needed to store results
t1list = []
t2list = []
t3list = []
ablist = []
pulseenergylist = []
mobilitylist = []
mobilitydeconvlist = []
chargelist = []

os.chdir(folder)

savefolder=os.path.join('./analysisresults/')

if not os.path.exists(savefolder):
    os.makedirs(savefolder)

#main program block
for f in os.listdir("."):
    filenameconditions = 'AVG' not in f and 'decay' not in f and 'res' not in f and '.png' not in f and 'fit' not in f and '.csv' in f and 'J_' in f and '_1_p0' in f and 'long' not in f
    if filenameconditions:
        baseFileName = savefolder + f[:f.index('_1_p0')]
        
        averagedData, pulseEnergy, P0 = readdata(f)

        pulseEnergy = float(pulseEnergy)*lightReachingSample
        pulseenergylist.append(pulseEnergy)
        #also must output charge/QD, if necessary        
        
        averagedData = subtractOffset(averagedData)
        
        #save averaged data that has had its offset fixed
        saveArray(baseFileName+'_combined.csv', averagedData)        
        
        #find max signal point. One may want to "bin" the data beforehand to reduce the influence of noise
        mobilityIndex = findmaxormin(averagedData)
        dP = averagedData[mobilityIndex, 1]
        print(mobilityIndex)        
        
                #to do: fit lifetime here
        if exponential == 1:
            guess = [a, t1, offset]
            popt, pcov, params = fitSingle(averagedData[mobilityIndex:,0],averagedData[mobilityIndex:,1], guess)      
        elif exponential == 2:
            guess = [a, t1, b, t2, offset]
            popt, pcov, params = fitDouble(averagedData[mobilityIndex:,0],averagedData[mobilityIndex:,1], guess)  
            t2list.append(popt[3])
            ablist.append(popt[0]/popt[2])
        elif exponential == 3:
            guess = [a, t1, b, t2, c, t3, offset]
            popt, pcov, params = fitTriple(averagedData[mobilityIndex:,0],averagedData[mobilityIndex:,1], guess)  
            
            if popt[3]>popt[5]:
                t3list.append(popt[3])
                t2list.append(popt[5])
            else:
                t2list.append(popt[3])
                t3list.append(popt[5])
            ablist.append(popt[0]/popt[2])
        t1list.append(popt[1])
        #save fit data and fit params 
        fitArray = generateFitData(exponential, averagedData[mobilityIndex:,0], *popt)
        saveArray(baseFileName+'_lifetimeFit.csv', fitArray)
        saveFitParams(baseFileName+'_fitParams.txt', params)
        
        #mobility
        correctedP0=P0-P0offset
        phimu, I0 = mobility(dP, pulseEnergy,folder,correctedP0)
        mobilitylist.append(phimu)       


        #this part of program does deconvolution
        deconvolvedData = deconvolve(averagedData, responseTime)
        
        deconvolvedDataBinned = binData(deconvolvedData, 50)  
        saveArray(baseFileName+'_combined_deconvolved.csv', deconvolvedDataBinned)
        
#        plt.plot(deconvolvedData[:,0], deconvolvedData[:,1])
#        plt.plot(deconvolvedDataBinned[:,0], deconvolvedDataBinned[:,1])
#        plt.plot(averagedData[:,0], averagedData[:,1])
#        plt.show()
#        plt.close()
        
        mobilityDeconvIndex = findmaxormin(deconvolvedDataBinned[1:])+1
        print('mobilityDeconvIndex='+str(mobilityDeconvIndex))
        dPDeconv = deconvolvedDataBinned[mobilityDeconvIndex,1]
        print('dPDeconv = '+str(dPDeconv))
        phimudeconv, ignore = mobility(dPDeconv, pulseEnergy,folder,correctedP0)
        mobilitydeconvlist.append(phimudeconv)
        charge=chargePerQD(I0, Fa, radius, packingFraction, thickness)
        chargelist.append(charge)
        #chargelist.append(0)
        #make plots
        plt.figure(1)
        plt.plot(averagedData[:,0]/1e-9, -averagedData[:,1]/P0, label=str(format(charge, '.2f')))
        plt.figure(2)
        plt.plot(deconvolvedDataBinned[:,0]/1e-9, -deconvolvedDataBinned[:,1]/P0,label=str(format(pulseEnergy/1e-6, '.0f') + ' $\mathrm{\mu J}$'))
        
        plt.figure(3)
        plt.plot(averagedData[:,0]/1e-9, averagedData[:,1])
        plt.plot(fitArray[:,0]/1e-9, fitArray[:,1],'r-')
        plt.xlabel('Time (ns)')
        plt.ylabel('$\mathrm{\Delta}$P (V)')
        plt.savefig(baseFileName+'_fit.png')
        plt.close()
        print()
        #write stuff to file
if exponential == 1:
    summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list]).T
elif exponential == 2:
    summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list, ablist]).T
elif exponential == 3:
    summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list, t3list, ablist]).T

saveArray(savefolder + 'summary.csv', summary)

pulseenergylistuJ=np.array(pulseenergylist)/1e-6

#make more plots
plt.figure(1)
plt.xlabel('Time (ns)')
plt.ylabel('$\mathrm{\Delta P/P_0}$')
plt.xlim(left=0)
plt.legend()
plt.savefig(savefolder+'decays.png')
plt.close()

plt.figure(2)
plt.xlabel('Time (ns)')
plt.ylabel('$\mathrm{\Delta P/P_0}$')
plt.xlim(left=0)
plt.legend()
plt.savefig(savefolder+'deconvolved_decays.png')
plt.close()

plt.figure(4)
plt.plot(chargelist, mobilitylist, 'o')
plt.xlabel('Charge per Quantum Dot')
plt.ylabel('$\phi \Sigma \mu$ ($\mathrm{cm^2V^{-1}s^{-1}}$)')
plt.savefig(savefolder+'mobilities.png')
plt.close()

plt.figure(5)
plt.plot(chargelist, mobilitydeconvlist, 'o')
plt.xlabel('Charge per Quantum Dot')
plt.ylabel('$\phi \Sigma \mu$ ($\mathrm{cm^2V^{-1}s^{-1}}$)')
plt.savefig(savefolder+'deconvolved_mobilities.png')
plt.close()

    
if exponential == 1:
    plt.figure(6, figsize=(8,7))
    plt.plot(chargelist, np.array(t1list)/1e-9, 'o')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel('Lifetime (ns)')
    plt.savefig(savefolder+'lifetimes.png')
    plt.close()
elif exponential == 2:
    plt.figure(6,figsize=(8,7))
    plt.plot(chargelist, np.array(t1list)/1e-9, 'o', label=r'$\tau_1$')
    plt.plot(chargelist, np.array(t2list)/1e-9, 'o', label=r'$\tau_2$')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel('Lifetime (ns)')
    plt.legend()
    plt.savefig(savefolder+'lifetimes.png')
    plt.close()
    plt.figure(7)
    plt.plot(chargelist, ablist, 'o')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
    plt.savefig(savefolder+'timeconstantratios.png')
    plt.close()
elif exponential == 3:
    plt.figure(6, figsize=(8,7))
    plt.plot(chargelist, np.array(t1list)/1e-9, 'o', label=r'$\tau_1$')
    plt.plot(chargelist, np.array(t2list)/1e-9, 'o', label=r'$\tau_2$')
    plt.plot(chargelist, np.array(t3list)/1e-9, 'o', label=r'$\tau_3$')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel('Lifetime (ns)')
    plt.yscale('log', nonposy='clip')
    #plt.xscale('log', nonposy='clip')
    #plt.legend()
    plt.savefig(savefolder+'lifetimes.png')
    plt.close()
    plt.figure(7)
    plt.plot(chargelist, ablist, 'o')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
    plt.savefig(savefolder+'timeconstantratios.png')
    plt.close()

        
        #to do: deconvolution and fit mobility again, and plot both datas
    