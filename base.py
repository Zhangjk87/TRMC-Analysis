# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:38:40 2015

@author: Daniel
"""

from cavityparams import *
import matplotlib as mpl
#hack; if run with default backend through command prompt get crash at end of program when making plots. PDF seems to work on my system. This should make no difference even if using Spyder because none of these plots are interactive.
mpl.use('pdf')
from functions import *
from mobility import *
from scipy.constants import *
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
import matplotlib.pyplot as plt
import math
import sys
#to make fonts from plots look normal

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

try:
    folder=sys.argv[1]
    print('folder = ', folder, '\n')
except:
    print('no folder command line argument\n')
    folder = input('path containing TRMC data: ')

sys.path.insert(0, folder)
from sampleparams import *
from resonanceparams import *

#exponential=1

#create lists that are needed to store results
t1list = []
t2list = []
t3list = []
a1list = []
a2list = []
a3list = []
#ablist = []
pulseenergylist = []
photondensitylist = []
mobilitylist = []
mobilitydeconvlist = []
chargelist = []

os.chdir(folder)

savefolder=os.path.join('./analysisresults/')

saveArrays = False
deconvolutionon=True
filteron=True


if not os.path.exists(savefolder):
    os.makedirs(savefolder)

#main program block
for f in os.listdir("."):
    filenameconditions = 'AVG' not in f and 'decay' not in f and 'res' not in f and '.png' not in f and 'fit' not in f and '.csv' in f and 'J_' in f and '_p0' in f
    if filenameconditions:
        baseFileName = savefolder + f[:f.index('_p0')]
        print(f)
        averagedData, pulseEnergy, P0 = readdata(f)
        

        #averagedData = trim(averagedData, -9e9,3e-6)


        pulseEnergy = float(pulseEnergy)*lightReachingSample
        pulseenergylist.append(pulseEnergy)
        freq = c/(wavelength) #Hz
        photonenergy = h*freq
        spotarea=pi*(spotdiameter/2)**2#m^2
        photonDensity=pulseEnergy/photonenergy/spotarea*Fa/100**2 #photons/cm^2
        #also must output charge/QD, if necessary        
        photondensitylist.append(photonDensity)
        averagedData = subtractOffset(averagedData)
        
        #save averaged data that has had its offset fixed
        #saveArray(baseFileName+'_combined.csv', averagedData)        
        if filteron==True:
        #filter data                
            cutoff=15e6
            fs=20e9
            order=5
            
            filteredData=np.zeros(np.shape(averagedData))
            filteredData[:,0] = np.copy(averagedData[:,0])
            filteredData[:,1] = butter_lowpass_filter(averagedData[:,1], cutoff, fs, order)
        else:
            filteredData =np.copy(averagedData)
        mobilityDeconvData =np.copy(filteredData)# #trim(filteredData, -9e9,3e-6)#
        
        #find max signal point. One may want to "bin" the data beforehand to reduce the influence of noise
        mobilityIndex = findmaxormin(filteredData)
        dP = filteredData[mobilityIndex, 1]
        print(mobilityIndex)        
        print('dp=', dP)
        

        
        #mobility
        correctedP0=P0-P0offset
        phimu, I0, PhiMuNormalization = mobility(dP, pulseEnergy,folder,correctedP0)
        mobilitylist.append(phimu)       
        
        normalizedPhiMuArray=np.zeros(np.shape(filteredData))
        normalizedPhiMuArray[:,0]=filteredData[:,0]
        normalizedPhiMuArray[:,1]=filteredData[:,1]*PhiMuNormalization
        if saveArrays:
            saveArray(baseFileName+'PhiMu_Normalized.csv', normalizedPhiMuArray)
        #this part of program does deconvolution
        print('arraylength=', np.shape(mobilityDeconvData))
        if deconvolutionon==True:
            deconvolvedData = deconvolve(mobilityDeconvData, responseTime)
        else:
            deconvolvedData=mobilityDeconvData
        print(deconvolvedData)
        #print(np.shape(deconvolvedData))
        deconvolvedDataBinned = deconvolvedData#binData(deconvolvedData, 1)
        deconvolvedDataBinnedNormalized = np.copy(deconvolvedDataBinned)
        deconvolvedDataBinnedNormalized[:,1] = deconvolvedDataBinnedNormalized[:,1] * PhiMuNormalization
        if saveArrays:
            saveArray(baseFileName+'_combined_deconvolved.csv', deconvolvedDataBinnedNormalized)
        
    
#        plt.plot(deconvolvedData[:,0], deconvolvedData[:,1])
#        plt.plot(deconvolvedDataBinned[:,0], deconvolvedDataBinned[:,1])
#        plt.plot(averagedData[:,0], averagedData[:,1])
#        plt.show()
#        plt.close()
        
        mobilityDeconvIndex = findmaxormin(deconvolvedDataBinned[1:])+1
        print('mobilityDeconvIndex='+str(mobilityDeconvIndex))
        dPDeconv = deconvolvedDataBinned[mobilityDeconvIndex,1]
        print('dPDeconv = '+str(dPDeconv))
        phimudeconv, ignore, ignore2 = mobility(dPDeconv, pulseEnergy,folder,correctedP0)
        mobilitydeconvlist.append(phimudeconv)
        charge=chargePerQD(I0, Fa, radius, packingFraction, thickness)
        chargelist.append(charge)
        #chargelist.append(0)
       
       
                    #to do: fit lifetime here
        fitdeconv=True
        if fitdeconv is True:
            if exponential == 1:
                guess = [a, t1, offset]
                popt, pcov, params = fitSingle(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,0],deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,1], guess)      
            elif exponential == 2:
                guess = [a, t1, b, t2, offset]
                popt, pcov, params = fitDouble(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,0],deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,1], guess)  
                if popt is not 0:
                    t2list.append(popt[3])
                    a2list.append(popt[2])
                else:
                    t2list.append(0)
                    a2list.append(0)
                #ablist.append(popt[0]/popt[2])
            elif exponential == 3:
                guess = [a, t1, b, t2, c, t3, offset]
                popt, pcov, params = fitTriple(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,0],deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,1], guess)  
                if popt is not 0:
                    if popt[3]>popt[5]:
                        t3list.append(popt[3])
                        a3list.append(popt[2])
                        t2list.append(popt[5])
                        a2list.append(popt[4])
                    else:
                        t2list.append(popt[3])
                        a2list.append(popt[2])
                        t3list.append(popt[5])
                        a3list.append(popt[4])
                else:
                        t3list.append(0)
                        a3list.append(0)
                        t2list.append(0)
                        a2list.append(0)
                #ablist.append(popt[0]/popt[2])
            if popt is not 0:
                t1list.append(popt[1])
                #def singleExp(x, a, tau, offset):
                    #return a*np.exp(-(x)/tau) + offset
                print('fit mu = ', singleExp(deconvolvedDataBinnedNormalized[mobilityDeconvIndex,0], popt[0], popt[1], popt[2]))
                a1list.append(singleExp(deconvolvedDataBinnedNormalized[mobilityDeconvIndex,0], popt[0], popt[1], popt[2]))
                #a1list.append(popt[0])
                fitArray = generateFitData(exponential, deconvolvedDataBinnedNormalized[mobilityDeconvIndex:,0], *popt)
                if saveArrays:
                    saveArray(baseFileName+'_lifetimeFit.csv', fitArray)
                saveFitParams(baseFileName+'_fitParams.txt', params)   
            else:
                t1list.append(0)
                a1list.append(0)
                    #save fit data and fit params 
    
        if fitdeconv is False:
            if exponential == 1:
                guess = [a, t1, offset]
                popt, pcov, params = fitSingle(filteredData[mobilityIndex:,0],filteredData[mobilityIndex:,1], guess)      
            elif exponential == 2:
                guess = [a, t1, b, t2, offset]
                popt, pcov, params = fitDouble(filteredData[mobilityIndex:,0],filteredData[mobilityIndex:,1], guess)
                t2list.append(popt[3])
                a2list.append(popt[2])
                #ablist.append(popt[0]/popt[2])
            elif exponential == 3:
                guess = [a, t1, b, t2, c, t3, offset]
                popt, pcov, params = fitTriple(filteredData[filteredData:,0],filteredData[mobilityIndex:,1], guess)
                
                if popt[3]>popt[5]:
                    t3list.append(popt[3])
                    a3list.append(popt[2])
                    t2list.append(popt[5])
                    a2list.append(popt[4])
                else:
                    t2list.append(popt[3])
                    a2list.append(popt[2])
                    t3list.append(popt[5])
                    a3list.append(popt[4])
                #ablist.append(popt[0]/popt[2])
            if popt is not 0:
                t1list.append(popt[1])
                #def singleExp(x, a, tau, offset):
                    #return a*np.exp(-(x)/tau) + offset
                print('fit mu = ', singleExp(filteredData[mobilityIndex,0], popt[0], popt[1], popt[2]))
                a1list.append(singleExp(filteredData[mobilityIndex,0], popt[0], popt[1], popt[2]))
                #a1list.append(popt[0])
                fitArray = generateFitData(exponential, filteredData[mobilityIndex:,0], *popt)
                if saveArrays:
                    saveArray(baseFileName+'_lifetimeFit.csv', fitArray)
                saveFitParams(baseFileName+'_fitParams.txt', params)   
            else:
                t1list.append(0)
                    #save fit data and fit params            

        print('donefit')
       
       #make plots
        plt.figure(1)
        plt.plot(filteredData[:,0]/1e-9, -filteredData[:,1]/P0, label=str(format(charge, '.2f')))
        plt.figure(2)
        plt.plot(deconvolvedDataBinned[1:,0]/1e-9, -deconvolvedDataBinned[1:,1]/P0,label=str(format(pulseEnergy/1e-6, '.0f') + ' $\mathrm{\mu J}$'))
        
        plt.figure(3)
        plt.plot(deconvolvedDataBinned[1:,0]/1e-9, deconvolvedDataBinned[1:,1]*PhiMuNormalization, 'g-')
        plt.plot(averagedData[:,0]/1e-9, averagedData[:,1]*PhiMuNormalization, '-b')
        plt.plot(filteredData[:,0]/1e-9, filteredData[:,1]*PhiMuNormalization, '-k')
        if popt is not 0:
            plt.plot(fitArray[:,0]/1e-9, fitArray[:,1],'r-')
        plt.xlabel('Time (ns)')
        plt.ylabel('$\mathrm{\Delta}$P (V)')
        plt.xlim(1, filteredData[-1,0]/1e-9)
        plt.xscale('log')
        plt.savefig(baseFileName+'_fit.png')
        plt.close()
    
        plt.figure(99)
        plt.plot(deconvolvedDataBinned[1:,0]/1e-9, deconvolvedDataBinned[1:,1]*PhiMuNormalization, 'g-')
        plt.plot(averagedData[:,0]/1e-9, averagedData[:,1]*PhiMuNormalization, '-b')
        plt.plot(filteredData[:,0]/1e-9, filteredData[:,1]*PhiMuNormalization, '-k')
        if popt is not 0:
            plt.plot(fitArray[:,0]/1e-9, fitArray[:,1],'r-')
        plt.xlabel('Time (ns)')
        plt.ylabel('$\mathrm{\Delta}$P (V)')
        plt.xlim(deconvolvedDataBinned[0,0]/1e-9, deconvolvedDataBinned[-1,0]/1e-9)
        plt.savefig(baseFileName+'_fit_deconvolved.png')
        plt.close()
        print()
        #write stuff to file

      

olddata=False        
if olddata:
    if exponential == 1:
        summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list]).T
    elif exponential == 2:
        summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list]).T
    elif exponential == 3:
        summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list, t3list]).T
else:
    if exponential == 1:
        summary = np.array([pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list]).T
    elif exponential == 2:
        summary = np.array([pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list, a2list, t2list]).T
    elif exponential == 3:
        summary = np.array([pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list, a2list, t2list, a3list, t3list]).T
        
print('writing summary')  
saveArray(savefolder + 'summary.csv', summary)
print('savedarray')
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
    #plt.figure(7)
    #plt.plot(chargelist, ablist, 'o')
    #plt.xlabel('Charge per Quantum Dot')
    #plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
    #plt.savefig(savefolder+'timeconstantratios.png')
    #plt.close()
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
    #plt.figure(7)
    #plt.plot(chargelist, ablist, 'o')
    #plt.xlabel('Charge per Quantum Dot')
    #plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
    #plt.savefig(savefolder+'timeconstantratios.png')
    #plt.close()

print('done')
quit()
        
        #to do: deconvolution and fit mobility again, and plot both datas
    
