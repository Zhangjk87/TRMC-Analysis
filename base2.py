import matplotlib as mpl

# hack; if run with default backend through command prompt get crash at end of program when making plots.
# PDF seems to work on my system. This should make no difference even if using Spyder because none of these
# plots are interactive.
mpl.use('pdf')
from functions import *
from mobility import *
from scipy.constants import *
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import t
import matplotlib.pyplot as plt
import scipy.signal
import math
import sys
import configparser
import csv


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
    return a * np.exp(-(x) / tau) + offset


def doubleExp(x, a, t1, b, t2, offset):
    return a * np.exp(-(x) / t1) + b * np.exp(-(x) / t2) + offset


def tripleExp(x, a, t1, b, t2, c, t3, offset):
    return a * np.exp(-(x) / t1) + b * np.exp(-(x) / t2) + c * np.exp(-(x) / t3) + offset


def readdata(folder, f):
    print(f)
    print(f[f.index('_') + 1:f.index('J')], 'J')  # current laser power
    pulseenergy = float(f[f.index('nm_') + 3:f.index('J')])
    p0 = f[f.index('p0_') + 3:f.index('V')]  # current laser power
    print('P0 = ' + p0 + 'V')
    p0val = float(p0)
    TRMCdata = np.genfromtxt(os.path.join(folder, f), delimiter=',', skip_header=4, usecols=(0, 1))
    # tempTRMCdata = np.genfromtxt(f[:f.index('_p0')]+ '_' ++'_p0_' + p0 + 'V.csv', delimiter=',', skip_header=4, usecols=(0,2))
    # TRMCdata[:,1] = np.add(TRMCdata[:,1], tempTRMCdata[:,1])
    # TRMCdata[:,1] = np.divide(TRMCdata[:,1], numberofarrays)
    # print('numberofarrays='+str(numberofarrays))
    return (TRMCdata, pulseenergy, p0val)


def findmaxormin(array):
    zerotime = np.argwhere(abs(array[:, 0]) - 2.5e-11 < 1e-17)  # 1e-16)
    maxlooktime = int(np.argwhere(abs(array[zerotime:, 0]) - 250e-9 > 0)[0])
    print(maxlooktime)
    if (abs(np.amax(array[zerotime:maxlooktime + zerotime, 1])) > abs(
            np.amin(array[zerotime:maxlooktime + zerotime, 1]))):
        return int(np.argmax(array[zerotime:maxlooktime + zerotime, 1]) + zerotime)
    elif (abs(np.amax(array[zerotime:maxlooktime + zerotime, 1])) < abs(
            np.amin(array[zerotime:maxlooktime + zerotime, 1]))):
        return int(np.argmin(array[zerotime:maxlooktime + zerotime, 1]) + zerotime)
    else:
        print("Error in findmaxormin(). Max not less than or greater than min array value.")
        raise RuntimeError('findmaxormin')
        # return None


def subtractOffset(array):
    zerotime = np.argwhere(abs(array[:, 0]) - 2.5e-11 < 1e-17)  # 1e-16)
    print(zerotime)
    datamean = np.mean(array[0:zerotime[0], 1])
    array[:, 1] = np.subtract(array[:, 1], datamean)
    return array


def saveArray(f, array):
    np.savetxt(f, array, delimiter=',')


def fitSingle(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(singleExp, xdata.T, ydata.T, p0=guess)

        print('a = ' + str(popt[0]))
        print('t1 = ' + str(popt[1]))
        print('offset = ' + str(popt[2]))

        # np.savetxt(f+'_fit.csv', popt, delimiter = ',')
        # np.savetxt(f+'_fit_cov.csv', popt, delimiter = ',')

        # from http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

        alpha = 0.05  # 95% confidence interval = 100*(1-alpha)

        n = len(ydata)  # number of data points
        p = len(popt)  # number of parameters

        dof = max(0, n - p)  # number of degrees of freedom

        # student-t value for the dof and confidence level
        tval = t.ppf(1.0 - alpha / 2., dof)
        params = []
        for i, p, var in zip(range(n), popt, np.diag(pcov)):
            sigma = var ** 0.5
            params.append('p{0}: {1} [{2}  {3}]'.format(i, p, p - sigma * tval, p + sigma * tval))
        # print(params)
        # with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        #    newfile.write('\n'.join(params))
        return (popt, pcov, params)
    except RuntimeError:
        print("Error - curve_fit failed")
        return (0, 0, 0)


def fitDouble(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(doubleExp, xdata.T, ydata.T, p0=guess)
        print('a = ' + str(popt[0]))
        print('t1 = ' + str(popt[1]))
        print('b = ' + str(popt[2]))
        print('t2 = ' + str(popt[3]))
        print('offset = ' + str(popt[4]))

        # np.savetxt(f+'_fit.csv', popt, delimiter = ',')
        # np.savetxt(f+'_fit_cov.csv', popt, delimiter = ',')

        # from http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

        alpha = 0.05  # 95% confidence interval = 100*(1-alpha)

        n = len(ydata)  # number of data points
        p = len(popt)  # number of parameters

        dof = max(0, n - p)  # number of degrees of freedom

        # student-t value for the dof and confidence level
        tval = t.ppf(1.0 - alpha / 2., dof)
        params = []
        for i, p, var in zip(range(n), popt, np.diag(pcov)):
            sigma = var ** 0.5
            params.append('p{0}: {1} [{2}  {3}]'.format(i, p, p - sigma * tval, p + sigma * tval))
        # print(params)
        # with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        #    newfile.write('\n'.join(params))
        return (popt, pcov, params)
    except RuntimeError:
        print("Error - curve_fit failed")
        return (0, 0, 0)


def fitTriple(xdata, ydata, guess):
    try:
        popt, pcov = curve_fit(tripleExp, xdata.T, ydata.T, p0=guess)
        print('a = ' + str(popt[0]))
        print('t1 = ' + str(popt[1]))
        print('b = ' + str(popt[2]))
        print('t2 = ' + str(popt[3]))
        print('c = ' + str(popt[4]))
        print('t3 = ' + str(popt[5]))
        print('offset = ' + str(popt[6]))

        # np.savetxt(f+'_fit.csv', popt, delimiter = ',')
        # np.savetxt(f+'_fit_cov.csv', popt, delimiter = ',')

        # from http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/

        alpha = 0.05  # 95% confidence interval = 100*(1-alpha)

        n = len(ydata)  # number of data points
        p = len(popt)  # number of parameters

        dof = max(0, n - p)  # number of degrees of freedom

        # student-t value for the dof and confidence level
        tval = t.ppf(1.0 - alpha / 2., dof)
        params = []
        for i, p, var in zip(range(n), popt, np.diag(pcov)):
            sigma = var ** 0.5
            params.append('p{0}: {1} [{2}  {3}]'.format(i, p, p - sigma * tval, p + sigma * tval))
        # print(params)
        # with open(f+'_fit_param_confidence_int.txt', 'w') as newfile:
        #    newfile.write('\n'.join(params))
        return (popt, pcov, params)
    except RuntimeError:
        print("Error - curve_fit failed")
        return (0, 0, 0)


def generateFitData(singleExponential, xdata, *popt):
    if singleExponential == 1:
        return (np.array([xdata, singleExp(xdata, *popt)]).T)
    elif singleExponential == 2:
        return (np.array([xdata, doubleExp(xdata, *popt)]).T)
    elif singleExponential == 3:
        return (np.array([xdata, tripleExp(xdata, *popt)]).T)
    else:
        raise RuntimeError('generateFitData')


def saveFitParams(f, params):
    with open(f + '_fit_param_confidence_int.txt', 'w') as newfile:
        newfile.write('\n'.join(params))


def deconvolve(array, cavityLifetime):
    # print('called')
    cavitydecayfft = np.fft.fft(np.exp(-array[:, 0] / cavityLifetime) / np.sum(np.exp(-array[:, 0] / cavityLifetime)))
    # datafft=np.fft.fft(signal.wiener(array[:,1],mysize=13))
    datafft = np.fft.fft(array[:, 1])
    fftconv = np.divide(datafft, cavitydecayfft)
    deconv = np.fft.ifft(fftconv)
    deconvolvedData = np.array([array[:, 0], deconv]).T
    return (deconvolvedData.real)


def binData(array, width):
    # Try to use averaging to make the data look nicer--http://stackoverflow.com/questions/21921178/binning-a-numpy-array
    ydata = array[:, 1]
    ydataavg = ydata[:(ydata.size // width) * width].reshape(-1, width).mean(axis=1)
    xdata = array[:, 0]
    xdataavg = xdata[:(xdata.size // width) * width].reshape(-1, width).mean(axis=1)
    return (np.array([xdataavg, ydataavg]).T)


def chargePerQD(I0, Fa, radius, packingFraction, thickness):
    # charges per QD = (Io*Fa*Volume_of_one_nanocrystal) / (packing_fraction*thickness)
    # packing_fraction = 0.6
    charge = (I0 * Fa * (4 / 3 * math.pi * radius ** 3) / (packingFraction * thickness))
    print('chargePerQD=' + str(charge))
    return (charge)


def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y

    
def mobility(dP, laserPower, folder, P0, sample_params_dict, Q, R0, responseTime, f0):

    config = configparser.RawConfigParser()
    config.read('config.ini')
    cavity_params_dict = dict(config.items('Cavity Params'))
    for keys in cavity_params_dict:
        cavity_params_dict[keys] = float(cavity_params_dict[keys])
    #correction for dP/P=ndV/V    
    dPCorrection = 1.42
    dP=dPCorrection*dP    
    
    pi=math.pi
    e0=epsilon_0#farads/meter
    q=e
    freq = c/(sample_params_dict['wavelength']) #Hz
    photonenergy = h*freq
    #print(photonenergy)#in J
    #J/photon, this is for 532nm right now
    #For now what this will do is just put dP in here, and it will calculate a mobility. Next step is build this into the fitting program for one click data analysis
    
    if dP<0:
        sign = 1
    else:
        sign = -1
        
    print('sign = ', sign)
    
    beta = cavity_params_dict['a']/cavity_params_dict['b']
    
    #dP=0.0014#volts
    #laserPower=    #J/pulse on power meter
    
    spotarea=pi*(sample_params_dict['spotdiameter']/2)**2#m^2
    #laserpower read in in joules
    #print(laserPower)
    I0=laserPower/photonenergy/spotarea
    print('I0 = '+str(I0))
    #print(format(I0, "e"))
    #careful with signs
    print('dP = ' + str(dP))
    K=sign*Q*(1/np.sqrt(R0)+sign*1)/(beta*e0*sample_params_dict['er']*f0*pi*cavity_params_dict['l']) #see absorption vs emission; I'm assume p0 is negative since that's what we get out of the detector, and a negative dP corresponds to an absorption with HP 462A amplifier
    print('K = ' + str(K))
    dG=sample_params_dict['illuminationfactor']*dP/P0/K
    #print(dG)
    
    phimu=dG/beta/I0/sample_params_dict['fa']/q/1e-4
    print('mu = ' + str(phimu))
    return(phimu,I0, dPCorrection*sample_params_dict['illuminationfactor']/P0/K/I0/sample_params_dict['fa']/q/beta/1e-4)        


def main(argv):

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 26}

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
    mpl.rcParams['axes.linewidth'] = 1.5  # set the value globally
    mpl.rcParams['lines.markersize'] = 12

    try:
        folder = sys.argv[1]
        print('folder = ', folder, '\n')
    except IndexError:
        print('no folder command line argument\n')
        folder = input('path containing TRMC data: ')

    # import sample params into dict
    config = configparser.RawConfigParser()
    config.read('config.ini')
    sample_params_dict = dict(config.items('Sample Params'))
    for keys in sample_params_dict:
        sample_params_dict[keys] = float(sample_params_dict[keys])
    sample_params_dict['lightReachingSample'] = 5.56*0.906*10**-1*sample_params_dict['od_magnitude']*0.92307525
    sample_params_dict['radius'] = 4.6e-9/2


    # import resonance params
    with open(os.path.join(folder, 'resonance_params.csv'), 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        next(reader, None)
        for row in reader:
            Q, R0, responseTime, f0 = list(map(float, row))

    # exponential=1

    # create lists that are needed to store results
    t1list = []
    t2list = []
    t3list = []
    a1list = []
    a2list = []
    a3list = []
    # ablist = []
    pulseenergylist = []
    photondensitylist = []
    mobilitylist = []
    mobilitydeconvlist = []
    chargelist = []

    savefolder = folder + '/analysisresults/'
    print(savefolder)
    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)

    saveArrays = False
    deconvolutionon = True
    filteron = True

    # main program block
    for file in os.listdir(folder):
        filenameconditions = 'AVG' not in file and 'decay' not in file and 'res' \
                             not in file and '.png' not in file and 'fileit' not in file \
                             and '.csv' in file and 'J_' in file and '_p0' in file
        if filenameconditions:
            baseFileName = savefolder + file[:file.index('_p0')]
            averagedData, pulseEnergy, P0 = readdata(folder, file)

            # averagedData = trim(averagedData, -9e9,3e-6)


            pulseEnergy = float(pulseEnergy) * sample_params_dict['lightReachingSample']
            pulseenergylist.append(pulseEnergy)
            freq = c / (sample_params_dict['wavelength'])  # Hz
            photonenergy = h * freq
            spotarea = pi * (sample_params_dict['spotdiameter'] / 2) ** 2  # m^2
            photonDensity = pulseEnergy / photonenergy / spotarea * sample_params_dict['fa'] / 100 ** 2  # photons/cm^2
            # also must output charge/QD, if necessary
            photondensitylist.append(photonDensity)
            averagedData = subtractOffset(averagedData)

            # save averaged data that has had its offset fixed
            # saveArray(baseFileName+'_combined.csv', averagedData)
            if filteron == True:
                # filter data
                cutoff = 15e6
                fs = 20e9
                order = 5

                filteredData = np.zeros(np.shape(averagedData))
                filteredData[:, 0] = np.copy(averagedData[:, 0])
                filteredData[:, 1] = butter_lowpass_filter(averagedData[:, 1], cutoff, fs, order)
            else:
                filteredData = np.copy(averagedData)
            mobilityDeconvData = np.copy(filteredData)  # #trim(filteredData, -9e9,3e-6)#

            # find max signal point. One may want to "bin" the data beforehand to reduce the influence of noise
            mobilityIndex = findmaxormin(filteredData)
            dP = filteredData[mobilityIndex, 1]
            print(mobilityIndex)
            print('dp=', dP)

            # mobility
            correctedP0 = P0 - sample_params_dict['p0offset']
            phimu, I0, PhiMuNormalization = mobility(dP, pulseEnergy, folder, correctedP0, sample_params_dict, Q, R0, responseTime, f0)
            mobilitylist.append(phimu)

            normalizedPhiMuArray = np.zeros(np.shape(filteredData))
            normalizedPhiMuArray[:, 0] = filteredData[:, 0]
            normalizedPhiMuArray[:, 1] = filteredData[:, 1] * PhiMuNormalization
            if saveArrays:
                saveArray(baseFileName + 'PhiMu_Normalized.csv', normalizedPhiMuArray)
            # this part of program does deconvolution
            print('arraylength=', np.shape(mobilityDeconvData))
            if deconvolutionon == True:
                deconvolvedData = deconvolve(mobilityDeconvData, responseTime)
            else:
                deconvolvedData = mobilityDeconvData
            print(deconvolvedData)
            # print(np.shape(deconvolvedData))
            deconvolvedDataBinned = deconvolvedData  # binData(deconvolvedData, 1)
            deconvolvedDataBinnedNormalized = np.copy(deconvolvedDataBinned)
            deconvolvedDataBinnedNormalized[:, 1] = deconvolvedDataBinnedNormalized[:, 1] * PhiMuNormalization
            if saveArrays:
                saveArray(baseFileName + '_combined_deconvolved.csv', deconvolvedDataBinnedNormalized)

            # plt.plot(deconvolvedData[:,0], deconvolvedData[:,1])
            #        plt.plot(deconvolvedDataBinned[:,0], deconvolvedDataBinned[:,1])
            #        plt.plot(averagedData[:,0], averagedData[:,1])
            #        plt.show()
            #        plt.close()

            mobilityDeconvIndex = findmaxormin(deconvolvedDataBinned[1:]) + 1
            print('mobilityDeconvIndex=' + str(mobilityDeconvIndex))
            dPDeconv = deconvolvedDataBinned[mobilityDeconvIndex, 1]
            print('dPDeconv = ' + str(dPDeconv))
            phimudeconv, ignore, ignore2 = mobility(dPDeconv, pulseEnergy, folder, correctedP0, sample_params_dict, Q, R0, responseTime, f0)
            mobilitydeconvlist.append(phimudeconv)
            charge = chargePerQD(I0, sample_params_dict['fa'],
                                 sample_params_dict['radius'],
                                 sample_params_dict['packingfraction'],
                                 sample_params_dict['thickness'])
            chargelist.append(charge)
            # chargelist.append(0)


            # to do: fit lifetime here
            exponential = sample_params_dict['exponential']
            fitdeconv = True
            if fitdeconv is True:
                if exponential == 1:
                    guess = [a, sample_params_dict['t1'], sample_params_dict['offset']]
                    popt, pcov, params = fitSingle(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 0],
                                                   deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 1], guess)
                elif exponential == 2:
                    guess = [sample_params_dict['a'], sample_params_dict['t1'], sample_params_dict['b'],
                             sample_params_dict['t2'], sample_params_dict['offset']]
                    popt, pcov, params = fitDouble(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 0],
                                                   deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 1], guess)
                    if popt is not 0:
                        t2list.append(popt[3])
                        a2list.append(popt[2])
                    else:
                        t2list.append(0)
                        a2list.append(0)
                        # ablist.append(popt[0]/popt[2])
                elif exponential == 3:
                    guess = [sample_params_dict['a'], sample_params_dict['t1'], sample_params_dict['b'],
                            sample_params_dict['t2'], sample_params_dict['c'], sample_params_dict['t3'],
                            sample_params_dict['offset']]
                    popt, pcov, params = fitTriple(deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 0],
                                                   deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 1], guess)
                    if popt is not 0:
                        if popt[3] > popt[5]:
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
                        # ablist.append(popt[0]/popt[2])
                if popt is not 0:
                    t1list.append(popt[1])
                    # def singleExp(x, a, tau, offset):
                    # return a*np.exp(-(x)/tau) + offset
                    print('fit mu = ',
                          singleExp(deconvolvedDataBinnedNormalized[mobilityDeconvIndex, 0], popt[0], popt[1], popt[2]))
                    a1list.append(
                        singleExp(deconvolvedDataBinnedNormalized[mobilityDeconvIndex, 0], popt[0], popt[1], popt[2]))
                    # a1list.append(popt[0])
                    fitArray = generateFitData(exponential, deconvolvedDataBinnedNormalized[mobilityDeconvIndex:, 0], *popt)
                    if saveArrays:
                        saveArray(folder + baseFileName + '_lifetimeFit.csv', fitArray)
                    saveFitParams(baseFileName + '_fitParams.txt', params)
                else:
                    t1list.append(0)
                    a1list.append(0)
                    # save fit data and fit params

            if fitdeconv is False:
                if exponential == 1:
                    guess = [a, sample_params_dict['t1'], sample_params_dict['offset']]
                    popt, pcov, params = fitSingle(filteredData[mobilityIndex:, 0], filteredData[mobilityIndex:, 1], guess)
                elif exponential == 2:
                    guess = [sample_params_dict['a'], sample_params_dict['t1'], sample_params_dict['b'],
                             sample_params_dict['t2'], sample_params_dict['offset']]
                    popt, pcov, params = fitDouble(filteredData[mobilityIndex:, 0], filteredData[mobilityIndex:, 1], guess)
                    t2list.append(popt[3])
                    a2list.append(popt[2])
                    # ablist.append(popt[0]/popt[2])
                elif exponential == 3:
                    guess = [sample_params_dict['a'], sample_params_dict['t1'], sample_params_dict['b'],
                            sample_params_dict['t2'], sample_params_dict['c'], sample_params_dict['t3'],
                            sample_params_dict['offset']]
                    popt, pcov, params = fitTriple(filteredData[filteredData:, 0], filteredData[mobilityIndex:, 1], guess)

                    if popt[3] > popt[5]:
                        t3list.append(popt[3])
                        a3list.append(popt[2])
                        t2list.append(popt[5])
                        a2list.append(popt[4])
                    else:
                        t2list.append(popt[3])
                        a2list.append(popt[2])
                        t3list.append(popt[5])
                        a3list.append(popt[4])
                        # ablist.append(popt[0]/popt[2])
                if popt is not 0:
                    t1list.append(popt[1])
                    # def singleExp(x, a, tau, offset):
                    # return a*np.exp(-(x)/tau) + offset
                    print('fit mu = ', singleExp(filteredData[mobilityIndex, 0], popt[0], popt[1], popt[2]))
                    a1list.append(singleExp(filteredData[mobilityIndex, 0], popt[0], popt[1], popt[2]))
                    # a1list.append(popt[0])
                    fitArray = generateFitData(exponential, filteredData[mobilityIndex:, 0], *popt)
                    if saveArrays:
                        saveArray(baseFileName + '_lifetimeFit.csv', fitArray)
                    saveFitParams(baseFileName + '_fitParams.txt', params)
                else:
                    t1list.append(0)
                    # save fit data and fit params

            print('donefit')

            # make plots
            plt.figure(1)
            plt.plot(filteredData[:, 0] / 1e-9, -filteredData[:, 1] / P0, label=str(format(charge, '.2f')))
            plt.figure(2)
            plt.plot(deconvolvedDataBinned[1:, 0] / 1e-9, -deconvolvedDataBinned[1:, 1] / P0,
                     label=str(format(pulseEnergy / 1e-6, '.0f') + ' $\mathrm{\mu J}$'))

            plt.figure(3)
            plt.plot(deconvolvedDataBinned[1:, 0] / 1e-9, deconvolvedDataBinned[1:, 1] * PhiMuNormalization, 'g-')
            plt.plot(averagedData[:, 0] / 1e-9, averagedData[:, 1] * PhiMuNormalization, '-b')
            plt.plot(filteredData[:, 0] / 1e-9, filteredData[:, 1] * PhiMuNormalization, '-k')
            if popt is not 0:
                plt.plot(fitArray[:, 0] / 1e-9, fitArray[:, 1], 'r-')
            plt.xlabel('Time (ns)')
            plt.ylabel('$\mathrm{\Delta}$P (V)')
            plt.xlim(1, filteredData[-1, 0] / 1e-9)
            plt.xscale('log')
            plt.savefig(baseFileName + '_fit.png')
            plt.close()

            plt.figure(99)
            plt.plot(deconvolvedDataBinned[1:, 0] / 1e-9, deconvolvedDataBinned[1:, 1] * PhiMuNormalization, 'g-')
            plt.plot(averagedData[:, 0] / 1e-9, averagedData[:, 1] * PhiMuNormalization, '-b')
            plt.plot(filteredData[:, 0] / 1e-9, filteredData[:, 1] * PhiMuNormalization, '-k')
            if popt is not 0:
                plt.plot(fitArray[:, 0] / 1e-9, fitArray[:, 1], 'r-')
            plt.xlabel('Time (ns)')
            plt.ylabel('$\mathrm{\Delta}$P (V)')
            plt.xlim(deconvolvedDataBinned[0, 0] / 1e-9, deconvolvedDataBinned[-1, 0] / 1e-9)
            plt.savefig(baseFileName + '_fit_deconvolved.png')
            plt.close()
            print()
            # write stuff to file

    olddata = False
    exponential = sample_params_dict['exponential']
    if olddata:
        if exponential == 1:
            summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list]).T
        elif exponential == 2:
            summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list]).T
        elif exponential == 3:
            summary = np.array([pulseenergylist, chargelist, mobilitylist, mobilitydeconvlist, t1list, t2list, t3list]).T
    else:
        if exponential == 1:
            summary = np.array(
                [pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list]).T
        elif exponential == 2:
            summary = np.array(
                [pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list, a2list,
                 t2list]).T
        elif exponential == 3:
            summary = np.array(
                [pulseenergylist, photondensitylist, chargelist, mobilitylist, mobilitydeconvlist, a1list, t1list, a2list,
                 t2list, a3list, t3list]).T

    print('writing summary')
    saveArray(savefolder + 'summary.csv', summary)
    print('savedarray')
    pulseenergylistuJ = np.array(pulseenergylist) / 1e-6

    # make more plots
    plt.figure(1)
    plt.xlabel('Time (ns)')
    plt.ylabel('$\mathrm{\Delta P/P_0}$')
    plt.xlim(left=0)
    plt.legend()
    plt.savefig(savefolder + 'decays.png')
    plt.close()

    plt.figure(2)
    plt.xlabel('Time (ns)')
    plt.ylabel('$\mathrm{\Delta P/P_0}$')
    plt.xlim(left=0)
    plt.legend()
    plt.savefig(savefolder + 'deconvolved_decays.png')
    plt.close()

    plt.figure(4)
    plt.plot(chargelist, mobilitylist, 'o')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel('$\phi \Sigma \mu$ ($\mathrm{cm^2V^{-1}s^{-1}}$)')
    plt.savefig(savefolder + 'mobilities.png')
    plt.close()

    plt.figure(5)
    plt.plot(chargelist, mobilitydeconvlist, 'o')
    plt.xlabel('Charge per Quantum Dot')
    plt.ylabel('$\phi \Sigma \mu$ ($\mathrm{cm^2V^{-1}s^{-1}}$)')
    plt.savefig(savefolder + 'deconvolved_mobilities.png')
    plt.close()

    if exponential == 1:
        plt.figure(6, figsize=(8, 7))
        plt.plot(chargelist, np.array(t1list) / 1e-9, 'o')
        plt.xlabel('Charge per Quantum Dot')
        plt.ylabel('Lifetime (ns)')
        plt.savefig(savefolder + 'lifetimes.png')
        plt.close()
    elif exponential == 2:
        plt.figure(6, figsize=(8, 7))
        plt.plot(chargelist, np.array(t1list) / 1e-9, 'o', label=r'$\tau_1$')
        plt.plot(chargelist, np.array(t2list) / 1e-9, 'o', label=r'$\tau_2$')
        plt.xlabel('Charge per Quantum Dot')
        plt.ylabel('Lifetime (ns)')
        plt.legend()
        plt.savefig(savefolder + 'lifetimes.png')
        plt.close()
        # plt.figure(7)
        # plt.plot(chargelist, ablist, 'o')
        # plt.xlabel('Charge per Quantum Dot')
        # plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
        # plt.savefig(savefolder+'timeconstantratios.png')
        # plt.close()
    elif exponential == 3:
        plt.figure(6, figsize=(8, 7))
        plt.plot(chargelist, np.array(t1list) / 1e-9, 'o', label=r'$\tau_1$')
        plt.plot(chargelist, np.array(t2list) / 1e-9, 'o', label=r'$\tau_2$')
        plt.plot(chargelist, np.array(t3list) / 1e-9, 'o', label=r'$\tau_3$')
        plt.xlabel('Charge per Quantum Dot')
        plt.ylabel('Lifetime (ns)')
        plt.yscale('log', nonposy='clip')
        # plt.xscale('log', nonposy='clip')
        # plt.legend()
        plt.savefig(savefolder + 'lifetimes.png')
        plt.close()
        # plt.figure(7)
        # plt.plot(chargelist, ablist, 'o')
        # plt.xlabel('Charge per Quantum Dot')
        # plt.ylabel(r'Ratio of fit amplitudes for $\tau_1/\tau_2$')
        # plt.savefig(savefolder+'timeconstantratios.png')
        # plt.close()

    print('done')


if __name__ == '__main__':
    main(sys.argv[1:])

# to do: deconvolution and fit mobility again, and plot both datas

