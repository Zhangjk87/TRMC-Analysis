# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:38:40 2015

@author: Daniel
"""

from cavityparams import *
from tempdata.sampleparams import *

from functions import *

from scipy.constants import *
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
import matplotlib.pyplot as plt
import math

from tempdata.resonanceparams import *

#create lists that are needed to store results
lifetimelist=[]
powerlist = []
mobilitylist = []
chargelist = []

#important parameters for the program
numberoffiles=2

folder = './tempdata/'

savefolder = './analysisresults/'

os.chdir(folder)

if not os.path.exists(savefolder):
    os.makedirs(savefolder)

#main program block
for f in os.listdir("."):
    filenameconditions = 'AVG' not in f and 'decay' not in f and 'p0' not in f and 'res' not in f and '.png' not in f and 'fit' not in f and '.csv' in f and 'J_' in f and '_1.csv' in f
    if filenameconditions:
        averagedData, pulseEnergy= readdata(f, numberoffiles)
        
        averagedData=subtractOffset(averagedData)
        
        #save averaged data that has had its offset fixed
        saveArray(savefolder + f[:f.index('_1.csv')]+'_combined.csv', averagedData)        
        
                #to do: fit lifetime here
        if singleExponential==True:
            guess=[a, t1, offset]
            popt, pcov, params = fitSingle(averagedData[:,0],averagedData[:,1], guess)      
        if singleExponential==False:
            guess=[a, t1, b, t2, offset]
            popt, pcov, params = fitSingle(averagedData[:,0],averagedData[:,1], guess)  
            
        #save fit data and fit params
        
        #find max signal point. One may want to "bin" the data beforehand to reduce the influence of noise
       
        mobilityIndex = findmaxormin(averagedData)
        print(mobilityIndex)
        
        
        
        #to do: fit mobility here.
        
        #to do: deconvolution and fit mobility again, and plot both datas
    