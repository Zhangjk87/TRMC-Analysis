# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:38:40 2015

@author: Daniel
"""

from cavityparams import *

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

os.chdir(folder)

#main program block
for f in os.listdir("."):
    filenameconditions = 'AVG' not in f and 'decay' not in f and 'p0' not in f and 'res' not in f and '.png' not in f and 'fit' not in f and '.csv' in f and 'J_' in f and '_1.csv' in f
    if filenameconditions:
        averagedData, pulseEnergy= readdata(f, numberoffiles)
        
        
        #to do: fit lifetime here        
        
        averagedData=subtractOffset(averagedData)
        saveArray(f[:f.index('_1.csv')+'averaged.csv', averagedData)        
        
        mobilityIndex = findmaxormin(averagedData)
        print(mobilityIndex)
        
        #to do: fit mobility here.
        
        #to do: deconvolution and fit mobility again, and plot both datas
    