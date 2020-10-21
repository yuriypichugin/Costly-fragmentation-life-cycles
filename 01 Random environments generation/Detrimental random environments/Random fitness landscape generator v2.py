# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 07:25:42 2017

@author: pichugin
"""
import numpy as np
#import pandas as pd
import LifeCycleSupplementary_v1_3_3 as LC
import matplotlib.pyplot as plt
from copy import deepcopy as deepcopy

NMax = 20
repeats = 100
Files = 100
for i in np.arange(Files):
    
#    NameFileB = 'B_profiles_Random_'+str(i)+'.txt'
#    NameFileD = 'D_profiles_Random_'+str(i)+'.txt'
#    ModelB = 3
#    ModelD = 3
#    Params = []
    
    NameFileB = 'B_profiles_Detrimental_'+str(i)+'.txt'
    NameFileD = 'D_profiles_Detrimental_'+str(i)+'.txt'
    ModelB = 7
    ModelD = 6
    Params = []
    
#    NameFileB = 'B_profiles_Beneficial_'+str(i)+'.txt'
#    NameFileD = 'D_profiles_Beneficial_'+str(i)+'.txt'
#    ModelB = 6
#    ModelD = 7
#    Params = []

    
#    NameFileB = 'B_profiles_Unimodal_10_'+str(i)+'.txt'
#    NameFileD = 'D_profiles_Unimodal_10_'+str(i)+'.txt'
#    ModelB = 8
#    ModelD = 8
#    Params = [10]
    
    BSet = np.zeros((repeats, NMax+1))
    DSet = np.zeros((repeats, NMax+1))
    
    for sample in np.arange(repeats):
        ID = sample + 10000*i        
        
#        BProfile = LC.BirthInit(3, NMax, [])
#        DProfile = LC.DeathInit(3, NMax, [])
        BProfile = LC.BirthInit(ModelB, NMax, Params)
        DProfile = LC.DeathInit(ModelD, NMax, Params)
        
        BProfile.insert(0, ID)       
        DProfile.insert(0, ID)
        
        BSet[sample,:] = deepcopy(BProfile)
        DSet[sample,:] = deepcopy(DProfile)
    
    np.savetxt(NameFileB, BSet, delimiter = ',')
    np.savetxt(NameFileD, DSet, delimiter = ',')