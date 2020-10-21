# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 09:57:28 2015

@author: pichugin
"""
import sys
import numpy as np
#import pandas as pd
import LifeCycleSupplementary_v1_3_3 as LC

"the main code"

##Parse parameter
#if len(sys.argv) < 2:
#    print("No parameter")
#    sys.exit()

job_id = sys.argv[1]
#job_id='0'

N_MAX=19


LossDelta=1
LossList=np.arange(N_MAX)


NameFileWrite = 'Results/DetrimentalRandom_loss_costs_'+job_id+'.txt'
File2Write = open(NameFileWrite, 'w')

""" write the first line with risk values """
LossLine = 'ID, '
for i in LossList:
    LossLine = LossLine+str(i)+', '
LossLine=LossLine[0:-2]+'\n'
File2Write.write(LossLine)

""" files with fitness landscapes """
NameFileB = '../Input data/B_profiles_Detrimental_'+job_id+'.txt'
NameFileD = '../Input data/D_profiles_Detrimental_'+job_id+'.txt'

FileB=open(NameFileB,'r')
FileD=open(NameFileD,'r')

Bline = FileB.readline()
Dline = FileD.readline()

while Bline!='':
    """ convert lines in files into arrays of floats """
    Bsplit=Bline.split(',')
    Dsplit=Dline.split(',')
    
    ID_str=Bsplit[0]
    
    Bsplit = Bsplit[1:len(Bsplit)]
    Dsplit = Dsplit[1:len(Dsplit)]    
    
    B=[]
    for Bentry in Bsplit:
        B.append(float(Bentry))
    
    D=[]
    for Dentry in Dsplit:
        D.append(float(Dentry))
    
    LC_string=ID_str+', '    
    for Loss in LossList:
        LC_ID=LC.PartList_Loss(N_MAX, Loss)
        """ Find the optimal life cycle """
        CurMaximum = -100
        for Reaction in range(len(LC_ID)): 
            LifeCycle=LC_ID[Reaction][1]
            ProliferationRate=LC.Eig(LifeCycle,B,D)
            if ProliferationRate>CurMaximum:
                CurMaximum=ProliferationRate
                BestLC=LC_ID[Reaction][0]
                
        LC_string = LC_string + str(BestLC) +', '
    
    LC_string = LC_string[0:-2]+'\n'
    File2Write.write(LC_string)

    Bline = FileB.readline()
    Dline = FileD.readline()


FileB.close()
FileD.close()
File2Write.close()