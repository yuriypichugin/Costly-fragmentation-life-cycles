# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 09:57:28 2015

@author: Yuriy Pichugin

The collection of functions used to calculate the fitness landscape, generate life cycle structure and calculate the proliferation rate of the life cycle under the given fitness landscape
"""


"""
The life cycle data structure.

The life cycle is encoded as a list of pathways.
Each pathway is represented by a tuple: (Fragmented_Size, [Fragments_Sizes], Fragmentation_Rate, Trigger_Flag)

(int) Fragmented_Size - the size of the fragmented group
(int list) Fragments_Sizes - the list of size of fragments
(float) Fragmentation_Rate - the rate of fragmentation (for free pathways), or the probability of fragmentation (for triggered pathways)
(int) Trigger_Flag - the flag indicating the motive of the fragmentation

TriggerFlag==0 triggered pathways
TriggerFlag==1 free pathways (can be used by ProjectionMatrix(LifeCycle, MMSizes, b, d))

The sum of Fragments_Sizes should not exceed Fragmented_Size.
For a fragmentation without loss, this sum is equal to Fragmented_Size.

The sum of Fragmentation_Rate values for all triggered pathways with the same Fragmented_Size should not exceed one.
The sum of Fragmentation_Rate values for all triggered pathways with the maximal Fragmented_Size should be equal to one.
"""

import numpy as np
from IntegerPartitions import revlex_partitions


"""
PART 1: calculation of the proliferation rate as the largest eigenvalue of the projection matrix
"""

def MinMaxSizes(LifeCycle):
    """
    calculates the minimal and maximal sizes of groups 
    that appears in the life cycle.
    This is later used for the projection matrix construction
       
    Arguments:
    (list of tuples) LifeCycle - life cycle data, see description of the format in the separate comment
       
    Return:
    The list of two integers: the first one is the minimal size of groups, the second one is the maximal size of groups in the course of a given life cycle
    """
    MaxSize=0
    for pathway in range(len(LifeCycle)):
        Fragmenting_Group_Size = LifeCycle[pathway][0]
        TriggerFlag = LifeCycle[pathway][3]
        if MaxSize < (Fragmenting_Group_Size - 1 + TriggerFlag):
            MaxSize = (Fragmenting_Group_Size - 1 + TriggerFlag)
    
    MinSize=MaxSize
    for pathway in range(len(LifeCycle)):
        Offsprings=LifeCycle[pathway][1]
        MinOffspring=min(Offsprings)
        if MinOffspring<MinSize:
            MinSize=MinOffspring
        
    return [MinSize, MaxSize]


def ProjectionMatrix(LifeCycle, MMSizes, b, d):
    """
    calculates the projection matrix of a given life cycle
    
    Arguments:
    
    (list of tuples) LifeCycle - life cycle data, see description of the format in the separate comment
    (list of ints) MMSizes - the list of the minimal and the maximl sizes of groups in a given life cycle
    (list of floats) b - the list of the growth rates
    (list of floats) d - the list of the death rates
    
    Return:
    (matrix) A - projection matrix corresponding to a given life cycle in a given fitness landscape
    """
    MatrixSize=MMSizes[1]-MMSizes[0]+1
    A=np.zeros((MatrixSize, MatrixSize))
    A=np.asmatrix(A)
    
    #set up the fitness landscape
    for i in range(MatrixSize):
        A[i,i]=A[i,i]-b[i+MMSizes[0]-1]
        A[i,i]=A[i,i]-d[i+MMSizes[0]-1]
        if i < (MatrixSize-1):
            A[i+1,i]=A[i+1,i]+b[i+MMSizes[0]-1]
    
    #include the impact of fragmentation pathways
    for pathway in range(len(LifeCycle)):
    
        Pathway=LifeCycle[pathway]
    
        if Pathway[3]==0:        #triggered pathways
            Rate=Pathway[2]*b[Pathway[0] - 2]
            Column=Pathway[0] - 2 - (MMSizes[0]-1)
        elif Pathway[3]==1:      #free pathways
            Rate=Pathway[2]
            Column=Pathway[0] - 1 - (MMSizes[0]-1)
        
        for i in range(len(Pathway[1])):
            Row=Pathway[1][i] - 1 - (MMSizes[0]-1)
            A[Row, Column] = A[Row, Column] + Rate
        
        if Pathway[3]==0:
            if Column < (MatrixSize-1):
                A[Column+1,Column]=A[Column+1,Column] - Rate
        elif Pathway[3]==1:
            A[Column,Column]=A[Column,Column] - Rate
        
    return A


def Eig(LifeCycle, b, d):
    """
    calculates the proliferation rate provided by a life cycle under given fitness landscape
    
    Arguments:
    
    (list of tuples) LifeCycle - life cycle data, see description of the format elsewhere
    (list of floats) b - the list of the growth rates
    (list of floats) d - the list of the death rates
    
    Return:
    (float) ProliferationRate - proliferation rate provided by the life cycle in a given fitness landscape
    """    
       
    #Find the minimal and the maximal size of groups emerging in the given life cycle
    MMSizes=MinMaxSizes(LifeCycle)
    
    #check 1: the length of b and d vectors should be at least MaxSize 
    if (len(b)<MMSizes[1]):
        print ("Error: vector b is too short for the life cycle")
        return float('NaN')
    if(len(d)<MMSizes[1]):
        print("Error: vector d is too short for the life cycle")
        return float('NaN')
        
        
    
    #Initialize the projection matrix
    A=ProjectionMatrix(LifeCycle, MMSizes, b, d)
    
    #Calculate the proliferation rate
    Eigenvalues=np.linalg.eig(A)[0]
    ProliferationRate=max(Eigenvalues.real)

    return ProliferationRate
    
    
"""
PART 2: Initialization of the fitness landscape in a form of birth and death rate vectors
"""
    
def BirthInit(model, NMax, Params):
    """
    Initialization of the birth rates vector
    
    Arguments:
    (int) model - the number of model used
    (int) NMax  - the size of vector (the maximal size of group existing in a population)
    (list of floats) Params - parameters of the fitness landscape
    
    Return:
    (list of floats) b - the vector of birth rates (per group size!)
    """
    
    b=[0]*NMax
    
    if model==1:
        """Small sizes scans model. Rates are equal to Params"""
        b[0]=1
        b[1]=Params[0]
        b[2]=Params[1]
    elif model==2:
        """ Power model. """
        M=Params[0]
        S=Params[1]
        for i in range(NMax):
            i_eff=i/(NMax-1.0)
            b[i]=(i+1)*(1+M*pow(i_eff,S))
    elif model==3:
        """ Random uniform model. """
        for i in range(NMax):
            b[i]=(i+1)*np.random.uniform(0.0, 1.0)
    elif model==4:
        """ Random exponential model """
        Lambda = Params[0]
        for i in range(NMax):
            b[i]=(i+1)*np.random.exponential(Lambda)
    elif model==5:
        """ Random gauss model """
        Mean = Params[0]
        Std = Params[1]
        for i in range(NMax):
            b[i]=(i+1)*np.random.normal(Mean, Std)
    elif model == 6:
        """ monotonically increasing random sequence """
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)
        for i in range(NMax):
            b[i]=(i+1)*true_bi_sorted[i]
    elif model == 7:
        """ monotonically decreasing random sequence """
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)[::-1]
        for i in range(NMax):
            b[i]=(i+1)*true_bi_sorted[i]
    elif model == 8:
        """ unimodal random sequence """
        LocMax = Params[0]
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)[::-1]
        LeftSlot = LocMax - 1 
        RightSlot = LocMax + 1 
        unimodal_true_bi = np.zeros(NMax)
        unimodal_true_bi[LocMax] = true_bi_sorted[0]
        b[LocMax] = (LocMax+1)*unimodal_true_bi[LocMax]
        for i in range(1, NMax):
#            print([LeftSlot, RightSlot])
            if (LeftSlot >= 0) and (RightSlot < NMax):
                coin = np.random.randint(2)
                if coin == 0:
                    position = LeftSlot
                    LeftSlot -= 1
                else:
                    position = RightSlot
                    RightSlot += 1
            elif (LeftSlot < 0) and (RightSlot < NMax):
                position = RightSlot
                RightSlot += 1
            elif (LeftSlot >= 0) and (RightSlot >= NMax):
                position = LeftSlot
                LeftSlot -= 1
            else:
                position = 0
                print(" profile generation went wrong ! ")
#            print(position)
            unimodal_true_bi[position] = true_bi_sorted[i]
            b[position] = (position+1)*unimodal_true_bi[position]
    else:
        """ default model, no effect of the group size on the growth of cells """
        for i in range(NMax):
            b[i]=i+1
    return b 
    
def DeathInit(model, NMax, Params):
    """
    Initialization of the death rates vector
    
    Arguments:
    
    (int) model - the number of model used
    (int) NMax  - the size of vector (the maximal size of group existing in a population)
    (list of floats) Params - parameters of the fitness landscape
    
    Return:
    (list of floats) d - the vector of death rates
    """
    d=[0]*NMax
    
    if model==1:
        """Small sizes scans model. Death is determined by params """
        d[0]=0
        d[1]=Params[0]
        d[2]=Params[1]
    elif model==2:
        """ Power model. """
        M=Params[0]
        S=Params[1]
        for i in range(NMax):
            i_eff=i/(NMax-1.0)
            d[i]=-M*pow(i_eff,S)
    elif model==3:
        """ Random uniform model. """
        for i in range(NMax):
            d[i]=np.random.uniform(0.0, 1.0)
    elif model==4:
        """ Random exponential model """
        Lambda = Params[0]
        for i in range(NMax):
            d[i]=np.random.exponential(Lambda)
    elif model==5:
        """ Random gauss model """
        Mean = Params[0]
        Std = Params[1]
        for i in range(NMax):
            d[i]=np.random.normal(Mean, Std)
    elif model == 6:
        """ monotonically increasing random sequence """
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)
        for i in range(NMax):
            d[i]=di_sorted[i]
    elif model == 7:
        """ monotonically decreasing random sequence """
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)[::-1]
        for i in range(NMax):
            d[i]=di_sorted[i]
    elif model == 8:
        """ unimodal random sequence """
        LocMax = Params[0]
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)
        LeftSlot = LocMax - 1 
        RightSlot = LocMax + 1 
        unimodal_di = np.zeros(NMax)
        unimodal_di[LocMax] = di_sorted[0]
        d[LocMax] = unimodal_di[LocMax]
        for i in range(1, NMax):
#            print([LeftSlot, RightSlot])
            if (LeftSlot >= 0) and (RightSlot < NMax):
                coin = np.random.randint(2)
                if coin == 0:
                    position = LeftSlot
                    LeftSlot -= 1
                else:
                    position = RightSlot
                    RightSlot += 1
            elif (LeftSlot < 0) and (RightSlot < NMax):
                position = RightSlot
                RightSlot += 1
            elif (LeftSlot >= 0) and (RightSlot >= NMax):
                position = LeftSlot
                LeftSlot -= 1
            else:
                position = 0
                print(" profile generation went wrong ! ")
#            print(position)
            unimodal_di[position] = di_sorted[i]
            d[position] = unimodal_di[position]
    else:
        """ default model, no effect of the group size on the death of cells """
        for i in range(NMax):
            d[i]=0
    return d


"""
PART 3: Construction of a list of all partitions up to certain size
"""

def PartList(N):
    """
    Calculation of the partition list in a fragmentation without losses.
    The sum of offspring sizes is equal to the size of the fragmented group
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for Size in range(2, N+2): 
        #create generator of partitions of Size (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(Size)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            LifeCycle.append([(Size, PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    LC_ID=[]
    for counter in range(len(LifeCycle)):
        LC_ID.append([counter,LifeCycle[counter]])

    return LC_ID
    
def PartList_Loss(N, Loss):
    """
    Calculation of the partition list in a fragmentation with given losses.
    This is intended to replace PartList_WL().
    The sum of offspring sizes is equal to the size of the fragmented group minus Loss
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    (int) Loss - the number of cells lost in a fragmentation
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for SizeOff in range(2, N+2-Loss): 
        #create generator of partitions of SizeOff (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(SizeOff)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            LifeCycle.append([(SizeOff+Loss, PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    LC_ID=[]
    for counter in range(len(LifeCycle)):
        LC_ID.append([counter,LifeCycle[counter]])

    return LC_ID

def PartList_PropLoss(N, Loss):
    """
    Calculation of the partition list in a fragmentation with given proportional losses.
    This case represent the fragmentation of linear filaments
    The sum of offspring sizes is equal to the size of the fragmented group minus Loss*(#parts-1)
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    (int) Loss - the number of cells lost at each break of the chain
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for SizeOff in range(2, N+2): 
        #create generator of partitions of SizeOff (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(SizeOff)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            Num_of_parts=len(PartList[Reaction])            
            LifeCycle.append([(SizeOff+Loss*(Num_of_parts - 1), PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    pre_LC_ID=[]
    for counter in range(len(LifeCycle)):
        pre_LC_ID.append([counter,LifeCycle[counter]])

    LC_ID=[]
    for ID in pre_LC_ID:
        if ID[1][0][0]<(N+2):
            LC_ID.append(ID)
    
    return LC_ID


def PartList_WL(N):
    """
    Calculation of the partition list in a fragmentation with losses
    The sum of offspring sizes is one less than the size of the fragmented group
        
    Arguments:
    (int) N - the maximal size of group presented in population.
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for Size in range(2, N+1): 
        #create generator of partitions of Size (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(Size)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            LifeCycle.append([(Size+1, PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    LC_ID=[]
    for counter in range(len(LifeCycle)):
        LC_ID.append([counter,LifeCycle[counter]])
    
    return LC_ID
    
    
"""
PART 4: Calculation of the proliferation rate by the pure strategy by simplified formula.
"""   
def F(GRate, birth, death, size, MinimalSize):
    """calculates the F terms for the grwoth rate computation
        F_i=(1+(lambda+d_i)/b_i)
        GRate [double] is the assigned population growth rate
        birth [double vector] is the vector of cell birth rates
        death [double vector] is the vector of group death rates
        size [int] is the index of F
        MinimalSize [int] is the minimal size of groups existing in a population
    """
    FV=1
    for i in range(MinimalSize-1, size-1):
        FV=FV*(1+(GRate+death[i])/birth[i])
    
    return FV
    
def CharExp(GRate, birth, death, k):
    """calculates the characteristic expression for a pure life cycle 
        GRate [double] is the assigned population growth rate
        birth [double vector] is the vector of cell birth rates
        death [double vector] is the vector of group death rates
        k     [int vector] is the fragmentation pattern vector
    """
    MinimalSize=len(k)
    for i in range(len(k)):
        if k[i]!=0 and (i+1)<MinimalSize:
            MinimalSize=i+1
        
#    MinimalSize=min(NProgenitor-NOffspring, NOffspring)    
    Result=0
    for i in range(len(k)):
        Result=Result+k[i]*F(GRate, birth, death, i+1, MinimalSize)

    Result = -Result #to be consistent with the paper
    
    return (Result)

def Bisection(birth, death, k):
    """calculates the growth rate for a life cycle using the pathway (NProgenitor) -> (NOffspring) + (NProgenitor - NOffspring)
        by solving CharExp(GRate)=0 with bisection
        
        birth [double vector] is the vector of cell birth rates
        death [double vector] is the vector of group death rates
        k     [int vector] is the fragmentation pattern vector
    """
    MinimalSize=len(k)
    for i in range(len(k)):
        if k[i]!=0 and (i+1)<MinimalSize:
            MinimalSize=i+1
    DA=np.array(death)
    BA=np.array(birth)
    LowGR=-min(BA[MinimalSize-1:(len(k)-1)]+DA[MinimalSize-1:(len(k)-1)])
#    LowGR=0.0
#    LowGR=-death[NProgenitor-1]
#    LowGR=-max(death)
#    LowGR=min([0,-min(death)])
    LowCE=CharExp(LowGR, birth, death, k)
    
    HighGR=1000000.0
    HighCE=CharExp(HighGR, birth, death, k)
        
    if LowCE*HighCE<=0:
        while HighGR-LowGR>0.0000001*(abs(HighGR)+abs(LowGR)):
            TestGR=(LowGR+HighGR)/2.0
            TestCE=CharExp(TestGR, birth, death, k)
            if TestCE*LowCE<=0:
                HighGR=TestGR
                HighCE=TestCE
            else:
                LowGR=TestGR
                LowCE=TestCE
        return (LowGR+HighGR)/2
    else:
        return 0  

"""
PART 5: Generation of random stochastic life cycle
"""   


def Generate_Random_Stochastic_LC(Nmax):
    """
    Generates life cycle as a list of tuples
    
    Arguments:
    Nmax - the maximal size of groups EXISTING in the population
    
    Return:
    LifeCycle - life cycle structure
    """
    TrigProbMagnitude=np.random.uniform(0.0, 1.0, Nmax+1)     #the sum of probabilities of all rates of triggered fragmentation occurring with groups of a each size
    
    LifeCycle=[]
    #set up triggered reactions
    for Size in range(2, Nmax+2): 
        #create generator of partitions of Size (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(Size)
        #convert generator into a list 
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        
        #compute fragmentation probabilities
        Probs=np.random.exponential(1, len(PartList))
        if Size < Nmax+1:
            Norm=TrigProbMagnitude[Size]/sum(Probs[1:])
        else:                                           #the combined fragmentation probability at the maximal size is equal to one
            Norm=1/sum(Probs[1:])
        Probs=Probs*Norm
        
        #Construct biological reactions data
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            RData=(Size, PartList[Reaction], Probs[Reaction], 0)        
            LifeCycle.append(RData)
    
    return LifeCycle