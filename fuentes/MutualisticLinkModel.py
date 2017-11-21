# -*- coding: utf-8 -*-
"""
Mutualistic Link Model for weighted mutualistic networks

The script has been written in Python 2.7, it requires numpy, scipy, 
pylab and matplotlib modulles to be installed.

Please cite:
    
    doi:    . (XXXX)
    
------------
Description
------------
The code provided includes the latest version of the Mutualistic Link Model 
generator (SimNet), a function for reading network.txt files (ReadWeb), a 
function for computing strength/degree vectors from interaction matrices (SumAP),
a function for computing strength/degree distributions from interaction matrices
(Distribution), a function for computing cumulative distributions 
(CumulativeDistribution), a function for averaging distributions (AverageDistribution),
a function for running multiple trials of the model (NsimDist) and two output
funtions: one for visual comparison of an empirical and a simulated interaction
matrix (CompareMatrix) and other for comparison of empirical and average simulated
distributions (CompareDistributions).

Input and output from the functions is described on each function's documentation.

The text files with the empirical matrices and this file must be in the same 
directory. Matrices text files can easily be obtained from the datafiles 
available on the web using the excel exporting tool, columns have to be separated
with tabs (\t) and rows with line jumps (\n).

"""
__author__ = 'Manuel Jiménez-Martín (manuel.jimenez@bec.uned.es) and Mary Luz Mouronte (maryluz.mouronte@ufv.es)'
__version__ = '1.0'

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
import glob


import random
from matplotlib.colors import LogNorm
from matplotlib.pyplot import imshow
from scipy import stats


#mlml 04/10/2017

import os
import os.path as path

#_____________________________________________________________________________
# INPUT DATA 
def ReadWeb(web_file):
    '''
    Reads interaction matrix of an empirical network from 'web_file' file. 
    It must be a .txt archive. The format must be plain text, with columns
    separated by tabs '\t' and rows separated by line jumps '\n'.
    Returns interaction matrix R as an array.
    '''
    archive = open(str(web_file),'r')
    List = archive.read()
    
    List = List.split('\n')
    
    List = [line.split('\t') for line in List]   # change '\t' for ' ' or ',' 
    #print(List)
    
    ListInt= [[float(y) for y in x ] for x in List if len(x)!=1]

    
    R = array(ListInt)
    R=R[~np.all(R==0,axis=1),:] # deletes columns of zeros
    R=R[:,~np.all(R==0,axis=0)] # deletes rows of zeros
    return R
    
#_____________________________________________________________________________    
# MODEL NETWORK BUILDER


def SimNet2(Nl,w,Na,Np,fl):

    '''
    Simulates a mutualistic network with w total weight, Na importers, 
    Np exporters and fl percentage of forbidden links
    Returns the weighted interaction matrix
    '''
    edge_list=[(0,0),(0,1),(1,0)]
    edge_list_weight= [(0,0,0.01),(0,1,0.01),(1,0,0.01)]
    Amax = 4
    Pmax = 4
    # generate random forbidden links matrix a priori
    FL = np.zeros((Na,Np),dtype=bool)   
    while sum(FL)/float(Na*Np) < fl:      
        a = random.choice(range(Na))     
        p = random.choice(range(Np))      
        if (a,p) not in edge_list and not FL[a,p]:
            FL[a,p] = 1

    # link aggregation process generates a link list          
    
    W = np.zeros((Na,Np),dtype=float)
    while ((sum(elem[2] for elem in edge_list_weight) < w) or (np.count_nonzero(W)<Nl)):
        if random.random() < ( (Na**2-Na)/(2.*w) )/(1.*Amax) and Amax < Na:
            A = Amax
            Amax += 1
        else:
            A = random.choice(edge_list)[0]
        if random.random() < ( (Np**2-Np)/(2.*w) )/(1.*Pmax) and Pmax < Np:
            P = Pmax
            Pmax += 1
        else:
            P = random.choice(edge_list)[1]
        if not FL[A,P]:  
            if (sum(elem[2] for elem in edge_list_weight) < w):
            	w_link= random.uniform(0.1,100)#Peso aleatorio entre 0 y 100
            else:
                w_link= random.uniform(0.01,0.1)
           
            edge_list.append((A,P))
            edge_list_weight.append((A,P,w_link))
            W[A,P] = W[A,P]+w_link
            print("Links in simulated network:")
            print(np.count_nonzero(W))
            
   
    return W



#_____________________________________________________________________________    
# DATA MANIPULATION & ANALYSIS FUNCTIONS

def SumAP(W,s):
    '''
    Returns A and P strength/degree-vectors from interaction matrix W
    (s==1 for strength and s==0 for degree)
    '''
    if s: # if s==1
        A = W.copy()
    else:
        A = W>0
    (N,M) = shape(A)
    SA = [sum(A[i,:]) for i in range(N)]    #for i in range(N):    
    SP = [sum(A[:,i]) for i in range(M)]    #    SA[i]=sum(A[i,:])
    return SA, SP 

def Pack(W):
    '''
    Returns W with rows and columns ordered from generalist (upper-left corner)
    to specialist (bottom-right corner)
    '''
    SA, SP = SumAP(W,0)
    W1 = W[argsort(SA)[::-1],:]
    return W1[:,argsort(SP)[::-1]] 
    
def Distribution(W,s):
    '''
    Returns A and P strength/degree distribution from interaction matrix W
    (x=s, y=P(s) if s==1 
     x=k, y=P(k) if s==0)
    '''
    SA, SP = SumAP(W,s)
    SA.sort()
    SP.sort()
    while SA[0] == 0:
        SA.pop(0)
    while SP[0] == 0:
        SP.pop(0)
    xA = []
    yA = []
    xP = []
    yP = []
    while len(SA) > 0:
        n = SA.count(SA[-1])
        xA.append(SA[-1])
        yA.append(n)
        SA = SA[:-n]
    while len(SP) > 0:
        n = SP.count(SP[-1])
        xP.append(SP[-1])
        yP.append(n)
        SP = SP[:-n]
    NA = sum(yA)
    NP = sum(yP)
    yA = [y/(1.*NA) for y in yA]
    yP = [y/(1.*NP) for y in yP]
    return xA, yA, xP, yP 
    
def CumulativeDistribution(W,s):
    '''
    Returns A and P strength/degree cumulative distributions from int. matrix W
    (x=s, y=P(s) if s==1 
     x=k, y=P(k) if s==0)
    '''
    xA, yA, xP, yP = Distribution(W,s)
    yA = [sum(yA[:i+1]) for i in range(len(xA))]
    yP = [sum(yP[:i+1]) for i in range(len(xP))]
    return xA, yA, xP, yP 

def AverageDistribution(X,Y,n):
    '''
    X,Y are x,y vectors of n distributions.
    For every value x in X, the function averages the Y values corresponding
    to every ocurrence of x  in X. Returns the average distribution meanX,
    meanY and the standard error .
    '''
    meanX = []
    meanY = []
    EY = []    
    # sort the distribution values in ascending order of X
    X,Y=(list(t) for t in zip(*sorted(zip(X,Y)))) 
    while len(X)>0:
        val = X[0]
        num = X.count(val)
        meanX.append(val)
        aux = list(zeros(n-num))
        meanY.append(mean(Y[0:num]+aux))
        EY.append(std(Y[0:num]+aux))
        X=X[num:]
        Y=Y[num:]
    return array(meanX[::-1]), array(meanY[::-1]), array(EY[::-1])  

#_____________________________________________________________________________    
# SIMULATION FUNCTIONS      

def NsimDist(Nl, w,Na,Np,fl,num_trials,s):
    '''
    Generates num_trials interaction matrices with parameters w, Na, Np and fl
    Returns average strength/degree distributions with standard errors and 
    average number of links (s==1 for strength and s==0 for degree).
    '''
    XA = []
    YA = []
    XP = []
    YP = []
    meanE = 0.
    for i in range(num_trials):
        W = SimNet2(Nl, w,Na,Np,fl)
        xA, yA, xP, yP = Distribution(W,s)
        XA += xA
        YA += yA
        XP += xP
        YP += yP
        meanE += sum(W>0)
    # average and acumulate distributions and statistical errors for both classes
    meanXA, meanYA, eYA = AverageDistribution(XA,YA,num_trials)
    meanYA = [sum(meanYA[:i+1]) for i in range(len(meanYA))]    
    eYA = [sqrt(sum(eYA[:i+1]**2)) for i in range(len(meanYA))]
    
    meanXP, meanYP, eYP = AverageDistribution(XP,YP,num_trials)
    meanYP = [sum(meanYP[:i+1]) for i in range(len(meanYP))]
    eYP = [sqrt(sum(eYP[:i+1]**2)) for i in range(len(meanYP))]
    
    meanE = meanE/(1.*num_trials)
    return array(meanXA[::-1]), array(meanYA[::-1]), array(eYA[::-1]), \
           array(meanXP[::-1]), array(meanYP[::-1]), array(eYP[::-1]),  meanE

#_____________________________________________________________________________    
# OUTPUT FUNCTIONS

def CompareMatrix(web_name,fl):
    '''
    Loads an empirical network from file 'web_name.txt' and runs the model
    with the same parameters w, Na, Np with a fl percentage of forbidden
    links. Represents and returns both matrices.
    '''
    R = ReadWeb(web_name)   
    (Na,Np) = shape(R)
    w = sum(R)
    Nl=np.count_nonzero(R) 
    print("Enlaces R")
    print(Nl)
    W = SimNet2(Nl, w,Na,Np,fl)
    print("Enlaces W")
    print(np.count_nonzero(W))
    R = array(Pack(R))
    W = array(Pack(W))
    fig, ax = subplots(1,2,sharex=False,sharey=False)
    m1=ax[0].imshow(R,interpolation='none',vmin=1,vmax=max(R.max(),W.max()),norm=LogNorm(),cmap='jet')
    m2=ax[1].imshow(W,interpolation='none',vmin=1,vmax=max(R.max(),W.max()),norm=LogNorm(),cmap='jet')
    cax = fig.add_axes([0.935,0.3,0.02,0.6]) 
    colorbar(m2,cax=cax,norm=LogNorm(),cmap='jet')    
    m1.axes.get_xaxis().set_ticks([])
    m2.axes.get_xaxis().set_ticks([])
    m1.axes.get_yaxis().set_ticks([])
    m2.axes.get_yaxis().set_ticks([])
    m1.axes.set_ylabel('I-importer')
    m1.axes.set_xlabel('E-exporter')
    m1.axes.set_title('Empirical int. matrix')
    m2.axes.set_ylabel('I-importer')
    m2.axes.set_xlabel('E-exporter')
    m2.axes.set_title('Simulated int. matrix')
    

    fileFig = web_name.replace(".txt","_"+str(fl)+"FigCMm1.png")
    myFig = m1.get_figure()  
    pathFile = fileFig         # Get the figure
    myFig.savefig(pathFile)  # Save to file
    
    fileFig = web_name.replace(".txt","_"+str(fl)+"FigCMm2.png")
    myFig = m2.get_figure()  
    pathFile = fileFig         # Get the figure
    myFig.savefig(pathFile)  # Save to file
    
    return R, W

def CompareDistributions(web_name,num_trials,fl):
    '''
    Loads an empirical network from file 'web_name.txt' and performs num_trials
    simulations of the model with same w, Na, Np as the empirical network and
    fl percentage of forbidden links. Represents the empirical and the average 
    simulated strength and degree distributions and number of links.
    '''
    R = ReadWeb(web_name)   
    (Na,Np) = shape(R)
    w = sum(R)
    Nl=np.count_nonzero(R) 
    W = SimNet2(Nl,w,Na,Np,fl)
    R = array(Pack(R))
    W = array(Pack(W))
    fig, ax = subplots(1,2,sharex=False,sharey=False)
    # strength distributions    
    sub = ax[1]
    sa, psa, sp, psp = CumulativeDistribution(R,1)
    xsa, ysa, esa, xsp, ysp, esp, me = NsimDist(Nl, w,Na,Np,fl,num_trials,1)    
    sub.loglog(sa,psa,'ro',sp,psp,'gs')
    sub.loglog(xsa,ysa,'r',xsp,ysp,'g')
    esa[esa>=ysa] = ysa[esa>=ysa]*.999 # prevents errorbars <= 0
    esp[esp>=ysp] = ysp[esp>=ysp]*.999 
    sub.fill_between(xsa,ysa+esa,ysa-esa,alpha=0.3,facecolor='gray',linewidth=0)
    sub.fill_between(xsp,ysp+esp,ysp-esp,alpha=0.3,facecolor='gray',linewidth=0)
    sub.set_ylim(0.8*min(min(psa),min(psp)),1.5)
    sub.set_xlim(0.9,1.5*max(max(sa),max(sp)))
    sub.set_xlabel('s')
    sub.set_ylabel('CDF(s)')
    # degree distributions
    sub = ax[0]
    ka, pka, kp, pkp = CumulativeDistribution(R,0)
    xka, yka, eka, xkp, ykp, ekp, me = NsimDist(Nl, w,Na,Np,fl,num_trials,0)    
    sub.loglog(ka,pka,'r*',kp,pkp,'gx') 
    sub.loglog(xka,yka,'r*',xkp,ykp,'gx') 
    eka[eka>=yka] = yka[eka>=yka]*.999 # prevents errorbars <= 0
    ekp[ekp>=ykp] = ykp[ekp>=ykp]*.999 
    sub.fill_between(xka,yka+eka,yka-eka,alpha=0.3,facecolor='gray',linewidth=0)
    sub.fill_between(xkp,ykp+ekp,ykp-ekp,alpha=0.3,facecolor='gray',linewidth=0)
    sub.set_ylim(0.8*min(min(pka),min(pkp)),1.5)
    sub.set_xlim(0.9,1.5*max(max(ka),max(kp)))
    sub.set_xlabel('k')
    sub.set_ylabel('CDF(k)')
    # empirical and average simulated number of links
    t = '$E_{emp}= '+str(sum(R>0))+'$\n'+'$E_{sim}= '+str(round(me))+'$'     
    sub.text(1.04,min(min(pka),min(pkp)),t)
    fig.tight_layout()
    

    fileFig = web_name.replace(".txt","_"+str(fl)+"FigCD.png")
    myFig = sub.get_figure()  
    pathFile = fileFig         # Get the figure
    myFig.savefig(pathFile)  # Save to file
    return 
    
#_____________________________________________________________________________    

def ExecuteExperiment(i):
#Executes experiment i

    global nameFile

    print("Path:"+nameFile)
    ReadWeb(nameFile)
    CompareDistributions(nameFile,30,i*0)
    CompareMatrix(nameFile,i*0)

for nameFile in glob.glob(os.getcwd() + "\\"+"..\datos\RedAdy*.txt", recursive=False): 
    print("nameFile")
    print(nameFile)
    [ExecuteExperiment(k) for k in range(1, 2)]