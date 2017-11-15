# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 12:27:54 2017

@author: maryluz.mouronte
"""

"""
Editor de Spyder

Este es un archivo temporal
"""

import pandas as pd
import numpy as np
import os
import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
import operator
import random
from random import randint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from os import listdir
import shutil
import time
from collections import OrderedDict
import copy
import glob



def ConstruirGrafo(fileName):

 
    print("Building graph")
    
    pathFile = os.getcwd() + "\\" + "..\\datos\\" + fileName
                    
    with open(pathFile) as f:
        lines = f.readlines()
        
    linesN=[]
    for line in lines:
       line= line.replace('\x00','')
       if (len(line) != 0):
           linesN.append(line)
    lines = linesN

    myList = [line.strip().split() for line in lines]
    # [['a', 'b'], ['a', 'c'], ['b', 'd'], ['c', 'e']]
    
    
    g = nx.Graph()
    #g.add_edges_from(myList)
    g.add_weighted_edges_from(myList)
    
    return g

def CalcularAdyacencia (Grafo, fileName):

    pathFile = os.getcwd() + "\\" + fileName
    if path.exists(pathFile):
        os.remove(pathFile)

    df = nx.to_pandas_dataframe(Grafo)
    df.to_csv(pathFile, header=None, index=None, mode='w', sep='\t')

def LeerFicheros():
 
    
    Header = ['year','origin','destiny','imp','exp']
    fileName = '..\datos\Red_V.csv'
    pathFile = os.getcwd() + "\\" + fileName
    print('Leyendo fichero')
    
    
    datosRed = pd.read_table(pathFile, 'engine=python', delimiter=';', header=0, encoding = "ISO-8859-1", names=Header)
    print(datosRed.head(5))
    anios = datosRed["year"]
    anios = anios.drop_duplicates()
    for i in anios:
        #Calculamos la red de negocios
        datosRedCom = datosRed[datosRed.year == i]
        maxCom = datosRedCom["imp"].max()
        print("Maximo")
        print(maxCom)
        print("Antes")
        print(datosRedCom.head(5))
        datosRedCom.loc[datosRedCom['imp']>0, 'imp']= round((datosRedCom["imp"]/maxCom)*100,2)
        print("Despues")
        print(datosRedCom.head(5))
        datosRedCom = datosRedCom.drop(labels="year", axis=1)
        datosRedCom = datosRedCom.drop(labels="exp", axis=1)
        datosRedCom = datosRedCom.drop_duplicates()
        

        fileName = '..\datos\RedCom'+str(i)+'.txt'
        pathFile = os.getcwd() + "\\" + fileName
        if path.exists(pathFile):
           os.remove(pathFile)  
        datosRedCom.to_csv(pathFile, header=None, index=None, mode='w', sep=' ')
        g=ConstruirGrafo('..\datos\RedCom'+str(i)+'.txt')
        fileName = '..\datos\RedAdyCom'+str(i)+'.txt'
        CalcularAdyacencia (g, fileName)
    
        
    
    #for pathFile in glob.glob(os.getcwd() + "\\"+"..\datos\RedCom*.txt", recursive=True): 
    #    os.remove(pathFile)  
    

        
print("Calculating network")
LeerFicheros()  








