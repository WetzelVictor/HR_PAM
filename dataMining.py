#-*- coding: utf-8 -*-
from dataset_builder import utils
from dataset_builder import audio_engine
import numpy as np
import matplotlib.pyplot as plt

#%% Computing dataStructure
datasetFolder = 'output/extracted/'
dataStruct = utils.build_data_structure(datasetFolder)
dataStruct = audio_engine.compute_features(dataStruct)
Nx = len(dataStruct["filepath"])

#%% Plots
#Labels = ['plectre','index','pouce']
#Labels = ['nylon','acier']
#Labels = ['0deg','45deg','90deg']
#Labels = ['plexiglas','epicea']
Labels = ['plectre','fil']

plotX = 'PerceptualSharpnessMean'
plotY = 'PerceptualSpreadMean'

# A FAIRE
classExclude = ['nylon']
# A FAIRE

# log ou None
plotType = None

# INIT
XX = {}
YY = {}

# LOOP
for tag in Labels:
    XX[tag] = []
    YY[tag] = []

for i in xrange(Nx):
    curLabel = dataStruct['assigned_label'][i]
    
    if True in np.isin(curLabel,classExclude):
        continue
    
    for tag in Labels:
        if tag in curLabel:
            XX[tag].append(dataStruct[plotX][i])
            YY[tag].append(dataStruct[plotY][i])

# PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLO
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(12,7))  
for tag in Labels:
    if plotType == 'log':
        plt.scatter(10*np.log(XX[tag]), 10*np.log(YY[tag]))
    elif plotType is None:
        plt.scatter(XX[tag],YY[tag])

               
plt.xlabel(plotX,fontsize=20)
plt.ylabel(plotY,fontsize=20)
plt.legend(Labels,fontsize=18)

#%%
plotX = 'SpectralSlopeMean'
plotY = 'SpectralCentroidMean'

# A FAIRE
classExclude = ['noise','fil']
# A FAIRE

# log ou None
plotType = None

# INIT
XX = {}
YY = {}

# LOOP
for tag in Labels:
    XX[tag] = []
    YY[tag] = []

for i in xrange(Nx):
    curLabel = dataStruct['assigned_label'][i]
    
    if not "sinus" in curLabel or "accel" in curLabel or not "repetabilite" in curLabel:
        continue
    
    for tag in Labels:
        if tag in curLabel:
            XX[tag].append(dataStruct[plotX][i])
            YY[tag].append(dataStruct[plotY][i])

# PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLO
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(12,7))  

plt.scatter(XX[tag],YY[tag])
plt.errorbar(np.mean(XX[tag]),np.mean(YY[tag]),xerr=np.std(XX[tag]), yerr=np.std(YY[tag]),fmt='o',ecolor='g',markersize=20,capsize=20,mec='b',mfc='r',capthick=2, dash_capstyle='butt')

plt.title('Moyenne et ecart-type des descripteurs %s et %s'%(plotX,plotY),fontsize=20)
plt.xlabel(plotX,fontsize=20)
plt.ylabel(plotY,fontsize=20)
plt.xlim([-0.000001,-0.000001])
plt.ylim([300,500])