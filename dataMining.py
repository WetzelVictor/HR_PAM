#-*-encoding:UTF8-*-
from dataset_builder import utils
from dataset_builder import audio_engine
import numpy as np
import matplotlib.pyplot as plt

#%% Computing dataStructure
datasetFolder = 'dataset'
dataStruct = utils.build_data_structure(datasetFolder)
dataStruct = audio_engine.compute_features(dataStruct)
Nx = len(dataStruct["filepath"])

#%% Plots
Labels = ['P1','P2','P3']
#Labels = ['nylon','acier']
plotX = 'SpectralRolloffMean'
plotY = 'SpectralCentroidMean'

# A FAIRE
classesIncludeOnly = ['noise','fil']
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
    
    if not "noise" in curLabel and not "accel" in curLabel:
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
