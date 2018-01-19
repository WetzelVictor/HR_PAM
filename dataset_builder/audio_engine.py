#-*-encoding:UTF8-*-
import librosa as lib
import yaafelib as yaf
import numpy as np

def compute_features(dataStruct):
    """ This function takes a data structure dictionnaire, and renders several
    audio features as spectral rolloff, spectral slope etc... and store the
    data into the datastructure.

    Args:
        - dataStruct: dictionnaire containing filepath, labels, and list of
          classes

    Returns:
        - dataSet: same as dataStruct, with the given spectral features
    """
    
    ### --- INIT --- ###
    # DSP settings
    Nwin_bin = 1024
    Hop_bin = round(Nwin_bin)
    
    # Const
    Nex = len(dataStruct["filepath"]) # Number of files

    # Listing audio features
    features_yaafe = ['SpectralFlatness',
                      'SpectralRolloff',
                      'SpectralSlope',
                      'SpectralDecrease',
                      'SpectralVariation',
                      'SpectralFlux']
    features_libro = ['Loudness', 
                      'SpectralCentroid', 
                      'SpectralContrast', 
                      'SpectralRolloff',
                      'SpectralBandwidth'] 
    
    dataStruct["SpectralFeatures"] = features_yaafe + features_libro
    
    # New fields
    dataStruct["signal"] = []
    dataStruct["sRate"] = []
    
    # Creating three fields per descriptor: full temporal vector, mean, and
    # standard deviation
    for f in dataStruct["SpectralFeatures"]:
        dataStruct[f] = []
        dataStruct[f + 'Mean'] = []
        dataStruct[f + 'Std'] = []
        dataStruct[f + 'Max'] = []


    ### --- Compute Feature --- ###
    print('\t \t \t Feature Extraction')
    # Computing the set of features
    for curFile in range(Nex):         
        print('%s' % dataStruct["filepath"][curFile])

        # Loading signal
        curSignal, curSRate = lib.load(dataStruct["filepath"][curFile], 
                                       mono=True, 
                                       offset=0)
            
        # Storing signal data
        dataStruct["signal"].append(curSignal)
        dataStruct["sRate"].append(curSRate)
        
        """ YAAFE Extraction """
        # Create YAAFE extraction engine
        fp = yaf.FeaturePlan(sample_rate=curSRate)
        
        # Formatting string for DSP
        for f in features_yaafe:
            fp.addFeature(f+': '+f+' blockSize='+str(Nwin_bin)+\
                          ' stepSize='+str(Hop_bin))
        
        engine = yaf.Engine()
        engine.load(fp.getDataFlow())
        features = engine.processAudio(curSignal.astype('float64')\
                                       .reshape((1, curSignal.shape[0])))
        
        # Computing mean and std for each 
        for key, val in sorted(features.items()):
            dataStruct[key].append(val)
            dataStruct[key + 'Mean'].append(np.mean(val))
            dataStruct[key + 'Std'].append(np.std(val))
            dataStruct[key + 'Max'].append(np.max(val))
        
        """ Librosa extraction """
        # Add the specific features from Librosa
        dataStruct["Loudness"].append(lib.feature.rmse(curSignal))

        # Compute the spectral centroid. [y, sr, S, n_fft, ...]
        dataStruct["SpectralCentroid"].append(
                                    lib.feature.spectral_centroid(curSignal))
        
        # Compute spectral contrast [R16] , sr, S, n_fft, ...])	
        dataStruct["SpectralContrast"].append(
                                    lib.feature.spectral_contrast(curSignal))
        
        # Compute roll-off frequency
        dataStruct["SpectralRolloff"].append(
                                    lib.feature.spectral_rolloff(curSignal))
        
        # Compute Bandwidth
        dataStruct["SpectralBandwidth"].append(
                                    lib.feature.spectral_bandwidth(curSignal))

        # Computing mean and std for each 
        for f in features_libro:
            val = dataStruct[f][-1]
            dataStruct[f + 'Mean'].append(np.mean(val))
            dataStruct[f + 'Std'].append(np.std(val))
            dataStruct[f + 'Max'].append(np.max(val))
 
    ### --- Formatting --- ###
    return dataStruct
