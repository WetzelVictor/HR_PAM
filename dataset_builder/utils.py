#-*-encoding:UTF8-*-

""" IMPORT STATEMENTS """
# lib
import sys
import os

def build_data_structure(database_folder):
    """ builds a data structure that contains: a list of labels based on the
    names of folder and the filename; the assignement of the differents file to
    a label; the list of audio file's path.

    Args:
        - database_folder: folder that contain the first architechture

    Returns:
        - dataStruct [dict]["filenames","assigned_labels","labels"]:
    assigned_labels gives the index in "labels" for every label

    """
    # --- INIT ---
    class_names = []
    file_paths = []

    # Build raw list of classes
    for dirpath, dirnames, filenames in os.walk(database_folder):
        for filename in filenames:
            temp_str = dirpath + '/' + filename
            if not '%' in temp_str and filename.endswith('.wav'):
                file_paths.append(temp_str)
                # Raw label list from filepath
                class_raw = temp_str.split('/')
                for classe in class_raw:
                    temp = classe.split('_')
                    for name in temp:
                        if name not in class_names:
                            class_names.append(name)

    # Remove unnecessary labels
    try:                            
        class_names.remove(database_folder)
    except ValueError:
        print 'Couldn\'t remove root file name'
        
    class_names.remove('.wav')
    to_remove = []
    Nclass = len(class_names)

    # list undesirable names
    for i in xrange(Nclass):
        if class_names[i].endswith('.wav'):
            to_remove.append(class_names[i])
    # removes them
    for name in to_remove:
        class_names.remove(name)


    ''' -=- Building data structure -=- '''
    Nex = len(file_paths)
    assigned_class = [[] for i in xrange(Nex)]

    for i,filepath in enumerate(file_paths):
        for j,label in enumerate(class_names):
            if label in filepath:
                assigned_class[i].append(j)
    
    ''' -=- formatting -=- '''
    dataStruct = {}
    dataStruct["filepath"] = file_paths
    dataStruct["assigned_label_index"] = assigned_class
    dataStruct["labels"] = class_names
    dataStruct["assigned_label"] = []
    
    # Makes the list of every label (in str) for each file
    for i in xrange(Nex):
        N_assLab = len(dataStruct['assigned_label_index'][i])
        curLabel = [dataStruct['labels'][dataStruct['assigned_label_index'][i][j]] for j in xrange(N_assLab)]
        dataStruct["assigned_label"].append(curLabel)
        

    return dataStruct
