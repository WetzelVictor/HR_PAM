#-*-encoding:UTF-8-*-

# IMPORT STATEMENTS
import os
import matlab.engine

#%% Booting matlab's engine
eng = matlab.engine.start_matlab()

# current working directory
cwd = os.getcwd() + '/'

# You can specify your own path here
inputpath = cwd + 'data/repetabilite/'
outputpath = cwd + 'output/repetabilite/'

#%% Goes through the whole database
for dirpath, dirnames, filenames in os.walk(inputpath):
    for filename in filenames:
        if filename.endswith('.wav'):
            # Creates filepath for input file, and output file
            temp_input  = dirpath + '/' + filename
            temp_output = temp_input.replace(inputpath, outputpath) 
            dir_out = dirpath.replace(inputpath, outputpath) + '/'

            fig_name = temp_output.replace('.wav','.fig')
            
            # Applies matlab function to extract signal and noise
            eng.space_separation(temp_input, 
                                 dir_out, 
                                 filename,
                                 fig_name,
                                 nargout=0) 

            print('File ' + temp_output + ' processed')
            
            ### END OF THE LOOP ###


print('done')

