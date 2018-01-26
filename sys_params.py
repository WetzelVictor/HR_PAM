#-*-encoding:UTF-8-*-

# IMPORT STATEMENTS
import os
import matlab.engine

#%% Booting matlab's engine
eng = matlab.engine.start_matlab()

# current working directory
cwd = os.getcwd() + '/'

# You can specify your own path here
inputpath = cwd + 'output/extracted/mesures/plexiglas/acier'

#%% Goes through the whole database
for dirpath, dirnames, filenames in os.walk(inputpath):
    for filename in filenames:
        if filename.endswith('.mat'):
            # Creates filepath for input file, and output file
            temp_input  = dirpath + '/' + filename
            fig_name = temp_input.replace('.mat','_modes_struct_reg.fig')
            
            # Applies matlab function to extract signal and noise
            eng.plot_modes_struct(temp_input, 
                            filename, fig_name, nargout=0) 

            print('File ' + temp_input + ' processed')
            
            ### END OF THE LOOP ###


print('done')

