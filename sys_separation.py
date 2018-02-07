#-*-encoding:UTF-8-*-

# IMPORT STATEMENTS
import os
import matlab.engine

def replicate_datastructure(input_path, output_path):
    """ This function replicated the architechture contained in 'input_path' in
    a different folder located at 'output_path'
    """
    for dirpath, dirnames, filenames in os.walk(input_path):
        if '%' in dirpath:
            continue
        
        structure = os.path.join(output_path, dirpath[len(input_path):])
        
        if not os.path.isdir(structure):
            os.mkdir(structure)
        else:
            print("Folder does already exits!")

#%% Booting matlab's engine
eng = matlab.engine.start_matlab()

#%% current working directory
cwd = os.getcwd() + '/'

# You can specify your own path here
inputpath = cwd + 'data/'
outputpath = cwd + 'output/extracted/'

# This function is defined at the beginning of the script
# replicate_datastructure(inputpath, outputpath)

#%% Goes through the whole database
for dirpath, dirnames, filenames in os.walk(inputpath):
    for filename in filenames:
        if filename.endswith('.wav'):
            # Creates filepath for input file, and output file
            temp_input  = dirpath + '/' + filename
            temp_output = temp_input.replace(inputpath, outputpath) 
            dir_out = dirpath.replace(inputpath, outputpath) + '/'

            fig_name = temp_output.replace('.wav','.fig')
            
            # Filename to check if the file exists already
            temp_check = dirpath + '/' + 'noise_' + filename
            temp_check = temp_check.replace(inputpath, outputpath) 
            
            # Applies matlab function to extract signal and noise
            if os.path.isfile(temp_check):
                print('file {0} already exists'.format(temp_check))
            else:
                eng.space_separation(temp_input, 
                                 dir_out, 
                                 filename,
                                 fig_name,
                                 nargout=0) 
                print('File ' + temp_output + ' processed')
            ### END OF THE LOOP ###


print('done')

