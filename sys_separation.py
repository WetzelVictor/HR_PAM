import os
import matlab.engine

eng = matlab.engine.start_matlab()

cwd = os.getcwd() + '/'
inputpath = cwd + 'data/'
outputpath = cwd + 'output/extracted/'

# This piece of code lists the filepath of .wav files in an architechture
for dirpath, dirnames, filenames in os.walk(inputpath):
    for filename in filenames:
        if filename.endswith('.wav'):
            # Creates filepath for input file, and output file
            temp_input  = dirpath + '/' + filename
            temp_output = temp_input.replace(inputpath, outputpath) 
            dir_out = dirpath.replace(inputpath, outputpath) + '/'
            
            # Applies matlab function to extract signal and noise
            eng.space_separation(temp_input, dir_out, filename,nargout=0) 

            print('File ' + temp_output + ' processed')
            
            ### END OF THE LOOP ###

print('done')

