import os

inputpath = '../data/'
outputpath = 'output/extracted/'

# This piece of code lists the filepath of .wav files in an architechture
for dirpath, dirnames, filenames in os.walk(inputpath):
    for filename in filenames:
        if filename.endswith('.wav'):
            # Creates filepath for input file, and output file
            temp_input  = dirpath + filename
            temp_output = temp_input.replace('../data/','output/extracted/')
            
            # Applies matlab function to extract
            # 
            # Matlab function should take temp_input and temp_output. Matlab
            # will deal with every file input/output

