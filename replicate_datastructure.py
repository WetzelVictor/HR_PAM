import os

input_path = 'data/repetabilite/'
output_path = 'output/repetabilite/'

for dirpath, dirnames, filenames in os.walk(input_path):
    if not '%' in dirpath:
        structure = os.path.join(output_path, dirpath[len(input_path):])
    if not os.path.isdir(structure):
        os.mkdir(structure)
    else:
        print("Folder does already exits!")
