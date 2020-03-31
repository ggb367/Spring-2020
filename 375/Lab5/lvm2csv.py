import os

import numpy as np

file_names = []
directory = input("What is the directory?\n")

for path, subdirs, files in os.walk(directory):
    for name in files:
        file_names.append(os.path.join(path, name))
file_renames = [file.replace('lvm', 'csv') for file in file_names]
for dex in range(np.size(file_names)):
    os.rename(file_names[dex], file_renames[dex])
    data_in = open(file_renames[dex], 'rb').readlines()
    with open(file_renames[dex], 'wb') as outfile:
        outfile.writelines(data_in[22:])
print("Conversion Complete!")