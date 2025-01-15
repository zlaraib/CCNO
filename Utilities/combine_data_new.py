# reads all the data in "run*/datafiles/t_*" and combines them into a single file in "combined_data/t_*".
# run with: python3 combine_data_new.py

import numpy as np
import os

# list directories in sorted list for all directories in the current directory starting with "run"
dirnames = [d for d in os.listdir() if (os.path.isdir(d) and d.startswith("run"))]
dirnames.sort()

# create a list of filenames based on the files in the first dirname that start with "t_"
filenames = [f for f in os.listdir(dirnames[0]+"/datafiles") if f.startswith("t_")]

# create a new directory called "combined_data"
os.makedirs("combined_data", exist_ok=True)

for f in filenames:
    data = []
    for d in dirnames:
        data.append(np.genfromtxt(d+"/datafiles/"+f))

    nvalues = data[0].shape[1]

    combined_data = np.concatenate(data, axis=0)

    # remove regions of overlap based on the time value in the first column
    # remove entries with a time earlier than the largest time up to that point in the array
    # remove only individual elements, not the entire rest of the array
    i = 1
    while i < combined_data.shape[0]:
        if combined_data[i, 0] < combined_data[i - 1, 0]:
            combined_data = np.delete(combined_data, i, axis=0)
        else:
            i += 1

    print(f, combined_data.shape)
    
    # write the combined data file to a new directory called "combined_data"
    np.savetxt("combined_data/"+f, combined_data)

