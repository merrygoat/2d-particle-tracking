import numpy as np
from sys import argv

filename = argv[1]

outputfile = open("xyztrajectory.xyz", 'w')

data = np.loadtxt(filename)
numframes = np.amax(data[:, 5])

for frame in range(0, int(numframes) + 1):
    thisframe = data[np.where(data[:, 5] == frame), :]
    thisframe = thisframe[0]  # The where statement nests the data in another layer of array so unnnest it.
    outputfile.write(str(thisframe.shape[0]) + "\n")
    outputfile.write("Frame number " + str(frame) + "\n")
    for i in range(0, thisframe.shape[0]):
        if thisframe[i, 7] == 0:
            outputfile.write("A\t")
        else:
            outputfile.write("B\t")
        outputfile.write(str(thisframe[i, 0]) + "\t" + str(thisframe[i, 1]) + "\t" + str(thisframe[i, 2]) + "\n")

outputfile.close()