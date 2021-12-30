import numpy as np

f = np.load("figure_of_merit.npy")
lastEpochReached = False

idx = []
for i in range(len(f)):
    if f[i][1] == 0:
        lastEpochReached = True
        break
    else:
        idx.append(np.max(np.nonzero(f[i])))

lastEpoch = len(idx)-1
lastIter = idx[lastEpoch]

with open('lastEpIter.txt','w') as outFile:
    outFile.write("Last Epoch: %d, Last Iteration: %d" % (lastEpoch, lastIter))
print("Last Epoch: %d, Last Iteration: %d" % (lastEpoch, lastIter))