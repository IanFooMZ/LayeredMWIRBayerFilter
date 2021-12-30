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

# with open('lastEpIter.txt','w') as outFile:
#     outFile.write("Last Epoch: %d, Last Iteration: %d" % (lastEpoch, lastIter))
# print("Last Epoch: %d, Last Iteration: %d" % (lastEpoch, lastIter))

# Plot FOM trace
import matplotlib.pyplot as plt
trace = []
numEpochs = len(f); numIter = len(f[0])
upperRange = np.ceil(np.max(f))
fig = plt.figure(figsize=(6.4,3))
vlinestyle = {'color': 'gray', 'linestyle': '--', 'linewidth': 1}
for i in range(len(f)):
    trace = np.concatenate((trace,f[i]),axis=0)
    plt.vlines((i+1),0,upperRange,**vlinestyle)
plt.plot(np.linspace(0,10,len(trace)),trace)
# plt.plot(trace)
plt.vlines(0,0,upperRange,**vlinestyle)
plt.vlines(numEpochs,0,upperRange,**vlinestyle)
plt.title('Figure of Merit - Trace')
plt.xlabel('Epoch'); plt.ylabel('Figure of Merit')
plt.savefig('fom_trace.png')