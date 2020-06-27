#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pylab as plt
import sys

R = int(sys.argv[2])
N = float(sys.argv[1])
dat  = open("diff.bin", "rb").read()
dat = np.frombuffer(dat, np.float32)
name = ["atan", "sin", "cos", "tan", "acos", "asin"]
dat = dat.reshape(-1,len(name)).T
F,M = dat.shape
logdiff = np.log2(dat)
xaxis = np.arange(R) * N / R - N / 2
fig = plt.figure()
ax = [0,0,0,0,0,0]
MIN = -0.01
MAX = 0.01
for i in range(len(name)):
    ax[i]= fig.add_subplot(len(name),1,i+1)
    ax[i].plot(xaxis, dat[i])
    ax[i].set_ylabel(name[i])
    ax[i].axis(
        ymin=np.maximum(MIN,np.min(dat[i])),
        ymax=np.minimum(MAX, np.max(dat[i]))
        )
ax[-1].set_xlabel('angle(radian)')
plt.show()

