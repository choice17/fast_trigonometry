#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pylab as plt
import sys

R = int(sys.argv[2])
N = float(sys.argv[1])
dat  = open("diff.bin", "rb").read()
dat = np.fromstring(dat, np.float32)
dat = dat.reshape(-1,4).T
F,M = dat.shape
logdiff = np.log2(dat)
xaxis = np.arange(R) * 3.14 * N / R - 3.14 * N / 2
fig = plt.figure()
ax = [0,0,0,0]
name = ["atan", "sin", "cos", "tan"]
for i in range(4):
    ax[i]= fig.add_subplot(4,1,i+1)
    ax[i].plot(xaxis, dat[i])
    ax[i].set_ylabel(name[i])
ax[3].set_xlabel('angle(radian)')
plt.show()

