# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:29:29 2021

@author: lass_j
"""

try:
    import IPython
    shell = IPython.get_ipython()
    shell.enable_matplotlib(gui='qt')
except:
    pass

import numpy as np
import os.path
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit


def Gaussian(x,A,mu,sigma):
    return A*np.exp(-np.power(x-mu,2.0)/(2.0*sigma**2))

def NGaussian(x,*pars):
    if np.mod(len(pars),3)!=0:
        raise AttributeError('Number of parameters mismatch. Got '+str(len(pars)))
    return np.sum([Gaussian(x,*p) for p in np.array(pars).reshape((-1,3))],axis=0)



data = np.loadtxt(r'C:/Users/lass_j/Documents/McStasCAMEA/BackEnd/FullInstrument_v4.5_source_sample_pos_20210922_134818/ReuterStokes11_4.psd')[:,3]
fig,ax = plt.subplots()
ax.scatter(np.arange(1024),data)

mu,A = np.array([(92.53750000000008, 1069.3007009780672),
 (217.31875000000008, 1017.0551094846686),
 (334.5375000000001, 939.3927437512385),
 (451.7562500000002, 889.9712382845102),
 (571.2437500000001, 827.8413456977661),
 (686.9500000000003, 786.8920983110482),
 (810.9750000000001, 737.47059284432),
 (935.7562500000001, 682.4009153242513)]).T

guess = np.array([[a,m,10] for a,m in zip(A,mu)]).flatten()
    
fit,err = curve_fit(NGaussian, np.arange(1024), data, p0=guess)
X = np.arange(1024)

ax.plot(X,NGaussian(X,*fit))


for I,pars in enumerate(fit.reshape(-1,3)):
    start = pars[1]-3.0*pars[2]
    stop = pars[1]+3.0*pars[2]
    XX = np.linspace(start,stop)
    ax.plot(XX,Gaussian(XX,*pars),'r')
    ax.vlines([start,stop],0,1100,'r')
    print(I,int(np.floor(start)),int(np.ceil(stop)))




