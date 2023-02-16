# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:59:58 2023
"""


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(n,t):
        alpha=.01*(10-v)/(np.exp((10-v)/10)-1)
        beta=.125*(np.exp(-v/80))
        dndt=alpha*(1-n)-beta*n
        return dndt
voltage=np.linspace(11,210,10)
for v in voltage:

    # initial condition
    n0 = .1

    # time points
    #v=100
    t = np.linspace(0,11.,100)

    # solve ODE
    n = odeint(model,n0,t)

# plot results
    plt.plot(t,n**4)
    plt.xlabel('time')
    plt.ylabel('n(t)')
plt.show() 

#plot results as conductance
#     gkmax=50
#     plot.plot(t,gkmax*(n**4))
#     plt.xlabel('time')
#     plt.ylabel('gk')
# plt.show
