# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:59:58 2023
"""


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
# function that returns dy/dt
def model(psi,t):
    NernstNa=55
    Nernstk=-80
    C=1.0e-6
    gna=0e-3
    gk=0e-3
    if t<.002:
        gk=0e-3
        gna=0e-3
    elif t>=.002 and t<0.003:
        gk=0e-3
        gna=4e-3
    elif t>=.003:
        gk=1e-3
        gna=.1e-3
    dpsidt = (-gna*(psi-NernstNa)-(gk*(psi-Nernstk)))/C
    return dpsidt

# initial condition
psi0 = -70

# time points
t = np.linspace(0,.01,1000)

# solve ODE
psi = odeint(model,psi0,t)

# plot results
plt.plot(t,psi)
plt.xlabel('time')
plt.ylabel('dpsi(t)')
plt.show() """

#Try using similar concept to LIF model
def model(psi,t):
    NernstNa=55
    Nernstk=-80
    C=1.0e-6
    gna=0e-3
    gk=0e-3
    if t<.002:
        gk=0e-3
        gna=0e-3
    elif t>=.002 and t<0.003:
        gk=0e-3
        gna=4e-3
    elif t>=.003 and t<.01:
        gk=1e-3
        gna=.1e-3
    elif t>=.01 and t<.03:
        dpsidt=-80
    elif t>=.03 and t<.031:
        gk=0e-3
        gna=4e-3
    elif t>=.031:
        gk=1e-3
        gna=.1e-3 
    dpsidt = (-gna*(psi-NernstNa)-(gk*(psi-Nernstk)))/C
    return dpsidt

# initial condition
psi0 = -70

# time points
t = np.linspace(0,.1,1000)

# solve ODE
psi = odeint(model,psi0,t)

# plot results
plt.plot(t,psi)
plt.xlabel('time')
plt.ylabel('dpsi(t)')
plt.show() 