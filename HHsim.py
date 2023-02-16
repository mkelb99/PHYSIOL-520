# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 09:39:20 2023

@author: 14073

instructions are to play around/change something, explain why it makes sense
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def dXdT(t,X):
    #state variables
    v, m, h, n = X
    ''' Or
    v=X[0]
    m=X[1]
    h=X[2]
    n=X[3] '''
    #resting potential, conductivity, capacitance, applied current
    vNa=115   #nernst relative to rest
    vK=-12
    vL= 10.6
    gNa=120
    gK=36
    gL=.3
    Iapp=3
    Cm=1e-6
    #alphas and betas
    alphan=.01*(10-v)/(np.exp((10-v)/10)-1)
    betan=.125*(np.exp(-v/80))
    alpham=.1*(25-v)/(np.exp((25-v)/10)-1)
    betam=4*np.exp(-v/18)
    alphah=.07*np.exp(-v/20)
    betah=1/(np.exp((30-v)/10)-1)
    #computing currents
    Ina=(m**3)*h*gNa*(v-vNa)
    Ik=(n**4)*gK*(v-vK)
    Il=gL*(v-vL)
    #computing derivatives
    dv= (-Ina - Ik - Il + Iapp)/Cm
    dm= alpham*(1-m)-betam*m
    dh=alphah*(1-h)-betah*h
    dn=alphan*(1-n)-betan*n
    dX= [dv, dm, dh, dn]
    return dX

def dXdT2(t,X):
    #state variables
    v, m, h, n = X
    ''' Or
    v=X[0]
    m=X[1]
    h=X[2]
    n=X[3] '''
    #resting potential, conductivity, capacitance, applied current
    vNa=115   #nernst relative to rest
    vK=-12
    vL= 10.6
    gNa=120
    gK=36
    gL=.3
    Iapp=6
    Cm=1e-6
    #alphas and betas
    alphan=.01*(10-v)/(np.exp((10-v)/10)-1)
    betan=.125*(np.exp(-v/80))
    alpham=.1*(25-v)/(np.exp((25-v)/10)-1)
    betam=4*np.exp(-v/18)
    alphah=.07*np.exp(-v/20)
    betah=1/(np.exp((30-v)/10)-1)
    #computing currents
    Ina=(m**3)*h*gNa*(v-vNa)
    Ik=(n**4)*gK*(v-vK)
    Il=gL*(v-vL)
    #computing derivatives
    dv= (-Ina - Ik - Il + Iapp)/Cm
    dm= alpham*(1-m)-betam*m
    dh=alphah*(1-h)-betah*h
    dn=alphan*(1-n)-betan*n
    dX2= [dv, dm, dh, dn]
    return dX2


tspan = np.arange(0,40,.1)
tspan2 = np.arange(41,81,.1)
results= solve_ivp(dXdT, [0, 40], [0,0,.75,0],method='LSODA', t_eval=tspan)
results2 =solve_ivp(dXdT2, [41, 81], [results.y[0,399],results.y[1,399],results.y[2,399],results.y[3,399]],method='LSODA', t_eval=tspan2)

plt.plot(np.append(results.t,results2.t),np.append(results.y[0],results2.y[0]))
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.show() 

