#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 20:26:04 2022

@author: Brian

edited 2/26/2023 by Madison Kelberman
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp

# Function that represents set of ODEs for Noble Pacemaker simulation
def dXdT(t, X):
    
    # State Variables:
    V, m, h, n, d, f = X                      # V --> Membrane potential (mV)
                                              # m --> Na chnl m gate opn prob (uls)
                                              # h --> Na chnl h gate opn prob (uls)
                                              # n --> K chnl n gate opn prob (uls)
                                              # d --> Ca chnl d gate opn prob
                                              # f --> Ca chnl f gate opn prob

    F = 96.4867                               # Faradays (C/mmol)
    R = 8.3143                                # Gas Constant (J/mol*K)
    T = 310                                   # temp (Kelvin)
    z_ca = 2                                  # valence of calcium
    VFRT = (V*F)/(R*T)                        # shortcut term

    gNamax= 400 *.9                           #max sodium conductance from Noble adjusted 
    Ena = 40                                  #noble na nernst
    Ek = -100                                 #noble k nernst
    gK1factor= .65   
    P_CaL= .00054 * 15                        #adjusted permeability of calcium
    KmCa = .6                                                                                          
    gL= .075
    EL= -60
    Cm= 12
    Cai=.7
    Cao= 1800
    
    t_simstart= 10                  #start of stim cycle
    period=1000                     #stim freq period (ms)
    Deltat = 10                     #stim duration (ms)
    Iappamp = 50                       #applied current amplitude (uA/cm^2)
    
    t_cycle = t%period
    if (t >= t_simstart and t_cycle >=0 and t_cycle <=Deltat ):
        Iapp=Iappamp
    else:
        Iapp=0

    # SODIUM CHANNEL
	# Sodium channel current
    g_Na = m**3 * h * gNamax                 # Na channel var cond (eqn 12)
    i_Na = (g_Na) * (V - Ena)         # Na channel current (eqn 20)
	# Sodium channel m gate (Note rates are in 1/ms)
    alpha_m = (0.100 *((-1*V) - 48)) /\
        (np.exp(((-1*V) - 48) / 15) - 1)       # Na chnl m gt opn rt (eqn 18)
    beta_m = (0.120 * (V + 8)) /\
        (np.exp((V + 8) / 5) - 1)              # Na chnl m gt cls rt (eqn 19)
    dmdt = (alpha_m * (1-m)) - (beta_m * m)    # RoC Na channel m gate (eqn 13)
    # Sodium channel h gate (Note rates are in 1/ms)
    alpha_h = 0.170 * np.exp(((-1*V) - 90) / 20)  # Na ch h gt opn rt (eqn 16)
    beta_h = 1 / (np.exp(((-1*V) - 42) / 10) + 1) # Na ch h gt cls rt (eqn 17)
    dhdt = (alpha_h * (1-h)) - (beta_h * h)       # RoC Na chnl h gate (eqn 14)

    # POTASSIUM CHANNEL
	# Potassium channel current
    g_K1 = gK1factor * (1.2 * np.exp(((-1*V) - 90) / 50)) +\
        (0.015 * np.exp((V + 90) / 60))         # K chnl inwrd rct cond (eqn 5)
    g_K2 = 1.2 * n**4                           # K ch outwrd rect cond (eqn 6)
    i_K = (g_K1+g_K2) * (V - Ek)               # K channel current (eqn 10)
	# Potassium channel n gate (Note rates are in 1/ms)
    alpha_n = (0.0001 * ((-1*V) - 50)) /\
        (np.exp(((-1*V) - 50) / 10) - 1)        # K chnl n gt opn rate (eqn 8)
    beta_n = 0.002 * np.exp(((-1*V) - 90) / 80) # K ch n gate open rate (eqn 9)
    dndt = (alpha_n * (1-n)) - (beta_n * n)     # RoC K channel n gate (eqn 7)

    #L-type calcium channel
    #d gate
    d_inf= 1/(1 + np.exp((-V-10)/6.24))
    tau_d = d_inf * (1-np.exp((-V-10)/6.24)) / (.035 * (V+10) )
    alpha_d = d_inf/tau_d
    beta_d = (1 - d_inf)/tau_d
    dddt = (alpha_d * (1-d)) - (beta_d * d)
    #f gate
    f_inf= (1/(1 + np.exp((V+35.06)/8.6))) +\
        (.6 / (1 + np.exp((-V+50)/20)))
    tau_f = 1 / ((.0197 * np.exp(-(.0337 * (V+10)**2))) + .02)
    alpha_f = f_inf / tau_f
    beta_f = (1-f_inf)/tau_f
    dfdt = (alpha_f * (1-f))- (beta_f*f)
    #fCa gate
    f_ca= 1 / (1 + (Cai/KmCa)**2)
    #calcium current
    i_CaL = P_CaL* d * f * f_ca * (z_ca**2 * F * VFRT) *\
        (Cai * np.exp(z_ca * VFRT) - .34*Cao) /\
        (np.exp(z_ca *VFRT)-1)
    
	# LEAK CURRENT
    i_Leak = gL * (V - EL)                    # Leak current (eqn 3) 
    
    # MEMBRANE POTENTIAL
    dVdt =( -1 * (i_Na + i_K + i_CaL + i_Leak) + Iapp) /Cm
    
    RoC = dVdt, dmdt, dhdt, dndt, dddt, dfdt
    
    return RoC

# Set evaluation time points and then integrate
tspan = np.arange(0., 2000, 0.1)
results = solve_ivp(dXdT, [0., 2000], [-80, 0.1, 0.8, 0.01, .9, .5],
                    method = 'Radau', atol = 1e-9, t_eval = tspan)
plt.figure(1)    
plt.plot(results.t/1000,results.y[0])
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (mv)')

plt.figure(2)
plt.subplot(2,1,1)
plt.plot(results.t/1000,results.y[1],'-r',label='m gate')
plt.plot(results.t/1000,results.y[2],'-g',label='h gate')
plt.ylabel('m and h gate open prob')
plt.legend(loc='upper right')
plt.subplot(2,1,2)
plt.plot(results.t/1000,results.y[3],'-k',label='n gate')
plt.xlabel('Time (s)')
plt.ylabel('n gate open prob')
plt.legend(loc='lower right')

plt.figure(3)
plt.plot(results.t/1000,results.y[4],'-r',label='d gate')
plt.plot(results.t/1000,results.y[5],'-g',label='f gate')
plt.ylabel('d and f gate open prob')
plt.legend(loc='upper right')