#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 20:26:04 2022

@author: Brian
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp

# Function that represents set of ODEs for Noble Pacemaker simulation
def dXdT(t, X):
    
    # State Variables:
    V, m, h, n, d, f = X                # V --> Membrane potential (mV)
                                        # m --> Na chnl m gate opn prob (uls)
                                        # h --> Na chnl h gate opn prob (uls)
                                        # n --> K chnl n gate opn prob (uls)
                                        # d --> Ca chnl n gate opn prob (uls)
                                        # f --> Ca chnl f gate opn prob (uls)
                                        
    F = 96.4867                         # Faraday constant (C/mmole)
    R = 8.3143                          # Universal gas constant (J/(mole*K))
    T = 310                             # Temperature (K)
    z_Ca = 2                            # Calcium valance
    VFoRT = (V*F)/(R*T)                 # Calculated term to use later
    
    g_Na_max = 400                            # Max Na channel cond (mS/cm^2)
    E_Na = 40                                 # Na channel Nernst potentl (mV) 
    E_K = -100                                # K channel Nernst potentl (mV)
    P_CaL = 0.00054                           # L-type Ca chnl perm (cm/s)
    KmCa = 0.6                                # L-type Ca fCa gt half act (uM)
    g_L = 0.075                               # Max leak current cond (mS/cm^2)
    E_L = -60                                 # Leak chnl Nernst potentl (mV)
    C_m = 12                                  # Membrane capacitance (uF/cm^2)
    Ca_i = 0.70                               # Cytosolic Ca conc (uM)
    Ca_o = 1800                               # Extracellular Ca conc (uM)
    
    t_SimStart = 10                           # When the stim cycle starts (ms)
    Period = 1000                             # Stim frequency period (ms)
    Deltat = 10;                               # Stim duration (ms)
    Iapp_Amp = 50;                            # Applied current amp (uA/cm^2)
    
    t_Cycle = t%Period
    if (t >= t_SimStart and t_Cycle >= 0 and t_Cycle <= Deltat):
        I_App = Iapp_Amp
    else:
        I_App = 0

    # SODIUM CHANNEL
	# Sodium channel current
    g_Na = m**3 * h * g_Na_max                # Na channel var cond (eqn 12)
    # i_Na = (g_Na + 0.140) * (V - E_Na)         # Na channel current (eqn 20)
    i_Na = g_Na * (V - E_Na)         # Na channel current
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
    g_K1 = (1.2 * np.exp(((-1*V) - 90) / 50)) +\
        (0.015 * np.exp((V + 90) / 60))         # K chnl inwrd rct cond (eqn 5)
    g_K2 = 1.2 * n**4                           # K ch outwrd rect cond (eqn 6)
    i_K = ((g_K1)+g_K2) * (V - E_K)     #modified to make 35% reduction in Ik1          # K channel current (eqn 10)
	# Potassium channel n gate (Note rates are in 1/ms)
    alpha_n = (0.0001 * ((-1*V) - 50)) /\
        (np.exp(((-1*V) - 50) / 10) - 1)        # K chnl n gt opn rate (eqn 8)
    beta_n = 0.002 * np.exp(((-1*V) - 90) / 80) # K ch n gate open rate (eqn 9)
    dndt = (alpha_n * (1-n)) - (beta_n * n)     # RoC K channel n gate (eqn 7)

    # L-Type Calcium Channel
    # Calcium d gate
    d_inf = 1 / (1 + np.exp((-V-10)/6.24));
    tau_d = d_inf * (1 - np.exp((-V-10)/6.24)) / (0.035 * (V + 10));
    alpha_d = d_inf / tau_d;
    beta_d = (1 - d_inf)/tau_d;
    dddt = (alpha_d * (1-d)) - (beta_d * d);
    # Calcium f gate
    f_inf = (1 / (1 + np.exp((V+35.06)/8.6))) +\
        (0.6 / (1 + np.exp((-V+50)/20)));
    tau_f = 1 / ((0.0197 * np.exp(-(0.0337 * (V+10)**2))) + 0.02);
    alpha_f = f_inf / tau_f;
    beta_f = (1 - f_inf)/tau_f;
    dfdt = (alpha_f * (1-f)) - (beta_f * f);
    # Calcium fCa gate
    f_Ca = 1 / (1 + (Ca_i/KmCa)**2);
    # Calcium channel current- edited for 50% inc in calcium perm/conductance
    i_CaL = (P_CaL) * d * f * f_Ca * (z_Ca**2 * F * VFoRT) *\
        (Ca_i * np.exp(z_Ca * VFoRT) - 0.34*Ca_o) /\
        (np.exp(z_Ca * VFoRT) - 1);
	# LEAK CURRENT
    i_Leak = g_L * (V - E_L)                    # Leak current (eqn 3) 
    
    # MEMBRANE POTENTIAL
    dVdt = (-1 * (i_Na + i_K + i_CaL + i_Leak) + I_App) / C_m             
    
    RoC = dVdt, dmdt, dhdt, dndt, dddt, dfdt
    
    return RoC

# Set evaluation time points and then integrate
tspan = np.arange(0., 2000., 0.1)
results = solve_ivp(dXdT, [0., 2000], [-80, 0.1, 0.8, 0.01, 0.9, 0.9],\
        method = 'Radau', atol = 1e-9, max_step = 0.25, t_eval = tspan)
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