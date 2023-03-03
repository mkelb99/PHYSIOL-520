#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 13:47:15 2022

@author: Brian
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp

# Function of the set of ODEs for Beeler-Reuter ventricular CM simulation
def dXdT(t, X):

    # State variables
    V, m, h, j, Cai, d, f, x1 = X   # V --> Membrane potential (mV)    
                                    # m --> Na channel m gate prob (uls)
                                    # h --> Na channel h gate prob (uls)
                                    # j --> Na channel j gate prob (uls)
                                    # Cai --> Cytosolic Ca2+ conc (uM??)
                                    # d --> Slow inwrd ch d gate prob (uls)
                                    # f --> Slow inwrd ch f gate prob (uls)
                                    # x1 --> Time-dep out K x1 gt prob (uls)
                                    
    # Model parameters                              
    Cm = 1                          # Tot membrne capactnce (uF/cm^2)
    g_Na = 4                        # Max Na chnl conduct (mS/cm^2)
    g_Nac = 0.003                   # Ohmic Na chnl cond (mS/cm^2)
    E_Na = 40                       # Na channel Nernst potentl (mV) 
    g_s = 0.09                      # Slw inwrd cur conduct (mS/cm^2)
    # Simulation parameters
    t_SimStart = 10                 # When the stim cycle starts (ms)
    Period = 1000                   # Period of the stim current (ms)
    t_Stim = 1                      # Stimulation duration (ms)
    I_StimAmp = 50                  # Amplitude of stim current (uA/cm^2)
    
    # STIM CURRENT
    t_Cycle = t%Period
    if (t >= t_SimStart and t_Cycle >= 0 and t_Cycle <= t_Stim):
        I_Stim = I_StimAmp
    else:
        I_Stim = 0
        
    # SODIUM CHANNEL
	# Sodium channel current
    i_Na = ((g_Na * m**3 * h * j) + g_Nac) *\
        (V - E_Na)                              # Na channel current (eqn 5)            
	# Sodium channel m gate
    alpha_m = ((-1) * (V + 47)) /\
        (np.exp(-0.1 * (V + 47)) - 1)           # Na chnl m gt opn rt (eqn 13)
    beta_m = 40 * np.exp(-0.056 * (V + 72))     # Na ch m gte cls rate (eqn 13)
    dmdt = (alpha_m * (1-m)) - (beta_m * m)     # RoC Na ch m gate (eqn 10-12)
    # Sodium channel h gate
    alpha_h = 0.126 * np.exp(-0.25 * (V + 77))  # Na chnl h gt opn rt (eqn 13)
    beta_h = 1.7 / (np.exp(-0.082 *\
        (V + 22.5)) + 1)                        # Na chnl h gt cls rt (eqn 13)
    dhdt = (alpha_h * (1-h)) - (beta_h * h)     # RoC Na ch h gate (eqn 10-12)
    # Sodium channel j gate
    alpha_j = (0.055 * np.exp(-0.25 *\
        (V + 78))) / (np.exp(-0.2 *\
        (V + 78)) + 1)                          # Na chnl j gt opn rt (eqn 13)
    beta_j = 0.3 / (np.exp(-0.1 * (V+32)) + 1)  # Na chnl j gte cls rt (eqn 13)
    djdt = (alpha_j*(1-j)) - (beta_j*j)         # RoC Na ch j gate (eqn 10-12)

    # SLOW INWARD CALCIUM CHANNEL
    # Slow inward channel Nernst, current and Ca handling
    E_s = -82.3 - (13.0287 * np.log(Cai))       # Slw inwrd Nernst ptn (eqn 7)
    i_s = g_s * d * f * (V - E_s)               # Slow inwrd Ca current (eqn 6)
    dCaidt = (-0.0000001 * i_s) +\
        (0.07 * (0.0000001 - Cai));             # RoC cytosolic Ca (eqn 9)
    # Slow inward channel d gate
    alpha_d = (0.095 * np.exp(-0.01 * (V - 5))) /\
        (1 + np.exp(-0.072 * (V - 5)))          # Slwin ch d gt opn rt (eqn 13)
    beta_d = (0.07 * np.exp(-0.017 * (V + 44))) /\
        (1 + np.exp(0.05 * (V + 44)))           # Slwin ch d gt cls rt (eqn 13)
    dddt = (alpha_d * (1 - d)) - (beta_d * d)   # RoC slw in d gate (eqn 10-12)
    # Slow inward channel f gate
    alpha_f = (0.012 * np.exp(-0.008 * (V+28))) /\
        (1 + np.exp(0.15 * (V + 28)));          # Slwin ch f gt opn rt (eqn 13)
    beta_f = (0.0065 * np.exp(-0.02 * (V + 30))) /\
        (1 + np.exp(-0.20 * (V + 30)))          # Slwin ch f gt cls rt (eqn 13)
    dfdt = (alpha_f * (1-f)) - (beta_f * f)     # RoC slw in f gate (eqn 10-12)

    # TIME-DEPENDENT OUTWARD POTASSIUM CHANNEL
    i_x1 = x1 * (0.8 * (np.exp(0.04 * (V + 77)) - 1)) /\
        np.exp(0.04 * (V + 35))                 # Time-dep outwrd K cur (eqn 4)
    # Time-dependent outward potassium x1 gate
    alpha_x1 = (0.0005 * np.exp(0.083 * (V + 50))) /\
         (np.exp(0.057 * (V + 50)) + 1)         # Td outK x1 gt opn rt (eqn 13)
    beta_x1 = (0.0013 * np.exp(-0.06 * (V + 20))) /\
        (np.exp(-0.04 * (V + 20)) + 1)          # Td outK x1 gt cls rt (eqn 13)
    dx1dt = (alpha_x1 * (1 - x1)) - \
        (beta_x1 * x1)                          # RoC Td outK x1 gt (eqn 10-12)
    
    # TIME-INDEPENDENT OUTWARD POTASSIUM CHANNEL
    i_K1 = 0.35 * ((4 * (np.exp(0.04 * (V+85)) - 1)) /\
        (np.exp(0.08 * (V+53)) + np.exp(0.04 * (V+53))) +\
        (0.2 * (V+23)) / (1 -\
        np.exp(-0.04 * (V+23))))                # Tm-indp outwrd K cur (eqn 2)
    

    # MEMBRANE POTENTIAL
    dVdt = (I_Stim -\
        (i_Na + i_s + i_x1 + i_K1)) / Cm;       # RoC of membrn potentl (eqn 8)
        
    RoC = dVdt, dmdt, dhdt, djdt, dCaidt, dddt, dfdt, dx1dt
    return RoC

# Set evaluation time points and then integrate
tspan = np.arange(0., 3000., 0.1)
results = solve_ivp(dXdT, [0., 3000],\
                    [-80, 0.1, 0.9, 0.9, 0.0001, 0.05, 0.9, 0.05],\
                    method = 'Radau', atol = 1e-9,\
                    max_step = 1, t_eval = tspan)
Cai0 = 0.0001;
d0 = 0.003;
f0 = 0.994;
x10 = 0.0001;
V0 = -84.624;
plt.figure(1)    
plt.plot(results.t/1000,results.y[0])
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (mv)')