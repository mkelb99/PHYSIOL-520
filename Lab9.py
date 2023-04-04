# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 09:20:07 2023

@author: 14073
"""
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp

def dXdT(t, X):
    
    ATP, ADP, Cr, CrP, Pi = X  #variables of interest
    
    #JATPASE0=.05e-3
    Vck= 1
    Km=30.0e-6
    Vm=0.1e-3
    Keq_ck=1
    
    JATPASE0 = Vm*ATP/(ATP+Km)  #variable JATPASE
    JOXPHOR=Vm*ADP/(ADP+Km)
    JCK = Vck*(CrP*ADP - Cr*ATP/Keq_ck)
    #if statement for JATPASE equaling itself or 0
    JATPASE= 5*JATPASE0/(1+(.1e-3)/ATP)
    if t< -300 or t>0.0:
        JATPASE = JATPASE0/(1+(.1e-3)/ATP)

        
    dATPdt= -JATPASE + JOXPHOR + JCK
    dADPdt= JATPASE - JOXPHOR - JCK
    dCrdt= JCK
    dCrPdt= -JCK
    dPidt= JATPASE - JOXPHOR
    
    RoC = dATPdt,dADPdt, dCrdt, dCrPdt, dPidt
    
    return RoC

 #initial condition
x0 = [10e-3, 0e-3, 40e-3, 10e-3, 2e-3]

tspan = np.arange(-400, 400, 0.1)
results = solve_ivp(dXdT, [-400, 400], x0,
                    method = 'LSODA', t_eval = tspan)

plt.plot(results.t,results.y[0],results.t,results.y[1], results.t, results.y[4])
#plt.plot(results.t,results.y[2])
plt.show()