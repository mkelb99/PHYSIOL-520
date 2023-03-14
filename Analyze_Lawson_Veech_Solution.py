# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 08:06:14 2023

@author: beardda
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import minimize


def CK_mb(x, pH, Mg_o, K_o, ATP, ADP, CrP, Pi, T):

#  INPUTS:
#  x - free ion concentrations (Mg, K)
#  pH
#  Mg_o, K_o: total ion concentrations
#  ATP, ADP, CrP, Pi: total reactant concentrations
# 
#  OUTPUT:
# E - mass balance error

  Mg, K = x

# Binding Constants
  if T == 38:
    Kh_ATP = 3.6199e-007
    Kh_ADP = 5.0335e-007
    Kh_Pi  = 2.4851e-007
    Kh_CrP = 3.53e-5
    Kmg_ATP  = 0.00014287
    Kmg_ADP  = 0.0010357
    Kmg_Pi   = 0.024585
    Kmg_CrP  = 0.0544
    Kk_ATP  = 0.12765
    Kk_ADP  = 0.16312
    Kk_Pi   = 0.4382
    Kk_CrP  = np.inf
    Kmg_HATP = 0.039516
    Kmg_HADP = 0.084579

  if T == 25:
    Kh_ATP = 3.6894e-007
    Kh_ADP = 5.149e-007
    Kh_Pi  = 2.2829e-007
    Kh_CrP = 3.27e-5
    Kmg_ATP  = 0.0001879
    Kmg_ADP  = 0.0013045
    Kmg_Pi   = 0.028442
    Kmg_CrP  = 0.058824
    Kk_ATP  = 0.12793
    Kk_ADP  = 0.16133
    Kk_Pi   = 0.43499
    Kk_CrP  = np.inf
    Kmg_HATP = 0.051933
    Kmg_HADP = 0.10162

  if T == 15:
    Kh_ATP = 3.7544e-007
    Kh_ADP = 5.2534e-007
    Kh_Pi  = 2.1292e-007
    Kh_CrP = 3.09e-5
    Kmg_ATP  = 0.00023669
    Kmg_ADP  = 0.001584
    Kmg_Pi   = 0.032155
    Kmg_CrP  = 0.062893
    Kk_ATP  = 0.12837
    Kk_ADP  = 0.16006
    Kk_Pi   = 0.4327
    Kk_CrP  = np.inf
    Kmg_HATP = 0.065326
    Kmg_HADP = 0.11854

  if T == 5:
    Kh_ATP = 3.8306e-007
    Kh_ADP = 5.3734e-007
    Kh_Pi  = 1.9774e-007
    Kh_CrP = 2.91e-5
    Kmg_ATP  = 0.00030401
    Kmg_ADP  = 0.0019546
    Kmg_Pi   = 0.036726
    Kmg_CrP  = 0.068493
    Kk_ATP  = 0.12902
    Kk_ADP  = 0.15887
    Kk_Pi   = 0.43056
    Kk_CrP  = np.inf
    Kmg_HATP = 0.083715
    Kmg_HADP = 0.14003

  H = 10**(-pH)

  P_ATP = 1.0 + H/Kh_ATP + Mg/Kmg_ATP + H*Mg/(Kmg_HATP*Kh_ATP) + K/Kk_ATP
  P_ADP = 1.0 + H/Kh_ADP + Mg/Kmg_ADP + H*Mg/(Kmg_HADP*Kh_ADP) + K/Kk_ADP
  P_CrP = 1.0 + H/Kh_CrP + Mg/Kmg_CrP + K/Kk_CrP
  P_Pi  = 1.0 + H/Kh_Pi + Mg/Kmg_Pi + K/Kk_Pi

  ATP1 = ATP/P_ATP
  ATP2 = ATP*(H/Kh_ATP)/P_ATP
  ATP3 = ATP*(Mg/Kmg_ATP)/P_ATP
  ATP4 = ATP*(H/Kh_ATP)*(Mg/Kmg_HATP)/P_ATP
  ATP5 = ATP*(K/Kk_ATP)/P_ATP

  ADP1 = ADP/P_ADP
  ADP2 = ADP*(H/Kh_ADP)/P_ADP
  ADP3 = ADP*(Mg/Kmg_ADP)/P_ADP
  ADP4 = ADP*(H/Kh_ADP)*(Mg/Kmg_HADP)/P_ADP
  ADP5 = ADP*(K/Kk_ADP)/P_ADP

  CrP1 = CrP/P_CrP
  CrP2 = CrP*(H/Kh_CrP)/P_CrP
  CrP3 = CrP*(Mg/Kmg_CrP)/P_CrP
  CrP5 = CrP*(K/Kk_CrP)/P_CrP

  Pi1 = Pi/P_Pi
  Pi2 = Pi*(H/Kh_Pi)/P_Pi
  Pi3 = Pi*(Mg/Kmg_Pi)/P_Pi
  Pi5 = Pi*(K/Kk_Pi)/P_Pi

  E = (Mg_o*1.0e6 - (Mg + ATP3 + ATP4 + ADP3 + ADP4 + CrP3 + Pi3)*1.0e6)**2 + \
        (K_o*1.0e6 - (K + ATP5 + ADP5 + CrP5 + Pi5)*1.0e6)**2
      
  return E

'''
MAIN PROGRAM:
'''

# Loading in data and assigning variables
data = np.loadtxt('Lawson_Veech_data.txt')
pH = data[:,5]
H = 10**(-pH)
Mg = data[:,4]*1e-3   # [Mg]
ATP = data[:,0]*1e-3  # [ATP]
ADP = data[:,2]*1e-3  # [ADP]
Cr  = data[:,1]*1e-3  # creatine
CrP = data[:,3]*1e-3  # cr-P
CrPi = data[:,6]*1e-3 # cr-P (init)
Mg_o = data[:,7]*1e-3 # total [Mg]
Pi   = data[:,8]*1e-3 # total [Pi]
Tris = data[:,9]*1e-3 # total [Tris]

# Dissociation constants used by Lawson and Veech
Kh_ATP   = 1.08e-007
Kh_ADP   = 1.20e-007
Kh_Pi    = 1.76e-007
Kh_CrP   = 3.16e-5
Kmg_ATP  = 7.19e-5
Kmg_ADP  = 7.58e-4
Kmg_Pi   = 0.0107
Kmg_CrP  = 0.050
Kmg_HATP = 0.0282
Kmg_HADP = 0.0309
Kh_tris  = 8.47e-9

# Equation (5.62)
Na = ATP + ADP + CrPi

# Equation (5.63)
K_pi = Pi

# Equation (5.64)
Cl_mgcl = 2*Mg_o;

# Equation (5.65)
P_ATP = 1 + H/Kh_ATP + Mg/Kmg_ATP + H*Mg/(Kmg_HATP*Kh_ATP)
P_ADP = 1 + H/Kh_ADP + Mg/Kmg_ADP + H*Mg/(Kmg_HADP*Kh_ADP)
P_CrP = 1 + H/Kh_CrP + Mg/Kmg_CrP
P_Pi  = 1 + H/Kh_Pi + Mg/Kmg_Pi
P_Tris  = 1 + H/Kh_tris

ATP1 = ATP/P_ATP
ATP2 = ATP*(H/Kh_ATP)/P_ATP
ATP3 = ATP*(Mg/Kmg_ATP)/P_ATP
ATP4 = ATP*(H/Kh_ATP)*(Mg/Kmg_HATP)/P_ATP

ADP1 = ADP/P_ADP
ADP2 = ADP*(H/Kh_ADP)/P_ADP
ADP3 = ADP*(Mg/Kmg_ADP)/P_ADP
ADP4 = ADP*(H/Kh_ADP)*(Mg/Kmg_HADP)/P_ADP

CrP1 = CrP/P_CrP
CrP2 = CrP*(H/Kh_CrP)/P_CrP
CrP3 = CrP*(Mg/Kmg_CrP)/P_CrP

Pi1 = Pi/P_Pi
Pi2 = Pi*(H/Kh_Pi)/P_Pi
Pi3 = Pi*(Mg/Kmg_Pi)/P_Pi

Tris1 = Tris/P_Tris
Tris2 = Tris*(H/Kh_tris)/P_Tris

# Equation (5.66)
H_surplus = Pi - Pi2 + 3*ATP - ATP2 - ATP4 + 3*ADP - ADP2 - ADP4  + CrP - CrP2 - Tris2
            
K_koh = np.zeros(len(H_surplus)) # initialize KOH vector
Cl_hcl = np.zeros(len(H_surplus)) # initialize HCl vector
for i in range(0,len(H_surplus)):
  if H_surplus[i] > 0:
    K_koh[i] = H_surplus[i]
  else:
    Cl_hcl[i] = -H_surplus[i]

I = 0.5*( Na + K_pi + Cl_mgcl + K_koh + Cl_hcl + \
    16*ATP1 + 9*ATP2 + 4*ATP3 + ATP4 + \
    9*ADP1 + 4*ADP2 + ADP3 + \
    4*CrP1 + CrP2 + \
    4*Pi1 + Pi2 + \
    Tris2)

K_kcl = 0.25 - I

'''
Now do new calculations acounting for K 
'''

K_o = K_pi + K_koh + K_kcl # total K 
Mg = np.zeros(len(H_surplus))
K = np.zeros(len(H_surplus))
for i in range(0,len(H_surplus)):
  bnds = ((0.0,1.0),(0.0,1.0))
  solution = minimize(CK_mb,[Mg_o[i]/10  ,    0.20],\
                    args=(pH[i],Mg_o[i],K_o[i],ATP[i],ADP[i],CrP[i],Pi[i],38.0),\
                    method='TNC',bounds=bnds)
  Mg[i] = solution.x[0]
  K[i] = solution.x[1]
    
# Calculate new binding polynomials based on estimaged [Mg] and [K]
Kh_ATP = 3.6199e-007
Kh_ADP = 5.0335e-007
Kh_Pi  = 2.4851e-007
Kh_CrP = 3.53e-5
Kmg_ATP  = 0.00014287
Kmg_ADP  = 0.0010357
Kmg_Pi   = 0.024585
Kmg_CrP  = 0.0544
Kk_ATP  = 0.12765
Kk_ADP  = 0.16312
Kk_Pi   = 0.4382
Kk_CrP  = 0.48978
Kmg_HATP = 0.039516
Kmg_HADP = 0.084579

P_ATP = 1 + H/Kh_ATP + Mg/Kmg_ATP + H*Mg/(Kmg_HATP*Kh_ATP)
P_ADP = 1 + H/Kh_ADP + Mg/Kmg_ADP + H*Mg/(Kmg_HADP*Kh_ADP)
P_CrP = 1 + H/Kh_CrP + Mg/Kmg_CrP
P_Pi  = 1 + H/Kh_Pi + Mg/Kmg_Pi
P_Tris  = 1 + H/Kh_tris

Kobs = ATP*Cr / (ADP*CrP) # observed Keq from data
Mgf = 10**(-np.linspace(1.5,6.8,25)) # range of [Mg] to plot

P_ATP = 1 + np.mean(H)/Kh_ATP + Mgf/Kmg_ATP + np.mean(H)*Mgf/(Kmg_HATP*Kh_ATP) + np.mean(K)/Kk_ATP
P_ADP = 1 + np.mean(H)/Kh_ADP + Mgf/Kmg_ADP + np.mean(H)*Mgf/(Kmg_HADP*Kh_ADP) + np.mean(K)/Kk_ADP
P_CrP = 1 + np.mean(H)/Kh_CrP + Mgf/Kmg_CrP + np.mean(K)/Kk_CrP
Keq = 3.0e8

plt.plot(-np.log10(Mgf),Keq*P_ATP*P_CrP/P_ADP)
plt.plot(-np.log10(Mg),Kobs/H,'o');
plt.xlabel('pMg = -log_{10}[Mg^{2+}]')
plt.ylabel('K_{obs}/[H^+] (M^{-1})')
plt.show()

