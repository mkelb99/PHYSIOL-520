import numpy as np 
import matplotlib.pyplot as plt 

Vm_Grid = 100

V_m = np.linspace(-120,50,Vm_Grid)

HHalpha_m = np.zeros((Vm_Grid,1))
HHbeta_m = np.zeros((Vm_Grid,1))
HHalpha_h = np.zeros((Vm_Grid,1))
HHbeta_h = np.zeros((Vm_Grid,1))
HHalpha_n = np.zeros((Vm_Grid,1))
HHbeta_n = np.zeros((Vm_Grid,1))

DNalpha_m = np.zeros((Vm_Grid,1))
DNbeta_m = np.zeros((Vm_Grid,1))
DNalpha_h = np.zeros((Vm_Grid,1))
DNbeta_h = np.zeros((Vm_Grid,1))
DNalpha_n = np.zeros((Vm_Grid,1))
DNbeta_n = np.zeros((Vm_Grid,1))

for i in range(len(V_m)):
    
    HHV_m = V_m[i] + 56
    # HH Sodium Opening and Closing Rates for m and h Gates
    HHalpha_m[i] = 0.10 * (25 - HHV_m) / (np.exp((25 - HHV_m)/25) - 1)
    HHbeta_m[i] = 4.00 * np.exp(-HHV_m/18)
    HHalpha_h[i] = 0.07 * np.exp(-HHV_m/20)
    HHbeta_h[i] = 1.00 / (np.exp((30 - HHV_m)/30) + 1)
    # HH Potassium Opening and Closing Rates for n Gate
    HHalpha_n[i] = 0.01 * (10 - HHV_m) / (np.exp((10 - HHV_m)/10) - 1)
    HHbeta_n[i] = 0.125 * np.exp(-HHV_m/80)
    
    DNV_m = V_m[i]
    # DN Sodium Opening and Closing Rates for m and h Gates
    DNalpha_m[i] = 0.10 * (-DNV_m - 48) / (np.exp((-DNV_m - 48)/15) - 1)
    DNbeta_m[i] = 0.12 * (DNV_m + 8) / (np.exp((DNV_m + 8)/5) - 1)
    DNalpha_h[i] = 0.17 * np.exp((-DNV_m -90)/20)
    DNbeta_h[i] = 1.00 / (np.exp((-DNV_m - 42)/10) + 1)
    # DN Potassium Opening and Closing Rates for n Gate
    DNalpha_n[i] = 0.0001 * (-DNV_m - 50) / (np.exp((-DNV_m - 50)/10) - 1)
    DNbeta_n[i] = 0.002 * np.exp((-DNV_m - 90)/80)
    
plt.figure(1) 
plt.subplot(2,1,1)       
plt.plot(V_m,HHalpha_m,'-g',label='HH alpha_m')
plt.plot(V_m,HHbeta_m,'--g',label='HH beta_m')
plt.ylabel('$\\alpha_m$ and $\\beta_m$ (1/ms)')
plt.ylim([0,15])
plt.legend(loc='upper right')
plt.subplot(2,1,2)    
plt.plot(V_m,DNalpha_m,'-b',label='Noble alpha_m')
plt.plot(V_m,DNbeta_m,'--b',label='Noble beta_m')
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('$\\alpha_m$ and $\\beta_m$ (1/ms)')
plt.ylim([0,15])
plt.legend(loc='upper right')

plt.figure(2)
plt.subplot(2,1,1)     
plt.plot(V_m,HHalpha_h,'-g',label='HH alpha_h')
plt.plot(V_m,HHbeta_h,'--g',label='HH beta_h')
plt.ylabel('$\\alpha_h$ and $\\beta_h$ (1/ms)')
plt.ylim([0,2])
plt.legend(loc='upper right')
plt.subplot(2,1,2)   
plt.plot(V_m,DNalpha_h,'-b',label='Noble alpha_h')
plt.plot(V_m,DNbeta_h,'--b',label='Noble beta_h')
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('$\\alpha_h$ and $\\beta_h$ (1/ms)')
plt.ylim([0,2])
plt.legend(loc='upper right')

plt.figure(3) 
plt.subplot(2,1,1)   
plt.plot(V_m,HHalpha_n,'-g',label='HH alpha_n')
plt.plot(V_m,HHbeta_n,'--g',label='HH beta_n')
plt.ylabel('$\\alpha_n$ and $\\beta_n$ (1/ms)')
plt.legend(loc='upper right')
plt.subplot(2,1,2)
plt.plot(V_m,DNalpha_n,'-b',label='Noble alpha_n')
plt.plot(V_m,DNbeta_n,'--b',label='Noble beta_n')
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('$\\alpha_n$ and $\\beta_n$ (1/ms)')
plt.legend(loc='upper right')