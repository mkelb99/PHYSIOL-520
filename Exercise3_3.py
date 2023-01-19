# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 14:17:44 2022

@author: beardda
"""
# Example3_4.py
import numpy as np
import matplotlib.pyplot as plt

#x1 = np.array([-0.090,-0.039,0,0.036,0.092,0.136,0.183])
#y1 = np.array([0.2810,0.5965,0.7184,0.6138,0.2939,0.0778,0.0259])
#plt.plot(x1,y1,'s')

def my_Gaussian(x):
    return M*np.exp(-x*x/(4*D*t))/(np.sqrt(4*np.pi*D*t))

D = 0.32e-6
M = 0.115
#t = 3000

#plots for various values of time (part A)
#times = np.linspace(600,6000,25)

#X = np.linspace(-0.35,0.35,500)

##make an empty list
#peaklist=[]

# for t in times:
#   Y = my_Gaussian(X)
#   peak=Y.max()
#   peaklist.append(peak)
#   plt.subplot(1,2,1)
#   plt.plot(X,Y)
#   plt.xlabel('x (cm)')
#   plt.ylabel('Concentration')

# plt.subplot(1,2,2)
# plt.plot(times,peaklist)
# plt.xlabel('Time (Seconds)')
# plt.ylabel('Peak Concentration Value')
# plt.show()


#Part B - measure of spread by eye
times=np.array([600,6000,10000])
X=np.linspace(-0.35,0.35,500)
for t in times:
    Y = my_Gaussian(X)
    plt.subplot(1,2,1)
    plt.plot(X,Y)
    plt.xlabel('x (cm)')
    plt.ylabel('Concentration')
    ##increase x resolution to better estimate
    #plt.xticks(np.linspace(-.4,.4,10))
    
#estimating between midpoint of max in min to get inflection based on visuals
inflection=np.array([.03*2,.08*2,.13*2])
plt.subplot(1,2,2)
plt.plot(times,inflection)
plt.xlabel('Time (seconds)')
plt.ylabel('Width')
plt.show()

