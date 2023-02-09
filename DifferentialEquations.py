# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:44:08 2023

@author: 14073
"""
import numpy as np
from scipy.integrate import odeint #this is the same as ODE15s
import matplotlib.pyplot as plt

#function that returns dy/dt
def model(x,t):
        v = x[0]
        y= x[1]
        k = 1
        m = 1
        return ( -(k/m)*y - .1*v, v) #the subtraction is friction (we see dampening)

#initial condition
x0=[0,1]

#time points
t=np.linspace(0,40,400)

#solve ODE
x = odeint(model,x0,t)

#plot results
plt.plot(t,x[:,1])
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()

# maple may be able to do it, matlab has the clumsy symbolic math toolbox