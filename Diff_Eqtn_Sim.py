# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 09:50:46 2023

@author: 14073
"""
import numpy as np
N=50 #number of x locations (make even)
L=.1 #cm
dx=L/N #cm (space step)
x=np.arange(dx/2,L,dx) #cm
D = 1e-5 #cm^2/sec
dt = .02 #sec (time step)

#Diffusion Operator 
A= np.zeros([50,50])
A[1,1:2] = np.multiply((D/(dx**2)),np.array[-1,1])
