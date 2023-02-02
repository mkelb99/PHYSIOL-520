#-*-coding: utf-8 -*-) 
#Simulation of random particle steps to left or right: Example4_1.py
#Description: diffusion in one-dimension starting from a square-wave distribution
 
import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 

N = 100 # number of x locations (make even)
L = 0.1 # cm
dx = L/N # cm (space step)
x = np.linspace(dx/2,L-dx/2,N) # cm
D = 1e-5 # cm^2 / sec
dt = 0.1 # sec (time step)

# Diffusion operator:
A = np.zeros([N,N])

A[0,[0,1]] = (D/(dx**2))*np.array([-1, 1])
for i in range(1,N-1):
  A[i,[i-1,i,i+1]] = (D/(dx**2))*np.array([1,-2, 1]);

A[N-1,[N-2,N-1]] = (D/(dx**2))*np.array([1, -1]);

# Compute Crank-Nicholson operator
K = np.matmul( inv(np.identity(N) - (dt/2)*A), (np.identity(N) + (dt/2)*A) )

# Initial condition:
Cleft = np.ones(int(N/2))
Cright = np.zeros(int(N/2))
C = np.hstack([Cleft,Cright])
C0 = C

# Diffusion for 100 time steps
for i in range(100): 
  C = np.matmul(K,C)
  
plt.plot(x,C0,x,C)  
