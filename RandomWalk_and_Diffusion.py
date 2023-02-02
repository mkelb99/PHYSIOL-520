#-*-coding: utf-8 -*-) 
#Simulation of random particle steps to left or right: Example4_1.py
#Description: diffusion in one-dimension starting from a square-wave distribution
# from Baylor book, chapter 4
 
import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 
from numpy.random import random as rng #import a psuedo-random-number generato

from numpy.random import seed
import datetime 

N = 50 # number of x locations (make even)
L = 0.1 # cm
dx = L/N # cm (space step)
x = np.linspace(dx/2,L-dx/2,N) # cm
Max_points = 10000
Min_points = 0
N_frames = 40

#x = np.linspace(dx/2,L-dx/2,N) # cm
YA = Max_points*np.ones(int(N/2),dtype=int)
YB = Min_points*np.ones(int(N/2),dtype=int)
Y = np.hstack([YA,YB])
#Y[5] = 10000
Y0 = 1*Y
Spread = np.zeros(N_frames,dtype=int)

seed(47) #seed(l) will start the random number generator in a reproducible way
print("Total number of particles at start = {} ".format(np.sum(Y)))
NN=4 #number of times to print out how the calculation is progressing
now0 = datetime.datetime.now() #get the date & time at the start of the calculation

#Calculate the particle profiles for a total number of iterations/frames = N frames
for i in range(0,N_frames):
#  print(i) # added by DAB to keep track of simulation
  Z=0*Y             #Z, the 'next' particle profile, starts with all zeros
  if (i!=0 and i%(N_frames/NN)==0):
    #print out the progress of calculation every 100/NN percent of the way
    now1 = datetime.datetime.now() #get the current time and print elapsed time
    print("{} out of {} frames; elapsed minutes = {: .2f}".format(i,N_frames,60*\
    (now1.hour-now0.hour)+(now1.minute-now0.minute )+(now1.second\
    -now0.second)/60))
  for j in range(N): #run through all the x-locations in the profile
    if j == 0: #treat the first x-location as a special case
      for k in range(Y[j]): #decide whether to move each particle at first loc
        if rng() < 0.6667: Z[j] += 1 #no, keep the particle at first loc
        else: Z[j+1] += 1 #yes, move the particle right one loc
    elif j == N-1: #also treat the last x-location as a special case
      for k in range(Y[j]): #decide whether to move each particle at last loc
        if rng() < 0.6667: Z[j] += 1 #no, keep the particle at last loc
        else: Z[j-1] += 1 #yes, move the particle left one loc
    else:
      for k in range(Y[j]): #decide how to move each particle at locationj
        rn = rng()
        if rn < 0.3333: Z[j-1] += 1 #move the particle left one loc
        elif rn > 0.6667: Z[j+1] += 1 #move the particle right one loc
        else: Z[j] += 1 #keep particle where it is
  Y = 1 * Z #Y now becomes the profile for the 'current' frame
  #get some statistics related to the particle spread of the current profile
  index_l=0;index_2=0
  for j in range(N):
    if index_l==0 and Z[j]<(0.75*Max_points): index_l=j
    if index_2==0 and Z[j]<(0.25*Max_points): index_2=j

  Spread[i] = index_2-index_l #record points for 3/4 to 1/4 fall of the Z profile 
  
# Now the diffusion equation code
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
for i in range(50): 
  C = np.matmul(K,C)
    

plt.figure() 
plt.plot(x,Y/Max_points,'.-k',markersize=4) #this is the last profile; plot as small black dots 
plt.plot(x,C0,x,C) 
plt.xlabel('x (cm)') 
plt.ylabel('Concentration') 
plt.title("Simulation of 1-D Diffusion Starting from a Square Wave") 
#plt.savefig("Example4_1a.png", dpi=200) 




