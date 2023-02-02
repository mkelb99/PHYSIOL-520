#-*-coding: utf-8 -*-) 
#Simulation of random particle steps to left or right: Example4_1.py
#Description: diffusion in one-dimension starting from a square-wave distribution
# from Baylor book, chapter 4
 
import numpy as np 
import matplotlib.pyplot as plt 
from numpy.random import random as rng #import a psuedo-random-number generato

from numpy.random import seed
import datetime 

X_locs = 2*15
Max_points = 10000
Min_points = 0
N_frames = 10

X = np.linspace(0,X_locs-1,X_locs,dtype=int)
YA = Max_points*np.ones(int(X_locs/2),dtype=int)
YB = Min_points*np.ones(int(X_locs/2),dtype=int)
Y = np.hstack([YA,YB])
#Y[15]=10000
Y0 = 1*Y
Spread = np.zeros(N_frames,dtype=int)

seed(4) #seed(l) will start the random number generator in a reproducible way
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
  for j in range(X_locs): #run through all the x-locations in the profile
    if j == 0: #treat the first x-location as a special case
      for k in range(Y[j]): #decide whether to move each particle at first loc
        if rng() < 0.6667: Z[j] += 1 #no, keep the particle at first loc
        else: 
            Z[j+1] += 1 #yes, move the particle right one locZ[j] += 1
            Z[j] += 1 #keep the first position value constant
    elif j == 1: #make sure the second location doesn't go back into the first
        for k in range(Y[j]):
            if rng() < 0.6667: Z[j] += 1 #keep the particle where it is
            else: Z[j+1] += 1 #move the particle to the right
    elif j == X_locs-2: #prevent the last bin from increasing in concentration
      for k in range(Y[j]): #decide whether to move each particle at the second to last loc
        if rng() < 0.6667: Z[j] += 1 #no, keep the particle at last loc
        else: Z[j-1] += 1 #yes, move the particle left one loc
    else:
      for k in range(Y[j]): #decide how to move each particle at locationj
        rn=rng()
        if rn < 0.3333: Z[j-1] += 1 #move the particle left one loc
        elif rn > .6667: Z[j+1] += 1 #move the particle right one loc
        else: Z[j] += 1
  Y = 1 * Z #Y now becomes the profile for the 'current' frame
  #get some statistics related to the particle spread of the current profile
  index_l=0;index_2=0
  for j in range(X_locs):
    if index_l==0 and Z[j]<(0.75*Max_points): index_l=j
    if index_2==0 and Z[j]<(0.25*Max_points): index_2=j

  Spread[i] = index_2-index_l #record points for 3/4 to 1/4 fall of the Z profile 
  


now2 = datetime.datetime.now() #get the time at the end of the calculation 
print(' {} out of {} frames; minutes of calculation = {: .2f} '.format(N_frames,\
 N_frames,60*(now2.hour-now0.hour )+(now2.minute-now0.minute )+\
 (now2.second-now0.second)/60)) 

#get some statistics related to the profile calculated in the last frame 
print("Total number of particles at end = {} ".format(np.sum(Y))) 
print("Number of particles in first half= {} ".format(np.sum(Y[0:int(X_locs/2)]))) 
print("Number of particles in second half= {} ".format(np.sum(Y[int(X_locs/2):X_locs])))
print("Spread of final profile (3/4 to 1/4) = {} ".format(index_2-index_l)) 



#plot the first and last particle profiles 
plt.figure() 
plt.plot(X,Y0,'blue') #this is the first profile; plot as connected curve in gray 
plt.plot(X,Y,'.k',markersize=4) #this is the last profile; plot as small black dots 
#plt.xlim(0,X_locs) #range ofx-locations to view 
#plt.ylim(0,1.2*Max_points) #vertical scale for number of particles 
plt.xlabel('x-location') 
plt.ylabel('Number of Particles') 
plt.title("Simulation of 1-D Diffusion Starting from a Square Wave") 
#plt.savefig("Example4_1a.png", dpi=200) 


"""
#also plot the frame-by-frame spread of the profiles 
plt.figure() 
plt.plot(range(N_frames),Spread,'.k',markersize= 1) 
plt.xlim(0,N_frames) 
plt.ylim(0,1.6*np.sqrt(N_frames)) 
plt.xlabel('Frame Number') 
plt.ylabel('Spread of Profile (3/4-way to 1/4-way down)') 
plt.title("Simulation of 1-D Diffusion Starting from a Square Wave") 
#plt.savefig("Example4_1b.png", dpi=200) 
"""





