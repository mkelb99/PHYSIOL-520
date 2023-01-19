# -*- coding: utf-8 -*-
# Example3_5.py
import numpy as np
import matplotlib.pyplot as plt

def my_Gaussian(x):
    # equation (3.1) from Baylor book
    return M*np.exp(-x*x/(4*D*t))/(np.sqrt(4*np.pi*D*t))

x1 = np.array([-0.090,-0.039,0,0.036,0.092,0.136,0.183])
y1 = np.array([0.2810,0.5965,0.7184,0.6138,0.2939,0.0778,0.0259])
plt.plot(x1,y1,'s')

# initial guess for parameters
t = 6000
D = 1.0e-6
M = 0.4

X = np.linspace(-0.50,0.50,500)
Y = my_Gaussian(X)
#plt.plot([],[],'',label='D values')
plt.plot(X,Y,label=D)
plt.xlabel('x (cm)')
plt.ylabel('concentration')
plt.show()

do = input('Would you like to take another guess? (y/n)')

while do == 'y':
    x = input('D={:.8}; enter value for D: '.format(D))
    D = float(x)

    x = input('M={:.8}; enter value for M: '.format(M))
    M = float(x)

    plt.plot(x1,y1,'s')
    Y = my_Gaussian(X)
    plt.plot(X,Y,label=D)

    plt.show()
    do = input('Would you like to take another guess? (y/n)')