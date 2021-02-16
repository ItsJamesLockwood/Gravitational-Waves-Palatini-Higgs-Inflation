# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 20:15:06 2021

@author: James
"""

''' A program for solving second order ODE and simulate the motion of pendulum
By Utkarsh Palwekar'''

global vi, dvi
def v2(phi,m):
    return phi**2 * m**2 /2
def dv2(phi,m):
    return phi * m**2 /2

def v4(phi,l):
    return phi**4 * l/4
def dv4(phi,l):
    return phi**3 * l

vi = v4
dvi = dv4
import numpy as np
from math import sqrt

import matplotlib.pyplot as plt

# Parameters
global m,l,M

'''
l = 1
m = 10**-0
M = 10**6
f0 = 10**0.5 * M
f1 = 0#10**(-3.5) * m
a0 = 1
'''
l = 1
m = 10**-0
M = 10**0
f0 = 10**0.5 * M
f1 = 0#10**(-3.5) * m
a0 = 1


init_vals = np.array([f0,a0,f1])
# Equations:
   
def V(u,t):
    f, a, f_t = u
    
    y = f_t
    a_t = np.sqrt(8*np.pi/(3*M**2) * (y**2 /2 + vi(f,m))) * a
    y_t = -dvi(f,m) - 3*a_t/a * y
    
    return np.array([ y, a_t, y_t]) 
    #(c*dv+k*dx)/mr, -(F+c*dv+k*dx)/m ])

def rk4(f, u0, t0, tf , n):
    t = np.linspace(t0, tf, n+1)
    u = np.array((n+1)*[u0])
    h = t[1]-t[0]
    for i in range(n):
        k1 = h * f(u[i], t[i])    
        k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
        k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
        k4 = h * f(u[i] + k3, t[i] + h)
        u[i+1] = u[i] + (k1 + 2*(k2 + k3) + k4) / 6
    return u, t


tmax = 1000
steps = 100*tmax
u, t  = rk4(V, init_vals , 0. , tmax , steps)
x1, x2, v1 = u.T
plt.plot(t, x1)
#plt.plot(t,1/t)
#plt.plot(t,x2)
plt.grid('on')
plt.show()