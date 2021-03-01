# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 20:15:06 2021

@author: James
"""

''' A program for solving second order ODE and simulate the motion of pendulum
By Utkarsh Palwekar'''

global vi, dvi
global m,l,M,PM,xi


def v2(phi,m):
    return phi**2 * m**2 /2
def dv2(phi,m):
    return phi * m**2 /2

def v4(phi,l):
    return phi**4 * l/4
def dv4(phi,l):
    return phi**3 * l

def vpal(chi,l):
    xi = 3.8*10**6 * 50**2 * l
    a = l * M**4 /4 / xi**2
    b= np.sqrt(xi)/M
    
    return a * np.tanh(b*chi)**4

def dvpal(chi,l):
    xi = 3.8*10**6 * 50**2 * l
    a = l *M**4 /4 / xi**2
    b= np.sqrt(xi)/M
    
    return 4*a*b* np.tanh(b*chi)**3 / np.cosh(b*chi)**2


vi = v4
dvi = dv4
import numpy as np
from math import sqrt

import matplotlib.pyplot as plt

# Parameters


l = 1
m = 10**-0
M = 10**6
f0 = 10**0.5 * M
f1 = 0#10**(-3.5) * m
a0 = 1
xi = 3.8*10**6 * 50**2 *l

'''
l = 10**-4
xi = 3.8*10**6 * 50**2 *l
m = 10**-0
M = 1024
PM = np.sqrt(8*np.pi)*M
f0 = 10**0.5 * M
f0= 1.5*PM/np.sqrt(xi)
f1 = 0#10**(-3.5) * m
f1 = -10**-20 * M**2
a0 = 1
'''

init_vals = np.array([f0,a0,f1])
# Equations:
   
def V(u,t):
    f, a, f_t = u
    y = f_t
    a_t = np.sqrt(8*np.pi/(3*M**2) * (y**2 /2 + vi(f,m))) * a
    y_t = -dvi(f,m) - 3*a_t/a * y
    if np.cosh(xi**.5/M * f)==np.inf:
        print("Infinite at",f)
    return np.array([ y, a_t, y_t]) 
    #(c*dv+k*dx)/mr, -(F+c*dv+k*dx)/m ])

def rk4(f, u0, t0, tf , n):
    t = np.linspace(t0, tf, n+1)
    u = np.array((n+1)*[u0])
    h = t[1]-t[0]
    last_print = 0
    for i in range(n):
        k1 = h * f(u[i], t[i])    
        k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
        k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
        k4 = h * f(u[i] + k3, t[i] + h)
        u[i+1] = u[i] + (k1 + 2*(k2 + k3) + k4) / 6
        if (int(i/n*100)%10==0) and int(i/n*100)!=last_print: 
            print(int(i/n*100),'%')
            last_print = int(i/n*100)
    return u, t


tmax = 1000000
t0 = 800000
steps = 1000*tmax
steps = 100000
u, t  = rk4(V, init_vals , t0 , tmax , steps)
x1, x2, v1 = u.T
power = 2/3
scale = t[-1] / t[-1]**power
#plt.plot(t**power * scale, x1)
plt.plot(t**1, x1)
#plt.plot(t,1/t)
#plt.plot(t,x2)
plt.grid('on')
plt.show()