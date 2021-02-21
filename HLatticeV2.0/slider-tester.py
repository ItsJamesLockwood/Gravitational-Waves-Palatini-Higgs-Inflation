# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 17:50:58 2021

@author: James
"""
import numpy as np
import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyl
from matplotlib.contour import QuadContourSet
from matplotlib.widgets import Slider

#Define display parameters
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
delta = 0.025

#Define model parameters
alpha = .5
beta = .5
x_bar, a, b, c = 2, 0, 1, .1
v = np.arange(0, 10, delta)
w = np.arange(0, 10, delta)

def compute_and_plot(ax, alpha):
    #Calculate grid values
    V, W = np.meshgrid(v,w)
    Z = (V**(beta))*(W**(1-beta))
    X = x_bar + a + b*Z
    U = alpha*np.log(V) + (1-alpha)*np.log(X) - c*(W+V)

    CS = QuadContourSet(ax, V, W, U, 200)
    pyl.clabel(CS, inline=1, fontsize=10)

# Plot
fig = pyl.figure()
pyl.title('Simplest default with labels')
ax = fig.add_subplot(221)
compute_and_plot(ax, alpha)

#Define slider for alpha
axcolor = 'lightgoldenrodyellow'
alpha_axis  = pyl.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
alpha_slider = Slider(alpha_axis, 'Amp', 0, 1, valinit=.5)

def update(ax, val):
    alpha = alpha_slider.val
    ax.cla()
    compute_and_plot(ax, alpha)
    pyl.draw()

alpha_slider.on_changed(lambda val: update(ax, val))

pyl.show()
plt.show()