# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:38:54 2021

@author: James
"""

from physUtils import *
from plotUtils import *

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab

import pandas as pd
import numpy as np


#%% Files
slice_test_p = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\tanh-slice-test%i_screen.log" 
slice_i = 2
slice_f = slice_test_p % slice_i

mayavi_file = slice_f
pw_field_number = 1
form = 'log'


fields_path = trim_name(mayavi_file) + '_whole_field_%i.%s'%(pw_field_number,form)
momenta_path = trim_name(mayavi_file) + '_momenta_%i.%s'%(pw_field_number,form)
slice_f_path=  trim_name(mayavi_file) + '_slice_f_%i.%s'%(pw_field_number,form)
slice_p_path=  trim_name(mayavi_file) + '_slice_p_%i.%s'%(pw_field_number,form)


#%% Get sim settings variables
try:
    strs, vals = sim_settings(mayavi_file)
    RESOLUTION = vals[6] 
    SKIP = int(vals[11]/vals[10] )
    BOXSIZE_H = vals[1] 
    Mpl = vals[14]
    CONF = vals[7]
    sim_vals = {'skip':SKIP, 'lh': BOXSIZE_H, 'mpl': Mpl, 'resolution':RESOLUTION,'conformal':CONF}
    for k,v in sim_vals.items():
        print(k,': ',v)
except FileNotFoundError:
    print("Did not update settings...")
    BOXSIZE_H = 15
    SKIP = 5
    RESOLUTION = 64

#%% import relevant files
data = import_screen(mayavi_file)
pw_data1, pw_data2 = import_pw(trim_name(mayavi_file) + '_pw_%i.%s'%(pw_field_number,form))

n_df = n_k(pw_data1,pw_data2,L=RESOLUTION)

my_rows = range(1,n_df.shape[0],2)


#%% Simple tests

sfdf = import_slice(slice_f_path,strict=False)
spdf = import_slice(slice_p_path,strict=False)

fdf = import_mesh(fields_path,)
pdf = import_mesh(momenta_path)


n= RESOLUTION
vals = np.linspace(1,n,n)
X,Y,Z = np.meshgrid(vals, vals, vals)
T = fdf.mesh[0]
#%%

def plot_a_slice(df,cut=-1,slice_ind=32, n=64):
    n= RESOLUTION
    vals = np.linspace(1,n,n)
    X,Y,Z = np.meshgrid(vals, vals, vals)
    T = fdf.mesh[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    #ax.set_zlabel('Z')
    
    #ax.countourf(X[:,:,slice_ind],Y[:,:,slice_ind])
    ax.countourf(T[:,:,slice_ind])
    ax.set_title("Heatmaps of the field perturbations for constant z=%i"%slice_ind)
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#scat = ax.scatter(X, Y, Z, c=T.flatten(), alpha=0.5)

#fig.colorbar(scat, shrink=0.5, aspect=5)
#plt.show()
#%% Mayavi tests
'''

s = T
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='x_axes',
                            slice_index=20,
                        )
mlab.outline()
mlab.show()

'''