# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:38:54 2021

@author: James
"""

from physUtils import *
from plotUtils import *

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
#from mayavi import mlab

import pandas as pd
import numpy as np

import plotly as pl
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import copy
pio.renderers.default = 'browser'

#%% Files
slice_test_p = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\tanh-slice-test%i_screen.log" 
cernboxR6a = r"D:\Physics\MPhys Project\DatasetArcive\CERNBox\Remote SIm 6a\data\Remote_SIM6_f5ec45f_screen.log"

slice_i = 2
slice_f = slice_test_p % slice_i

mayavi_file = slice_f
mayavi_file = cernboxR6a

pw_field_number = 1
form = 'log'

want_f = 'yes'
want_s = 'no'


if want_f=='yes':
    fields_path = trim_name(mayavi_file) + '_whole_field_%i.%s'%(pw_field_number,form)
    momenta_path = trim_name(mayavi_file) + '_momenta_%i.%s'%(pw_field_number,form)
if want_s=='yes':
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
'''
sfdf = import_slice(slice_f_path,strict=False)
spdf = import_slice(slice_p_path,strict=False)

fdf = import_mesh(fields_path,)
pdf = import_mesh(momenta_path)


n= RESOLUTION
vals = np.linspace(1,n,n)
X,Y,Z = np.meshgrid(vals, vals, vals)
T = fdf.mesh[0]
'''
#%%

def plotly_a_mesh(mesh_pair,n=64,cutoff=-3):
    a_val = mesh_pair[0]
    the_mesh = mesh_pair[1]
    coords = np.linspace(1,n,n)
    X,Y,Z = np.meshgrid(coords, coords, coords)
    xs = X.reshape(1,n**3)[0]
    ys = Y.reshape(1,n**3)[0]
    zs = Z.reshape(1,n**3)[0]
    T = the_mesh.reshape(1,n**3)[0]

    
    xs,ys,zs,T = mesh_cutoff(mesh_pair,n=n,cutoff=cutoff)
    print(xs.shape)
    print(ys.shape)
    
    
    '''
    axis_style=dict(showline=True)
    layout=dict(width=1200, height=900, title='My plot',
            xaxis=axis_style,
            yaxis=axis_style,
           hovermode='closest')
    '''
    
    fig = go.Figure(data=[go.Scatter3d(x=xs,y=ys,z=zs,
                   mode = 'markers',
                   marker=dict(
                       color=T,
                       colorbar=dict()))])
    
    fig.update_layout(margin=dict(l=0,r=0,b=0,t=0))
    fig.show()
    return fig

    '''
    cmap = plt.get_cmap('jet')
    
    cNorm = matplotlib.colors.Normalize(vmin=the_mesh.min(),vmax=the_mesh.max())
    scalarMap = cm.ScalarMappable(),norm=cNorm,cmap=camp)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter()
    '''
def matplot_a_mesh(mesh_pair,n=64,cutoff=-3):
    a_val = mesh_pair[0]
    the_mesh = mesh_pair[1]
    coords = np.linspace(1,n,n)
    X,Y,Z = np.meshgrid(coords, coords, coords)
    xs = X.reshape(1,n**3)[0]
    ys = Y.reshape(1,n**3)[0]
    zs = Z.reshape(1,n**3)[0]
    T = the_mesh.reshape(1,n**3)[0]

    
    xs,ys,zs,T = mesh_cutoff(mesh_pair,n=n,cutoff=cutoff)
    print(xs.shape)
    print(ys.shape)

    fig = plt.figure()    
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    scat = ax.scatter(xs, ys, zs, c=T, alpha=0.5)
    
    fig.colorbar(scat, shrink=0.5, aspect=5)
    plt.show()

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

'''
the_mesh = mesh_pair[1]
coords = np.linspace(1,n,n)
X,Y,Z = np.meshgrid(coords, coords, coords)
T = the_mesh

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

scat = ax.scatter(X, Y, Z, c=T.flatten(), alpha=0.5)

fig.colorbar(scat, shrink=0.5, aspect=5)
plt.show()
'''

def mesh_cutoff(mesh_pair, cutoff=-3,n=64):
    a_val = mesh_pair[0]
    the_mesh = mesh_pair[1]
    coords = np.linspace(1,n,n)
    X,Y,Z = np.meshgrid(coords, coords, coords)
    xs = X.reshape(1,n**3)[0]
    ys = Y.reshape(1,n**3)[0]
    zs = Z.reshape(1,n**3)[0]
    T = the_mesh.reshape(1,n**3)[0]
    
    mesh_avg = the_mesh.mean()
    mesh_perts = np.log10(np.abs(the_mesh - mesh_avg)).reshape(1,n**3)[0]
    
    
    x_keep, y_keep, z_keep, t_keep = [], [], [], [] 
    for i in range(len(xs)):
        if mesh_perts[i]>cutoff:
            x_keep.append(xs[i])
            y_keep.append(ys[i])
            z_keep.append(zs[i])
            t_keep.append(T[i])
    return np.array(x_keep),np.array(y_keep),np.array(z_keep),np.array(t_keep)



#%%




def template_pl():
    # Helix equation
    t = np.linspace(0, 20, 100)
    x, y, z = np.cos(t), np.sin(t), t
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=12,
            color=z,                # set color to an array/list of desired values
            colorscale='Viridis',   # choose a colorscale
            opacity=0.8
        )
    )])
    
    # tight layout
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.show()


mesh_pair = import_a_mesh(fields_path,selected=25)
#%%
m_temp = copy.deepcopy(mesh_pair)
k=128
m_temp[1] = m_temp[1][:k,:k,:k]
print(m_temp[1].shape)
print(mesh_pair[1].shape)

matplot_a_mesh(m_temp,n=k,cutoff=-2.2)


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