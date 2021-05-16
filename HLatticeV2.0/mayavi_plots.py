# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 03:22:22 2021

@author: James
"""

from mayavi import mlab
import numpy as np

from plotUtils import *
from physUtils import *

import os

x, y, z = np.ogrid[-2:2:100j, -2:2:100j, -2:2:100j]

np_file = "m_26.npy"
#mesh = np.load(np_file)

j = 12
wd = os.getcwd()
my_mesh = "sim16_m%i.npy"%j
my_mesh = "m%ih.npy" % 45
if not os.path.exists(os.getcwd()+'\\' +my_mesh):
    print("Saving mesh to file")
    _, mesh = import_a_mesh('data\\sim16_run2_whole_field_1.log', selected=j)
    np.save(my_mesh,mesh)
else:
    print("Loading preexisting file.")
    mesh = np.load(my_mesh)

coords = np.linspace(1,mesh.shape[0],mesh.shape[0])
X,Y,Z = np.meshgrid(coords,coords,coords)
#s = np.sin(x*y*z + x + y*z)/(x*y*z + x + y*z)
s = mesh
s = np.log10(np.abs(mesh-mesh.mean()))
print(s.shape)


def mayavi_plane(s):
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                plane_orientation='x_axes',
                                slice_index=20,
                            )
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                plane_orientation='y_axes',
                                slice_index=20,
                            )
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                plane_orientation='z_axes',
                                slice_index=20,
                            )
    mlab.outline()
    
    
def mayavi_cube(s,vmin=0.008):
    scf = mlab.pipeline.scalar_field(s)
    vmax = (-vmin) + 0.5 * (s.max()-(-vmin))
    vmax = s.max()
    vol = mlab.pipeline.volume(scf,vmin=-vmin, vmax= vmax)
    mlab.colorbar(vol)

def mayavi_iso(s,contours=[-5.6,-5.4]):
    magnitude = mlab.pipeline.extract_vector_norm(s)
    scal = mlab.pipeline.scalar_field(s)
    #mlab.pipeline.iso_surface(magnitude, contours=[-5.0])
    mlab.pipeline.iso_surface(scal, contours=contours,opacity=0.7)
    mlab.show()
    mlab.close(all=True)
    #mlab.outline()
    #vol = mlab.pipeline.volume(scf)
    #mlab.colorbar(vol)
   
def mayavi_altiso(s,contours=[-5.6,-5.4]):
    from enthought.mayavi.modules.iso_surface import IsoSurface
    
def mayavi_plot3d(s,res=128,cutoff=5.5):
    """Currently deprecated"""
    coords = np.linspace(1,res,res)

    X,Y,Z = np.meshgrid(coords,coords,coords)
    x_ = X.flatten()
    y_ = Y.flatten()
    z_ = Z.flatten()
    s_ = s.flatten()
    x2_, y2_, z2_, s2_ = [],[],[],[]
    print(len(s_),len(x_))
    '''
    for i in range(len(s_)):
        if s_[i]>-cutoff:
            x2_.append(x_[i])
            y2_.append(y_[i])
            z2_.append(z_[i])
            s2_.append(s_[i])
    '''
    x2_, y2_, z2_, s2_ = np.array(x2_),np.array(y2_),np.array(z2_),np.array(s2_)
    print("Length reduction: %i -> %i"%(len(x_), len(x2_)))
    #mlab.plot3d(x2_,y2_,z2_,s2_)
    #return x2_, y2_, z2_, s2_ 
    #scf = mlab.pipeline.scalar_scatter(x2_,y2_,z2_,s2_)
    scf = mlab.pipeline.scalar_field(s)
    
    vol = mlab.pipeline.volume(scf,vmin=0.008)
    mlab.colorbar(vol)

if __name__=="__main__":
    #mayavi_cube(s, vmin=5.80)
    #mayavi_plane(s)
    #mayavi_plot3d(4)
    #mayavi_iso(s)
    mayavi_altiso(s)
    #mlab.plot3d(X.flatten(),Y.flatten(),Z.flatten(),s.flatten())
    #mlab.volume_slice(s)
    #s1 = mlab.pipeline.scalar_field(s)
    #mlab.pipeline.volume(s1)
    #mlab.pipeline.scalar_field(s)
    
    
    
    
    
    mlab.show()
    