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
    
    
def mayavi_cube(s):
    scf = mlab.pipeline.scalar_field(s)
    vol = mlab.pipeline.volume(scf,vmin=0.008)
    mlab.colorbar(vol)
    
if __name__=="__main__":
    mayavi_cube(s)
    #mayavi_plane(s)
    #mlab.plot3d(X.flatten(),Y.flatten(),Z.flatten(),s.flatten())
    #mlab.volume_slice(s)
    #s1 = mlab.pipeline.scalar_field(s)
    #mlab.pipeline.volume(s1)
    #mlab.pipeline.scalar_field(s)
    
    
    
    
    
    mlab.show()
    