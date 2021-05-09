# -*- coding: utf-8 -*-
"""
Created on Wed May  5 11:50:40 2021

@author: James
"""
from physUtils import *
from plotUtils import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import matplotlib.animation as animation
import math as maths
import numpy as np
from PIL import Image

#%% Files
j = 9
inc = 3
mesh_tot=0
inds = [i for i in range(13,25)]
meshes = {};

for i in range(j,j+inc):
    path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\perturb_run%i_screen.log"
    path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\pert_mode_%i_screen.log"
    path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\pert_field_%i_screen.log"
    #path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\omp_test%i_screen.log"
    #path = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\romp_test%i_screen.log"
    
    #path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\omp_run%i_screen.log"
    
    #path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\pert_mode_all%i_screen.log"
    
    path_v = i
    path_f = path%path_v
    
    
    pertfile = path_f
    pw_field_number = 1
    form = 'log'
    #%%
    #%% Get sim settings variables
    try:
        strs, vals = sim_settings(pertfile)
        RESOLUTION = vals[6] 
        L = RESOLUTION
        SKIP = int(vals[11]/vals[10] )
        BOXSIZE_H = vals[1] 
        LH = BOXSIZE_H
        Mpl = vals[14]
    except FileNotFoundError:
        print("Did not update settings...")
        BOXSIZE_H = 15
        SKIP = 50 
        Mpl = 1
    #%% Import
    '''
    data = import_screen(pertfile)
    pw_data1, pw_data2 = import_pw(trim_name(pertfile) + '_pw_%i.%s'%(pw_field_number,form))
    
    n_df = n_k(pw_data1,pw_data2,L=RESOLUTION)
    
    my_rows = range(1,n_df.shape[0],2)
    '''
    
    fields_path = trim_name(pertfile) + '_whole_field_%i.%s'%(pw_field_number,form)
    momenta_path = trim_name(pertfile) + '_momenta_%i.%s'%(pw_field_number,form)
    slice_f_path=  trim_name(pertfile) + '_slice_f_%i.%s'%(pw_field_number,form)
    slice_p_path=  trim_name(pertfile) + '_slice_p_%i.%s'%(pw_field_number,form)
    
    #%% Analysis
    slice_f = import_slice(slice_f_path)
    mesh_tot += slice_f.mesh[0]
    meshes[i] = slice_f.mesh[0]
    if True:    
        fig,ax = plt.subplots()
        plt.axis('off')
        use_log = False
        if use_log:
            impl = ax.imshow(np.abs(slice_f.mesh[0]- slice_f.mesh[0].mean()),
                         cmap='hsv_r',
                         vmin = 1.e-12,
                         vmax = 1.e-6,
                         norm=matplotlib.colors.LogNorm())
        else:
            impl = ax.imshow(slice_f.mesh[0])
        #ax.set_title("Field values for slice in lattice")
        cutoff = [16,32,38,0,0,0,0,0,0,0]
        ax.set_title(r"Perturbations in XZ-slice for cutoff $|j_1|,|j_2|,|j_3|< %i$"%cutoff[i-j])
        #fig.savefig(trim_file_name(pertfile)+'.png')
        cbar = fig.colorbar(impl, ax=ax)

    #plt.close(fig)
    #plot_slices(slice_f)
mesh_tot /=1
slice_f.mesh[0] = mesh_tot
#plot_slices(slice_f)
if False:
    fig,ax = plt.subplots()
    plt.axis('off')
    plt.title("Sum of modes 13-24")
    ax.imshow(slice_f.mesh[0])
    #fig.savefig('pert_modes_13_24.png')
#plt.close(fig)

#%% Binary

def import_chk(path):
    ''' Assumed that the the chk-file is a non-binary (i.e. formatted) file.'''
    file = open(path,'r')
    
    res = int(file.readline())
    ns = int(file.readline())
    print(res,ns)
    meshes = {}
    momenta = {}
    temp_mesh, temp_momenta = {},{}

    for i in range(1,ns+1):
        meshes[i] = []
        momenta[i] = []
        temp_mesh[i] = []
        temp_momenta[i] = []
        
    l1 = file.readline().strip().split()   
    l2 = file.readline().strip().split()
    while l1!=[] and l2!=[]:
        for i in range(0,ns*res**2,ns):        
            for f in range(0,ns):
                meshes[f+1].append(float(l1[i+f]))
                momenta[f+1].append(float(l2[i+f]))
        l1 = file.readline().strip().split()   
        l2 = file.readline().strip().split()
    for i in range(1,ns+1):
        meshes[i] = np.array(meshes[i]).reshape(res,res,res).T
        momenta[i] = np.array(momenta[i]).reshape(res,res,res).T
    return meshes,momenta

              
from scipy.io import FortranFile
import struct

#chk = trim_name(pertfile) + '_fieldsB.log'

'''
with open(chk,'r') as file:
    res = file.readline()
    ns = file.readline()
    l1 = file.readline().strip()
'''

#f1, f2 = import_chk(chk)
                