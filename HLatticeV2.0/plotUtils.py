# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:38:40 2020

@author: James
"""
import os
import scipy
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib
import matplotlib.animation as animation
import matplotlib.ticker as tikr
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
from PIL import Image

import imageio

from physUtils import *

#%% Setup loggers
CWD_PATH = os.getcwd()
LOG_PATH = CWD_PATH + r"\importer_logs"
if not os.path.exists(LOG_PATH):
    os.mkdir(LOG_PATH)
LOG_FILE = datetime.now().strftime("importer_%y%m%d-%H_%M_%S.txt")
PATH = LOG_PATH +r"\\" +  LOG_FILE

logger = logging.getLogger()
logger.setLevel(logging.INFO)

fmter = logging.Formatter("%(asctime)s:%(levelname)s:%(name)s: %(message)s","%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(PATH)
fh.setLevel(logging.INFO)
fh.setFormatter(fmter)
logger.addHandler(fh)

#%% Utility definitions

def  trim_name(file):
    ind = file.index('_screen')
    return file[:ind]


def trim_file_name(file):
    ind1 = file.index('_screen')
    ind2 = file.rindex('\\')+1
    return file[ind2:ind1]

def max_field_number(screenfile):
    file_exists = True
    fld = 0
    while file_exists:        
        fld +=1
        path = trim_name(screenfile) + '_pw_%i.log'%fld
        file_exists = os.path.exists(path)
    return (fld-1)
    
def import_screen(file_name,print_logs=False):
    #old_col_names = [chr(i) for i in range(ord('a'),ord('a')+10)]
    fld = max_field_number(file_name)
    col_names = ['a',
    			'h',
    			'omega',
    			'pratio',
    			'kratio',
    			'gratio']
    means = ['mean%i'%i for i in range(1,fld+1)]
    rmss = ['rms%i'%i for i in range(1,fld+1)]
    col_names += means + rmss
    
    logger.info("Opening screen file... ")
    data = pd.read_csv(file_name,delim_whitespace=True,index_col=False,skiprows=1,header=None)
    with open(file_name,'r') as file:
        first_line = file.readline().strip().split()
    
    
    #Check if metric perturbations are switched on
    if first_line[2]=='E_f':
        print("colnames:",len(col_names))
        print("data:",len(data.columns))
        col_names = ['a','h','ef','pratio','kratio','gratio','kg','gg','eratio','rmsh']
        col_names += means + rmss
        data.columns = col_names
    else:
        print("colnames:",len(col_names))
        print("data:",len(data.columns))
        data.columns = col_names
    
    logger.info("Done.")
    logger.info("Formatting 'data': removing non numerical lines...")
    try:
    	data = data[~data.iloc[:,0].str.contains("[a-zA-Z]").fillna(False)]
    except AttributeError:
    	pass
    logger.info("Done")
    
    logger.info("Applying pd.to_numeric to 'data'")
    data = data.apply(pd.to_numeric)
    logger.info("Finished importing 'data'.")
    
    if print_logs==True:
        print(data.columns)
        print(data.shape)
        print(data.dtypes)
        print(data.iloc[:5])
        print(data.dtypes)
        
    return data


def import_GW(GW_file):
    file = open(GW_file,'r')
    
    a = []
    psf = []
    psdf = []
    
    while True:
        l1 = file.readline()
        l2 = file.readline()
        l3 = file.readline()
        
        
        if l1=='' or l2=='' or l3=='':
            break
        d1 = float(l1.strip())
        d2 = list(map(float,l2.strip().split()))
        d3 = list(map(float,l3.strip().split()))

        a.append(d1)
        psf.append(d2)
        psdf.append(d3)
    df1 = pd.DataFrame(psf)
    df2 = pd.DataFrame(psdf)
    
    df1['a'] = pd.Series(a)
    df2['a'] = pd.Series(a)
    file.close()
    return df1,df2

def import_pw(pw_file):
    file = open(pw_file,'r')
    global a,psf,psdf
    a = []
    psf = []
    psdf = []
    
    while True:
        global l1,l2,l3,d1,d2,d3
        l1 = file.readline()
        l2 = file.readline()
        l3 = file.readline()
        
        
        if l1=='' or l2=='' or l3=='':
            break
        d1 = float(l1.strip())
        d2 = list(map(float,l2.strip().split()))
        d3 = list(map(float,l3.strip().split()))
        
        a.append(d1)
        psf.append(d2)
        psdf.append(d3)

    df1 = pd.DataFrame(psf)
    df2 = pd.DataFrame(psdf)
    df1['a'] = pd.Series(a,dtype=np.float64)
    df2['a'] = pd.Series(a,dtype=np.float64)
    file.close()
    return df1,df2


def import_fields(fields_file,sep="SEPARATOR"):
    file = open(fields_file,'r')
    field_list = []

    l1 = file.readline()
    if l1.strip()!=sep:
        raise ValueError("Expected separator in first line. Instead found '"+sep+"'.")
        
    while True:
        l2 = file.readline()
        if l1=='' or l2=='':
            break
        la = file.readline()
        lb = file.readline()
        while la!=sep and lb!=sep:
            row = []
            if la=='' or lb=='':
                break
            da = list(map(int,la.strip().split()))
            db = list(map(float,lb.strip().split()))
            row.append(float(l2.strip()))
            row.append(da[0])
            row.append(da[1])
            row.append(db)
            
            field_list.append(row)
            la = file.readline()
            lb = file.readline()
        
        if lb=='SEPARATOR':
            raise ValueError("The second line was a separator when a list of coordinates was expected. Check file formatting.")
    
    field_df = pd.DataFrame(field_list,columns=['a','x','y','zs'])
    file.close()
    return field_df

def import_mesh(fields_file,sep="SEPARATOR"):
    '''
    Import field which has been stored directly as a mesh - rather than the convoluted line by line output that
    is read in by import_fields.
    '''
    file = open(fields_file,'r')
    field_list = []

    lres = file.readline().strip()
    if lres==sep:
        raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
    res = int(lres)
        
    l1 = file.readline().strip()    
    l2 = file.readline().strip()
    l3 = file.readline().strip()
    while l1!='':        
        if l1!=sep:
            raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
        if l2==sep or l3==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        mesh = list(map(float, l3.split()))
        mesh = np.array(mesh).reshape(res,res,res).T
        
        field_list.append([a,mesh])
        
        l1 = file.readline().strip()
        l2 = file.readline().strip()
        l3 = file.readline().strip()
    
    mesh_df = pd.DataFrame(field_list, columns=['a','mesh'])
    file.close()
    return mesh_df
    

def plot_fields(field_df,cond=['x', 1],use_FFT=False,use_contour=False):
    plane1 = field_df[field_df[cond[0]]==cond[1]]
    fvals = np.array(plane1.zs.values.tolist())
    
    if use_FFT:
        fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    
        
    if cond[0]=='x':
        anticond = 'y'
    elif cond[0]=='y':
        anticond = 'x'
    else:
        raise ValueError("Currently the condition '"+cond[0]+"' is not supported.")
    X, Y = np.meshgrid(np.linspace(1,plane1[anticond].max(),plane1.shape[0]),np.linspace(1,plane1[anticond].max(),plane1.shape[0]))
    
    fig = plt.figure()
    plt.subplots_adjust(left=0.25, bottom=0.25)    
    ax = fig.add_subplot(2,1,1)
    ax3 = fig.add_subplot(2,1,2,projection='3d')
    ax.set_aspect(1)
    if use_contour:
        ax.contourf(fvals)
    else:
        ax.imshow(fvals)
    pl3d = ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
    cbar = fig.colorbar(pl3d, ax=ax3)
    
    axcolor = 'lightgoldenrodyellow'
    slid_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    x_slider = Slider(slid_ax, cond[0].upper(), 1, plane1[anticond].max(), valinit=cond[1], valstep=1)
    
    def update(val):
        xv = x_slider.val
        plane1 = field_df[field_df[cond[0]]==xv]
        #print("Asked for update...")
        fvals = np.array(plane1.zs.values.tolist())
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    
        #fvals = np.log(np.abs(fvals))
        ax.cla()
        ax3.cla()
        if use_contour:
            ax.contourf(fvals)
        else:
            ax.imshow(fvals)
        ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
        fig.canvas.draw_idle()
    x_slider.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])    
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')     
    def reset(event):
        print("Reset requested...")
        x_slider.reset()
        
    button.on_clicked(reset)
    resetax._button = button
    return fvals

def plot_mesh(field_df,a_ind=0,cond=['x',0],use_FFT=False,use_contour=False):
    the_mesh = field_df.iloc[a_ind,:].mesh
    
    if cond[0]=='x':
        plane1 = the_mesh[cond[1],:,:]
    elif cond[0]=='y':
        plane1 = the_mesh[:,cond[1],:]
    elif cond[0]=='z':
        plane1 = the_mesh[:,:,cond[1]]
    else:
        raise ValueError("Coordinate '%s' not recognised."%cond[0])
    print(plane1.shape)
    fvals = plane1
    if use_FFT:
        fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    X, Y = np.meshgrid(np.linspace(1,plane1.shape[0],plane1.shape[0]),np.linspace(1,plane1.shape[1],plane1.shape[1]))
    fig = plt.figure()
    plt.subplots_adjust(left=0.25, bottom=0.25)    
    ax = fig.add_subplot(2,1,1)
    ax3 = fig.add_subplot(2,1,2,projection='3d')
    ax.set_aspect(1)
    #ax.contourf(fvals)
    if use_contour:
        ax.contourf(fvals)
    else:
        ax.imshow(fvals)
    pl3d = ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
    cbar = fig.colorbar(pl3d, ax=ax3)
    
    axcolor = 'lightgoldenrodyellow'
    slid_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    x_slider = Slider(slid_ax, cond[0].upper(), 1, plane1.shape[0], valinit=cond[1], valstep=1)
    
  
    def update(val):
        xv = int(x_slider.val)
        if cond[0]=='x':
            plane1 = the_mesh[xv-1,:,:]
        elif cond[0]=='y':
            plane1 = the_mesh[:,xv-1,:]
        elif cond[0]=='z':
            plane1 = the_mesh[:,:,xv-1]
        else:
            raise("Code is broken: this else should not be reachable.")
        #print("Asked for update...")
        fvals = plane1
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))
        
        ax.cla()
        ax3.cla()
        if use_contour:
            ax.contourf(fvals)
        else:
            ax.imshow(fvals)
        ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
        fig.canvas.draw_idle()
        
    x_slider.on_changed(update)
    
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])    
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')     
    def reset(event):
        print("Reset requested...")
        x_slider.reset()
        
    button.on_clicked(reset)
    resetax._button = button
    return fvals

def import_slice(slice_file,sep="SEPARATOR",strict=True):
    '''
    Import slice of field which has been stored directly as a mesh. 
    '''
    file = open(slice_file,'r')
    field_list = []

    lres = file.readline().strip()
    if lres==sep:
        raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
    res = int(lres)

    l1 = file.readline().strip()    
    l2 = file.readline().strip()
    l3 = file.readline().strip()
    counter = 0
    while l1!='':      
        counter +=1
        if l1!=sep:
            if strict:
                raise ValueError("Expected separator in second line. Instead found '"+l1[:1000]+"' at line %i."%counter) 
            else:
                print("WARNING: skipping at line %i due to misplaced separator..."%counter)

            while l1!='' and l1!=sep:
                l1 = file.readline().strip()
            l2 = file.readline().strip()
            l3 = file.readline().strip()
        if l2==sep or l3==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        mesh = list(map(float, l3.split()))
        try:
            mesh = np.array(mesh).reshape(res,res).T
        except ValueError:
            print("WARNING: The line for a=",a," is the wrong size (", len(mesh),"). Expected: ",res,"^2. Please check file: ",slice_file)
            break
        field_list.append([a,mesh])
        
        l1 = file.readline().strip()
        l2 = file.readline().strip()
        l3 = file.readline().strip()
    
    mesh_df = pd.DataFrame(field_list, columns=['a','mesh'])
    file.close()
    return mesh_df

def plot_slices(slice_df,a_ind=0,use_FFT=False,use_contour=False, title=''):
    the_mesh = slice_df.iloc[a_ind,:].mesh
    fvals = the_mesh
    fft_str=''
    if use_FFT:
        fft_str='\n(log of modulus of FFT)'
        fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    X, Y = np.meshgrid(np.linspace(1,fvals.shape[0],fvals.shape[0]),np.linspace(1,fvals.shape[1],fvals.shape[1]))
    fig = plt.figure()
    plt.subplots_adjust(left=0.25, bottom=0.25)    
    plt.tight_layout()
    ax = fig.add_subplot(2,1,1)
    ax3 = fig.add_subplot(2,1,2,projection='3d')
    ax.set_aspect(1)
    ax.set_title("Heatmaps of the %s perturbations for constant x "%title + fft_str)
    ax.set_xlabel("Y coordinate")
    ax.set_ylabel("Z coordinate")

    #ax.contourf(fvals)
    if use_contour:
        ax.contourf(fvals)
    else:
        ax.imshow(fvals)
    pl3d = ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
    cbar = fig.colorbar(pl3d, ax=ax3, pad=0.2)
    
    axcolor = 'lightgoldenrodyellow'
    
    slid_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    a_slider = Slider(slid_ax, 'A', 0, slice_df.shape[0]-1, valinit=a_ind, valstep=1)
    
    
    def update(val):
        av = int(a_slider.val) 
        the_mesh = slice_df.iloc[av,:].mesh
        fvals = the_mesh
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))
        
        ax.cla()
        ax.set_title("Heatmaps of the %s perturbations for constant x "%title + fft_str)
        ax.set_xlabel("Y coordinate")
        ax.set_ylabel("Z coordinate")

        ax3.cla()
        if use_contour:
            ax.contourf(fvals)
        else:
            ax.imshow(fvals)
        ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
        fig.canvas.draw_idle()
    plt.axis('off')
            
    a_slider.on_changed(update)
    
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])    
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')   
    
    playax = plt.axes([0.4, 0.025, 0.18, 0.04])    
    playButton = Button(playax, 'Play animation', color=axcolor, hovercolor='0.975')   
    
    def reset(event):
        print("Reset requested...")
        a_slider.reset()

    def draw_update(curr):
        the_mesh = slice_df.iloc[curr,:].mesh
        fvals = the_mesh
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))
        
        ax.cla()
        ax3.cla()
        ax.set_title("Heatmaps of the %s perturbations for constant x "%title + fft_str)
        ax.set_xlabel("Y coordinate")
        ax.set_ylabel("Z coordinate")
        if use_contour:
            ax.contourf(fvals)
        else:
            ax.imshow(fvals)
        ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
        
    def play(event):
        print("Play animation...")
        anim = animation.FuncAnimation(fig, draw_update, frames=fvals.shape[0],interval=50, repeat=False)
        fig.canvas.draw()
        
    button.on_clicked(reset)
    playButton.on_clicked(play)
    
    resetax._button = button
    playax._button = playButton
    
    return fvals

def import_energy(energy_file,sep="SEPARATOR",strict=True):
    '''
    Import energy of field in z axis, which has been stored directly as a mesh. 
    The x/y values have been averaged out for each z.
    '''
    file = open(energy_file,'r')
    PE_list = []
    KE_list = []
    GE_list = []
    
    lres = file.readline().strip()
    if lres==sep:
        raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
    res = int(lres)
        
    l1 = file.readline().strip()    
    l2 = file.readline().strip()
    l3 = file.readline().strip()
    l4 = file.readline().strip()
    l5 = file.readline().strip()
    counter = 0
    while l1!='' and l2!='' and l3!='' and l4!='' and l5!='':        
        counter +=1
        if l1!=sep:
            if strict:
                raise ValueError("Expected separator in second line. Instead found '"+l1[:1000]+"' at line %i."%counter) 
            else:
                print("WARNING: skipping at line %i due to misplaced separator..."%counter)
            while l1!='' and l1!=sep:
                l1 = file.readline().strip()
            l2 = file.readline().strip()
            l3 = file.readline().strip()
        if l2==sep or l3==sep or l4==sep or l5==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        mesh_PE = list(map(float, l3.split()))
        mesh_KE = list(map(float, l4.split()))
        mesh_GE = list(map(float, l5.split()))
        
        mesh_PE = np.array(mesh_PE)#.reshape(res,res).T
        mesh_KE = np.array(mesh_KE)
        mesh_GE = np.array(mesh_GE)
        
        PE_list.append([a,mesh_PE])
        KE_list.append([a,mesh_KE])
        GE_list.append([a,mesh_GE])
        
        l1 = file.readline().strip()
        l2 = file.readline().strip()
        l3 = file.readline().strip()
        l4 = file.readline().strip()
        l5 = file.readline().strip()
    
    df_PE = pd.DataFrame(PE_list, columns=['a','mesh'])
    df_KE = pd.DataFrame(KE_list, columns=['a','mesh'])
    df_GE = pd.DataFrame(GE_list, columns=['a','mesh'])
    
    file.close()
    return df_PE,df_KE,df_GE

def plot_energy(edf,a_ind=0,use_FFT=False,use_contour=False,title='',use_log=False):
    fig = plt.figure()
    plt.subplots_adjust(left=0.25, bottom=0.25)    
    ax = fig.add_subplot(1,1,1)
    ax.set_title("Evolution of the %s Energy for constant x/y"%title)
    ax.set_xlabel("Z coordinate")
    ax.set_ylabel("%s Energy"%title)
    fvals = edf.iloc[a_ind,:].mesh
    ymin = edf.mesh.apply(np.min).min()
    ymax = edf.mesh.apply(np.max).max()
    if use_log:
        ax.set_yscale('log')
    if use_FFT:
        fvals = np.log(np.absolute(np.fft.fft(fvals)))
        mesh_series = edf.mesh.apply(lambda x: np.log(np.absolute(np.fft.fft(x))))
        ymin = mesh_series.apply(np.min).min()
        ymax = mesh_series.apply(np.max).max()
    print("Y limits: ",ymin, ymax)
    ax.plot(np.arange(0,fvals.shape[0]),fvals)
    frame_box = ax.text(.9,.9,"a_index: %i\na: %.1f"%(a_ind,edf.a[a_ind]),horizontalAlignment='center',transform=ax.transAxes)
    
    axcolor = 'lightgoldenrodyellow'
    playax = plt.axes([0.4, 0.025, 0.18, 0.04])    
    playButton = Button(playax, 'Play animation', color=axcolor, hovercolor='0.975')
    
    ax.set_xlim([0,fvals.shape[0]])
    ax.set_ylim([ymin,ymax])
    
    
    slid_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    a_slider = Slider(slid_ax, 'A', 0, edf.shape[0]-1, valinit=a_ind, valstep=1)
    
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])    
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')   

    def update(val):
        av = int(a_slider.val) 
        fvals = edf.iloc[av,:].mesh
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))
        ax.cla()
        if use_log:
            ax.set_yscale('log')
        ax.set_xlim([0,fvals.shape[0]])
        ax.set_ylim([ymin,ymax])
        ax.set_title("Evolution of the %s Energy for constant x/y"%title)
        ax.set_xlabel("Z coordinate")
        ax.set_ylabel("%s Energy"%title)
        
        ax.plot(np.arange(0,fvals.shape[0]),fvals)

    a_slider.on_changed(update)
    def reset(event):
        print("Reset requested...")
        a_slider.reset()
    button.on_clicked(reset)
    resetax._button = button
    
    
    def draw_update(curr):
        fvals = edf.iloc[curr,:].mesh
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft(fvals)))
        
        ax.cla()
        if use_log:
            ax.set_yscale('log')

        ax.set_xlim([0,fvals.shape[0]])
        ax.set_ylim([ymin,ymax])
        ax.set_title("Evolution of the %s Energy for averaged x/y"%title)
        ax.set_xlabel("Z coordinate")
        ax.set_ylabel("%s Energy"%title)
        
        ax.plot(np.arange(0,fvals.shape[0]),fvals)
        frame_box = ax.text(.9,.9,"a_index: %i\na: %.1f"%(curr,edf.a[curr]),horizontalAlignment='center',transform=ax.transAxes)
        
    def play(event):
        print("Play animation...")
        print(edf.index)
        anim = animation.FuncAnimation(fig, draw_update, frames=edf.index,interval=50, repeat=False)
        fig.canvas.draw()
    
    playButton.on_clicked(play)
    playax._button = playButton
    return fvals

def slice_to_gif(slice_df,a_ind=0,a_max=-1,use_FFT=False,use_vmin=False, use_contour=False, title='',fps=-1,out='animated_slice',save=True,t_max=-1):
    the_mesh = slice_df.iloc[a_ind,:].mesh
    fvals = the_mesh
    fft_str=''
    if a_max ==-1 or a_max > slice_df.shape[0]:
        a_max = slice_df.shape[0]
    jump = 1
    if fps==-1 and t_max==-1:
        jump = 1
        fps = 24
    elif fps==-1:
        fps = 24
        jump = ceil(a_max/fps/t_max)
    elif t_max==-1:
        jump = 1
       
    if use_FFT:
        fft_str='\n(log of modulus of FFT)'
        fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    vmin, vmax = slice_df.mesh[a_ind].min(),slice_df.mesh[a_ind].max()
    if use_vmin:
        for j in slice_df.index[a_ind:a_max:jump]:
            v1 = slice_df.mesh[j].min()
            v2 = slice_df.mesh[j].max()
            if v1<vmin: vmin=v1
            if v2>vmax: vmax=v2
        
    print("Vmin,vmax: ",vmin,vmax)
    bounds = np.linspace(vmin,vmax,1000)
    
    X, Y = np.meshgrid(np.linspace(1,fvals.shape[0],fvals.shape[0]),np.linspace(1,fvals.shape[1],fvals.shape[1]))
    fig,ax = plt.subplots()
    fig.set_figheight(6)
    fig.set_figwidth(7.5)
    #plt.tight_layout()
    ax.set_aspect(1)
    ax.set_title("Heatmaps of the %s perturbations for constant x up to %d"%(title,slice_df.a[a_max-1]) + fft_str)
    ax.set_xlabel("Y coordinate")
    ax.set_ylabel("Z coordinate")

    #ax.contourf(fvals)
    if use_contour:
        img_plot = ax.contourf(fvals, vmin=vmin, vmax=vmax)
    else:
        img_plot = ax.imshow(fvals,  vmin=vmin, vmax=vmax)
    if use_vmin:
        cb = fig.colorbar(ScalarMappable(norm = img_plot.norm, cmap = img_plot.cmap),
                 boundaries=bounds,
                 ticks=np.linspace(vmin,vmax,10))
        
    def update(i):
        the_mesh = slice_df.iloc[i,:].mesh
        fvals = the_mesh
        
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))        
        ax.clear()
        if not use_vmin:
            vmin_loc, vmax_loc = slice_df.mesh[i].min(),slice_df.mesh[i].max()
    
            if use_contour:
                img_plot = ax.contourf(fvals, vmin=vmin_loc, vmax=vmax_loc)
            else:
                img_plot = ax.imshow(fvals,  vmin=vmin_loc, vmax=vmax_loc)
        else:
            if use_contour:
                img_plot = ax.contourf(fvals, vmin=vmin, vmax=vmax)
            else:
                img_plot = ax.imshow(fvals,  vmin=vmin, vmax=vmax)

        
        ax.set_title("Heatmaps of the %s perturbations for constant x up to %d"%(title,slice_df.a[a_max-1]) + fft_str)
        ax.set_xlabel("Y coordinate")
        ax.set_ylabel("Z coordinate")
        time_text = ax.text(.05,.87,'',transform = ax.transAxes, bbox=props,color='orange')
        time_text.set_text(time_template%(slice_df['a'][i],i))
        
    time_template = "Metric $a$= %.7f \nRow: %i "
    props = dict(boxstyle='round', facecolor='blue', alpha=0.3)
    time_text = ax.text(.05,.87,'',transform = ax.transAxes, bbox=props,color='orange')
    time_text.set_text(time_template%(slice_df['a'][a_ind],a_ind))

        
    anim = animation.FuncAnimation(fig, update,frames=slice_df.index[a_ind:a_max:jump], repeat=False)
    writer = animation.PillowWriter(fps=fps)
    if save:
        anim.save(out+'.gif', writer=writer)
    return anim


def perts_to_gif(slice_df,a_ind=0,a_max=-1,use_FFT=False,use_vmin=True, log_scale=False, use_contour=False, title='',fps=-1,out='animated_slice',save=True,t_max=-1):
    the_mesh = slice_df.iloc[a_ind,:].mesh
    fvals = np.abs(the_mesh - the_mesh.mean())
    fft_str=''
    if a_max ==-1 or a_max > slice_df.shape[0]:
        a_max = slice_df.shape[0]
    jump = 1
    if fps==-1 and t_max==-1:
        jump = 1
        fps = 24
    elif fps==-1:
        fps = 24
        jump = ceil(a_max/fps/t_max)
    elif t_max==-1:
        jump = 1
    if use_FFT:
        fft_str='\n(log of modulus of FFT)'
        fvals = np.log(np.absolute(np.fft.fft2(fvals)))
    vmin = np.abs(slice_df.mesh[a_ind]-slice_df.mesh[a_ind].mean()).min()
    vmax = np.abs(slice_df.mesh[a_ind]-slice_df.mesh[a_ind].mean()).max()
    if use_vmin:
        for j in slice_df.index[a_ind:a_max:jump]:
            v1 = np.abs(slice_df.mesh[j]-slice_df.mesh[j].mean()).min()
            v2 = np.abs(slice_df.mesh[j]-slice_df.mesh[j].mean()).max()
            if v1<vmin: vmin=v1
            if v2>vmax: vmax=v2
    if log_scale:
        norm = matplotlib.colors.LogNorm()
        fvals = np.log10(fvals)
        vmin, vmax = np.log10(vmin), np.log10(vmax)
    else:
        norm = matplotlib.colors.Normalize()

        
    print("Vmin,vmax: ",vmin,vmax)
    bounds = np.linspace(vmin,vmax,1000)
    
    X, Y = np.meshgrid(np.linspace(1,fvals.shape[0],fvals.shape[0]),np.linspace(1,fvals.shape[1],fvals.shape[1]))
    fig,ax = plt.subplots()
    fig.set_figheight(6)
    fig.set_figwidth(7.5)
    #plt.tight_layout()
    ax.set_aspect(1)
    ax.set_title("Heatmaps of the %s perturbations for constant x up to %d"%(title,slice_df.a[a_max-1]) + fft_str)
    ax.set_xlabel("Y coordinate")
    ax.set_ylabel("Z coordinate")

    #ax.contourf(fvals)
    if use_contour:
        img_plot = ax.contourf(fvals, vmin=vmin, vmax=vmax)
    else:
        img_plot = ax.imshow(fvals,  vmin=vmin, vmax=vmax)
    if use_vmin:
        cb = fig.colorbar(ScalarMappable(norm = img_plot.norm, cmap = img_plot.cmap),
                 boundaries=bounds,
                 ticks=np.linspace(vmin,vmax,10))
        
    def update(i):
        the_mesh = slice_df.iloc[i,:].mesh
        fvals = np.abs(the_mesh - the_mesh.mean())
        
        if use_FFT:
            fvals = np.log(np.absolute(np.fft.fft2(fvals)))        
        ax.clear()
        if log_scale:
            fvals = np.log(fvals)

        if not use_vmin:
            vmin_loc = np.abs(slice_df.mesh[i]-slice_df.mesh[i].mean()).min()
            vmax_loc = np.abs(slice_df.mesh[i]-slice_df.mesh[i].mean()).max()
            if log_scale:
                vmin_loc, vmax_loc = np.log10(vmin_loc), np.log10(vmax_loc)

                
            if use_contour:
                img_plot = ax.contourf(fvals, vmin=vmin_loc, vmax=vmax_loc)
            else:
                img_plot = ax.imshow(fvals,  vmin=vmin_loc, vmax=vmax_loc)
        else:
            if use_contour:
                img_plot = ax.contourf(fvals, vmin=vmin, vmax=vmax)
            else:
                img_plot = ax.imshow(fvals,  vmin=vmin, vmax=vmax)
        

        
        ax.set_title("Heatmaps of the %s perturbations for constant x up to %d"%(title,slice_df.a[a_max-1]) + fft_str)
        ax.set_xlabel("Y coordinate")
        ax.set_ylabel("Z coordinate")
        time_text = ax.text(.05,.87,'',transform = ax.transAxes, bbox=props,color='orange')
        time_text.set_text(time_template%(slice_df['a'][i],i))
        
    time_template = "Metric $a$= %.7f \nRow: %i "
    props = dict(boxstyle='round', facecolor='blue', alpha=0.3)
    time_text = ax.text(.05,.87,'',transform = ax.transAxes, bbox=props,color='orange')
    time_text.set_text(time_template%(slice_df['a'][a_ind],a_ind))

        
    anim = animation.FuncAnimation(fig, update,frames=slice_df.index[a_ind:a_max:jump], repeat=False)
    writer = animation.PillowWriter(fps=fps)
    if save:
        anim.save(out+'.gif', writer=writer)
    return anim



def import_a_mesh(fields_file,sep="SEPARATOR",selected=1,up_to=0,outfile="",save=True):
    if outfile=="":
        outfile="placeholder"
    
    # Selected: int, choose which field state to plot.
    if up_to==0: up_to=selected
    
    file = open(fields_file,'r')
    mesh_pair = []
    
    lres = file.readline().strip()
    if lres==sep:
        raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
    res = int(lres)
        
    l1 = file.readline().strip()    
    l2 = file.readline().strip()
    l3 = file.readline().strip()
    counter = 1
    while l1!='' and counter<=up_to:        
        if l1!=sep:
            raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
        if l2==sep or l3==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        
        if counter>=selected and counter<=up_to:
            print("Saving mesh %i..."%counter)
            mesh = list(map(float, l3.split()))
            mesh = np.array(mesh).reshape(res,res,res).T
            
            mesh_pair = [a, mesh]
            if save:
                np.save(outfile+"_"+str(counter)+".npy",mesh_pair[1])
            print("Done saving.")
        else:
            print("Counter: ",counter,"Skipping...")
        
        counter += 1
        l1 = file.readline().strip()
        l2 = file.readline().strip()
        l3 = file.readline().strip()
        
    
    file.close()
    return mesh_pair



    
def mesh_to_gif(fields_file,outfile='out_mesh',sep="SEPARATOR",start_from=1,up_to=1,use_perts=False,scientific=False,cutoff=-2.3):
    file = open(fields_file,'r')
    mesh_pair = []
    
    
    #Find vmin, vmax
    vmin,vmax = 0,0
    if scientific:
        lres = file.readline().strip()
        if lres==sep:
            raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
        res = int(lres)
            
        l1 = file.readline().strip()    
        l2 = file.readline().strip()
        l3 = file.readline().strip()
        counter = 1
        
        print("Going through file to find vmin, vmax...")
        while l1!='' and counter<=up_to:        
            if l1!=sep:
                raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
            if l2==sep or l3==sep:
                raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
            a = float(l2)
            
            if counter>=start_from:
                print("Scanning field %i..."%counter)
                mesh = list(map(float, l3.split()))
                mesh = np.array(mesh).reshape(res,res,res).T
                if not use_perts:
                    if counter==1:
                        vmin = mesh.min()
                        vmax = mesh.max()
                    else:
                        v1,v2 = mesh.min(),mesh.max()
                        if v1<vmin: vmin=v1
                        if v2>vmax: vmax=v2
                else:
                    if counter==1:
                        avg_mesh = mesh.mean()
                        vmin = np.log10(np.abs((mesh-avg_mesh))).min()
                        vmax = np.log10(np.abs((mesh-avg_mesh))).max()
                    else:
                        avg_mesh = mesh.mean()
                        v1,v2 = np.log10(np.abs((mesh-avg_mesh))).min(), np.log10(np.abs((mesh-avg_mesh))).max()
                        if v1<vmin: vmin=v1
                        if v2>vmax: vmax=v2
            else:
                print("Skipping field %i."%counter)
            counter += 1
            l1 = file.readline().strip()
            l2 = file.readline().strip()
            l3 = file.readline().strip()
        print("Vmin, vmax: ",vmin,vmax)
        file.close()
        file = open(fields_file,'r')

        
    lres = file.readline().strip()
    if lres==sep:
        raise ValueError("Separator found in first line. Expected resolution. Check that file format corresponds to requirements for 'import_mesh'.")
    res = int(lres)
        
    l1 = file.readline().strip()    
    l2 = file.readline().strip()
    l3 = file.readline().strip()
    counter = 1     
    
    while l1!='' and counter<=up_to:        
        if l1!=sep:
            raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
        if l2==sep or l3==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        
        if counter<start_from:
            print("Skipping mesh %i..."%counter)
        else:
            print("Saving mesh %i..."%counter,end=' ')
            mesh = list(map(float, l3.split()))
            mesh = np.array(mesh).reshape(res,res,res).T
                
            mesh_pair = [a, mesh]
                
            print("Done saving. Plotting...")
            if scientific:
                mesh_plot(mesh_pair, out=outfile, use_perts=use_perts,counter=counter,cutoff=cutoff,vs=[cutoff,vmax])
            else:
                mesh_plot(mesh_pair, out=outfile, use_perts=use_perts,counter=counter,cutoff=cutoff)
        counter += 1
        l1 = file.readline().strip()
        l2 = file.readline().strip()
        l3 = file.readline().strip()
    
    
    file.close()
    return mesh_pair
    
def normalise(x,vmin,vmax):
    return (x-vmin)/(vmax-vmin)

def unnormalise(x,vmin,vmax):
    return (vmax-vmin)*x+vmin

def mesh_plot(mesh_pair,out='out_mesh',counter=1,use_perts=False,cutoff=-2.5,vs=[],show=False):
    a_val = mesh_pair[0]
    the_mesh = mesh_pair[1]
    n= the_mesh.shape[0]
    
    coords = np.linspace(1,n,n)
    X,Y,Z = np.meshgrid(coords, coords, coords)
    xs = X.reshape(1,n**3)[0]
    ys = Y.reshape(1,n**3)[0]
    zs = Z.reshape(1,n**3)[0]
    ts = the_mesh.reshape(1,n**3)[0]
    
    mesh_avg = the_mesh.mean()
    mesh_perts = np.log10(np.abs(the_mesh - mesh_avg)).reshape(1,n**3)[0]
    
    if use_perts:
        ts = mesh_perts

    x_keep, y_keep, z_keep, t_keep = [], [], [], [] 
    for i in range(len(xs)):
        if mesh_perts[i]>cutoff:
            x_keep.append(xs[i])
            y_keep.append(ys[i])
            z_keep.append(zs[i])
            t_keep.append(ts[i])
            
    fig = plt.figure()    
    fig.set_figwidth(10)
    fig.set_figheight(8)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(r"Perturbations produced in the Palatini Higgs model (cutoff: $\log_{10}{|\chi/\chi_0|} >%.2f$)($a=%.3f$)"%(cutoff,mesh_pair[0]))
    
    print("Remaining points after cutoff: ",len(x_keep))        
    if len(x_keep)==0:
        fig.savefig(out + '%i.jpeg'%counter)
        plt.close(fig)
        return
    

    # Plot at different transparences
    scats=[]
    if len(vs)!=2:
        v1,v2 = min(t_keep), max(t_keep)
    else:
        v1,v2 = vs[0],vs[1]
    #print(v1,v2)
    bounds = [0.33,0.66]
    alphas = [0.25,0.5,0.75]
    x_seps, y_seps, z_seps, t_seps = [], [], [], []
    for _ in alphas:
        x_seps.append([])
        y_seps.append([])
        z_seps.append([])
        t_seps.append([])
    for i in range(len(t_keep)):
        found = False
        for b in range(len(bounds)):
            if t_keep[i]<unnormalise(bounds[b],v1,v2):
                x_seps[b].append(x_keep[i])
                y_seps[b].append(y_keep[i])
                z_seps[b].append(z_keep[i])
                t_seps[b].append(t_keep[i])
                found = True
                break
        if not found:
            x_seps[len(bounds)].append(x_keep[i])
            y_seps[len(bounds)].append(y_keep[i])
            z_seps[len(bounds)].append(z_keep[i])
            t_seps[len(bounds)].append(t_keep[i])
    for t in t_seps:
        #print("Length: ", len(t))
        pass
    
    if len(vs)!=2:
        #print("Scattering at alpha: ",end=' ')            
        for a in range(len(alphas)):
            #print(alphas[a],'...',end=' ')
            scats.append(ax.scatter(x_seps[a], y_seps[a], z_seps[a], c=t_seps[a], alpha=alphas[a]))
        fig.colorbar(scats[0],ax=ax)
        #print('\n')
    elif vs[0]<=vs[1]:
        #print("Scattering at alpha: ",end=' ')            
        for a in range(len(alphas)):
            #print(alphas[a],'...',end=' ')
            scats.append(ax.scatter(x_seps[a], y_seps[a], z_seps[a], c=t_seps[a], alpha=alphas[a], vmin=v1,vmax=v2))
        fig.colorbar(scats[0],ax=ax)
        #print('\n')

    
    fig.savefig(out + '%i.jpeg'%counter)
    if not show:
        plt.close(fig)



def imgs_to_gif(template_path,img_type='.jpeg', indices=[],fps=1):
    if len(indices)==0:
        print("No indices given")
        return
    images = []
    for i in sorted(indices):
        file_path = template_path+str(i)+img_type
        if os.path.exists(file_path):
            print("Processing index: ",i)
            images.append(imageio.imread(file_path))
    imageio.mimsave(template_path+'GIF'+ str(min(indices)) + '_' + str(max(indices)) + '.gif', images,fps=fps)

def import_metric(path):
    file = open(path,'r')
    field_list = []

    counter = 1
    res = int(file.readline().strip())
    l1 = file.readline().strip()    
    l2 = file.readline().strip()    

    while l1!='' and l2!='':        
        mesh = list(map(float, l2.split()))
        
        mesh = np.array(mesh).reshape(res,res,res).T
        
        field_list.append([float(l1),mesh])
        counter +=1
        l1 = file.readline().strip()
        l2 = file.readline().strip()
    
    mesh_df = pd.DataFrame(field_list, columns=['a','mesh'])
    file.close()
    return mesh_df    

def split_metric(path,screen_file,mode='h'):
    modes = {'h': 'metric_h',
             'p': 'metric_p',
             'f1': 'q_k',
             'f2': 'qdot_k',
             'fd':'fld_density',
             'fp': 'fld_pressure'}
    destination = trim_name(screen_file)+'_MESHES'
    if not os.path.isdir(destination):
        print("Making new folder...")
        os.mkdir(destination)
    else:
        print("Folder already found!")
    file = open(path,'r')
    
    counter = 1
    lres = file.readline()
    l1 = file.readline()   
    l2 = file.readline()

    while l1!='' and l2!='':                
        print("\rProcessing mesh: ",counter,end='')
        output = open(destination + r'\\%s_%i.log' %( modes[mode], counter), 'w')
        output.write(lres)
        output.write(l1)
        output.write(l2)
        
        counter +=1
        l1 = file.readline()
        l2 = file.readline()
        output.close()
    print("Done! Closing input file and exiting...")
    file.close()
    
def metric_folder_to_gif(first_file,out="",dim=2,inds=[],fps=2,use_perts=True,cutoff=3,scientific=False):  
    folder_ind = first_file.rindex('\\')  
    folder = first_file[:folder_ind]
    extension_ind = first_file.rindex('.')
    extension = first_file[extension_ind:]
    
    v_ind = extension_ind-1
    try:
        while type(int(first_file[v_ind]))==int:
            v_ind -=1
    except ValueError:
        pass
    
    first = int(first_file[v_ind+1:extension_ind])
    path = first_file[: v_ind+1] # Note: subtract 1 in order to remove counter
    
    counter = first
    plt.ioff()

    images = []
    if dim==2:
        if scientific:
            vmin,vmax = 0,0
            while os.path.exists(path + str(counter) + extension):
                dest = path + str(counter) + extension
                print("\rScanning for minimum/maximum in image %i"%counter,dest,end='; ')
                temp_df = import_metric(dest)
                if use_perts:
                    if counter==first:
                        vmin = np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())).min()
                        vmax = np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())).max()
                    v1 = np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())).min()
                    v2 = np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())).max()
                    if v1<vmin: vmin = v1s
                    if v2>vmax: vmax = v2
                else:
                    if counter==first:
                        vmin = temp_df.mesh[0][1,:,:].min()
                        vmax = temp_df.mesh[0][1,:,:].max()
                    v1 = temp_df.mesh[0][1,:,:].min()
                    v2 = temp_df.mesh[0][1,:,:].max()
                    if v1<vmin: vmin = v1
                    if v2>vmax: vmax = v2
                counter += 1
            counter = first
        while os.path.exists(path + str(counter) + extension):
            dest = path + str(counter) + extension
            print("\rProcessing:",dest,end='; ')
            temp_df = import_metric(dest)
            print("\rPlotting image %i..."%counter,end='')
            fig,ax = plt.subplots()
            if scientific:
                if use_perts:
                    img = ax.imshow(np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())),vmin=vmin,vmax=vmax)
                else:
                    img = ax.imshow(temp_df.mesh[0][1,:,:],vmin=vmin,vmax=vmax)
            else:
                if use_perts:
                    img = ax.imshow(np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())))
                else:
                    img = ax.imshow(temp_df.mesh[0][1,:,:])
            ax.set_title("$a=%.4f$"%temp_df.a[0])
            fig.colorbar(img,ax=ax)
            fig.savefig(folder + '\\dump.jpeg')
            plt.close(fig)
            print("\rAdding image %i..."%counter,end ='')
            images.append(imageio.imread(folder + '\\dump.jpeg'))
            counter += 1
    elif dim==3:
        while os.path.exists(path + str(counter) + extension):
            dest = path + str(counter) + extension
            print("\rProcessing in 3D:",dest,end='; ')
            temp_df = import_metric(dest)
            print("Plotting in 3D...",end=' ')
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            if use_perts:
                ax.imshow(np.log10(np.abs(temp_df.mesh[0][1,:,:]-temp_df.mesh[0][1,:,:].mean())))
            else:
                ax.imshow(temp_df.mesh[0][1,:,:])
            ax.set_title("$a=%.4f$"%temp_df.a[0])
            fig.savefig(folder + '\\dump.jpeg')
            plt.close(fig)
            print("Adding image...",end =' ')
            images.append(imageio.imread(folder + '\\dump.jpeg'))
            counter += 1
        


    print("\nDone! Saving to GIF...")
    if out=="":
        v=1
        gif_out = path+'_GIF'+ str(first) + '_' + str(counter-1) + '.gif'
        while os.path.exists(gif_out):
            gif_out = path+'_GIF'+ str(first) + '_' + str(counter-1) +'_v%i'%v + '.gif'
            v += 1
        imageio.mimsave(gif_out, images,fps=fps)
    else:
        imageio.mimsave(folder+'\\'+out+'_GIF'+ str(first) + '_' + str(counter-1) + '.gif', images,fps=fps)
        


path2 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\dump_test_screen.log"
some_path = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\dump_full_metric.log"
#%% Animating the results
def animate_potential(data,plt_obj=0,t_max=25,fps=25, Mpl=1):
    # Set plot object
    if plt_obj ==0:
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            print("ERROR: unable to find matplotlib (?!). Please check it is installed. Exiting function...")
            return
        plt_obj = plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.set_figwidth(8)
        fig.set_figheight(6)
    
    elif type(plt_obj)==list:
        if len(plt_obj)!=2:
            print("ERROR: list provided; expected 2 arguments (fig,ax), instead found %i."%len(plt_obj))
            return
        fig, ax = plt_obj
    else:
        print("ERROR: did not recognise plotting object.")
        return
        
    #t_max given in seconds
    anim_running = True
    
    #Work out the interval and fps to optimise the animation
    interval = 1000/fps  #Time per frame: 1 second / frames per second
    step = ceil(data.shape[0]/t_max / fps) #How large a step between rows to ensure fps and t_max are met
    print("Animation: fps: ",fps,", milliseconds: ",interval,", step: ",step,", data.size: ",data[::step].shape[0])
    
    phi_max = np.abs(data['mean1']).max() * 1.02 #Provide 2% extra margin
    xs = np.linspace(0,phi_max,1000)
    ys = palatiniV(xs)
    v_phi_interp = scipy.interpolate.PchipInterpolator(xs, ys, extrapolate=False)
    
    def onClick(event):
        nonlocal anim_running
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True
    '''        
    xp0 = data['mean1'][0]
    yp0 = palatiniV(xp0)
    def anim_data():
        nonlocal step
        for i in range(0,data.shape[0],step):
            #print(i)
            yield [data['mean1'][i],data['a'][i],i]
            
    def animate(anim_data):
        xp = anim_data[0]
        yp = palatiniV(xp)
        p1.set_data(xp,yp)
        time_text.set_text(time_template%(anim_data[1],anim_data[2]))
    
    p0, = ax.plot(xs,ys)
    p1, = ax.plot(xp0,yp0,'ro')
    time_template = "Metric $a$= %.7f \nRow: %i "
    time_text = ax.text(0.05,0.87,'', transform=ax.transAxes)
    time_text.set_text(time_template%(data['a'][0],0))
    
    anim = animation.FuncAnimation(fig,animate,anim_data,interval=interval,repeat=False,save_count=50)
    fig.canvas.mpl_connect('button_press_event', onClick)
    '''    

    phis = np.abs(data['mean1'])
    phi2s = np.sqrt(data['mean1']**2 + data['rms1']**2)
    
    xp1 = phis[0]
    yp1 = palatiniV(xp1)
    xp2 = phi2s[0]
    yp2 = palatiniV(xp2)
    
    def anim_data():
        nonlocal step
        for i in range(0,phis.shape[0],step):
            #print(i)
            yield [phis[i],phi2s[i],data['a'][i],i]
            
    def animate(anim_data):
        xp1 = anim_data[0]
        xp2 = anim_data[1]
        yp1 = palatiniV(xp1)
        yp2 = palatiniV(xp2)
        p1.set_data(xp1,yp1)
        p2.set_data(xp2,yp2)
        
        time_text.set_text(time_template%(anim_data[2],anim_data[3]))
    
    p0, = ax.plot(xs,ys)
    p1, = ax.plot(xp1,yp1,'ro')
    p2, = ax.plot(xp2,yp2,'go')
    time_template = "Metric $a$= %.7f \nRow: %i "
    time_text = ax.text(0.05,0.87,'', transform=ax.transAxes)
    time_text.set_text(time_template%(data['a'][0],0))
    ax.set_title(r"Potential for $\langle \phi \rangle(t)$ and $\sqrt{\langle \phi^2 \rangle(t)}$")
    ax.set_xlabel("Field value (units of $M_{P}=$%d)"%Mpl)
    ax.set_ylabel("Potential ($M_{P}=$%d)"%Mpl)
    print(anim_running)
    anim = animation.FuncAnimation(fig,animate,anim_data,interval=interval,repeat=True,save_count=50)
    fig.canvas.mpl_connect('button_press_event', onClick)
        


#%% Convert datafiles to CSV

def df_to_csv(*args, path='', transposed=False):
    dfs = []
    for a in args:
        dfs.append(a.copy())
    print("Length of dfs: ",len(dfs))
    if path=='':
        path = "test.csv"
    if transposed:
        for i in range(len(dfs)):
            dfs[i] = dfs[i].T
    if len(dfs)==1:
        dfs[0].to_csv(path+'_csv.txt', index=False)
    else:
        i = 1
        for d in dfs:
            d.to_csv(path+'_'+str(i)+'_csv.txt', index=False)
            i+=1
        
#%% Get information from the sim_settings file.
def sim_settings(filefile,delim=':'):
    # Relies on finding the colon ':' in any given line.
    sim_file = trim_name(filefile) + '_sim_settings.log'
    file = open(sim_file,'r')
    descriptors = []
    values = []
    
    l1 = file.readline().strip()
    while l1!='':
        pos = l1.find(delim)
        if pos!=-1:
            descriptors.append(l1[:pos])
            val_str = l1[pos+1:].strip()
            if int(float(val_str)) == float(val_str):
                values.append(int(float(val_str)))
            else:
                values.append(float(val_str))
        else:
            descriptors.append(l1)
            values.append('n/a')
        l1 = file.readline().strip()
    return (descriptors,values)
 


#%% Floquet data from Python
flfile = r"C:\Users\James\Downloads\floquet_data.txt"
import re, pprint
'''
with open(flfile,'r') as f:
    fl_data = f.read()
results = re.findall(r'\{\}', fl_data)   
res2 = re.findall('\[[^\]]*\]|\{[^\}]*\)|\"[^\"]*\"|\S+',fl_data) 
pprint.pprint(results)
  '''
def make_pure_qkdf(pw_data1,data,L=64,LH=1,process=True,Mpl=1):
    dim1 = (data.shape[0]-1)
    dim2 = (pw_data1.shape[0]-1)
    chkpt_ratio = (dim1-dim1%dim2) / dim2

    H0 = data.h[0]
    kmin = data.a[0] *H0
    
    # 2pi / n * j  * n * h /LH ~= 
    ks = k_list(pw_data1,L=L) * L *H0/ LH
    a_list = pw_data1['a']

    qks = pw_data1.drop('a',axis=1)
    etot = 3* data['h'][::int(chkpt_ratio)]**2 * Mpl**2 * (data['omega'][::int(chkpt_ratio)]+1)
    etot.reset_index(drop=True,inplace=True)

    vals = np.sqrt(qks)    
    if process:
        
        vals = vals.multiply(a_list,axis=0) 
        vals /= ks**2.5 
        
        vals = vals.multiply(etot**(1/2),axis=0) * np.sqrt(2*np.pi**2)
    vals.columns = ks
    vals['a'] = pw_data1['a']
    vals['etot'] = etot
    
    return vals

def plot_qk_t(pw_data1):
    fig,ax = plt.subplots()
    colors = np.flip(cm.magma(np.linspace(0.1,1,pw_data1.shape[1])),axis=0)
    for i in range(len(pw_data1.columns)-1):
        ax.plot(pw_data1.a,pw_data1.iloc[:,i]/pw_data1.iloc[0,i],c=colors[i])
    ax.set_yscale('log')
    ax.set_title(r"Perturbations in the field $\varphi_k$")

def plot_qk(pw_data1,data,L=64,use_metric=False,jump=5,LH=0.1,save=False,filefile='',use_xlog=True,plot=True,Mpl=1):
    if plot:
        fig, ax = plt.subplots()
    dim1 = (data.shape[0]-1)
    dim2 = (pw_data1.shape[0]-1)
    chkpt_ratio = (dim1-dim1%dim2) / dim2
    
    H0 = data.h[0]
    kmin = data.a[0] *H0
    
    # 2pi / n * j  * n * h /LH ~= 
    ks = k_list(pw_data1,L=L) * L *H0/ LH
    a_list = pw_data1['a']
    
    qks = pw_data1.drop('a',axis=1)
    t_mat = np.ones(qks.shape) * (chkpt_ratio * qks.index.values.reshape(-1,1) +1) # 1_matrix * time
    
    
    
    
    #data_indices = [data['a'][data['a'].eq(pw_data.a[i])].index[j] for i in range(pw_data1.shape[0]) for j in range(data['a'][data['a'].eq(pw_data.a[i])].index)]
    if data.columns[2] == 'ef':
        etot = data['ef'][::int(chkpt_ratio)]
    else:
        etot = 3* data['h'][::int(chkpt_ratio)]**2 * Mpl**2 * (data['omega'][::int(chkpt_ratio)]+1)
    etot.reset_index(drop=True,inplace=True)
    
    #Colours 
    colors = np.flip(cm.viridis(np.linspace(0.1,1,t_mat.shape[0])),axis=0)
    
    #Plot
    vals = np.sqrt(qks) 
    #vals /= t_mat
    p = 1
    vals = vals.multiply(a_list**p,axis=0) 
    vals /= ks**2.5 
    
    j = 1/2 *1
    vals = vals.multiply(etot**j,axis=0) * np.sqrt(2*np.pi**2)
    
    m = 1
    if use_xlog and plot:
        ax.set_xscale('log')
        ks /= kmin
    else:
        ks = np.log10(ks/kmin)
    if not plot:
        return vals,ks,t_mat
        
    for i in range(0,vals.shape[0]-1,jump):
        #print(i)
        #ax.plot(ks/kmin, np.log(vals.iloc[i,:])/a_list[i]**m,color=colors[i])
        ax.plot(ks, np.log(vals.iloc[i,:])/t_mat[i]**m,color=colors[i])
    #ax.set_yscale('log')
    ax.set_ylabel(r'$C \cdot \mu_k$')
    ax.set_xlabel(r'$k/k_{min}$')
    ax.set_title(r"Growth index $\mu_k$ for LH=%.3f"%LH)
    if save:
        fig.savefig(trim_file_name(filefile)+"_mu_k.png")
    return vals,ks,t_mat

def plot_spectogram(arr_2d):
    pw_da_data

#%% Main
if __name__=="__main__":
    path = r'D:\Physics\MPhys Project\DatasetArcive\Remote tests\rtanh-math-test12_screen.log'
    data = import_screen(path)
    pw1, pw2 = import_pw(trim_name(path) + '_pw_1.log')
    df_to_csv(pw1,pw2,path=trim_file_name(path))
    #animate_potential(data,t_max=50)
    

    