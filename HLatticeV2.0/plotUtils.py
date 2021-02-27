# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:38:40 2020

@author: James
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

def import_screen(file_name,print_logs=False):
    #old_col_names = [chr(i) for i in range(ord('a'),ord('a')+10)]
    col_names = ['a',
    			'h',
    			'omega',
    			'pratio',
    			'kratio',
    			'gratio',
    			'mean1',
    			'mean2',
    			'rms1',
    			'rms2']
    
    data = pd.read_csv(file_name,delim_whitespace=True,skiprows=1,names=col_names,index_col=False)
    try:
    	data = data[~data.iloc[:,0].str.contains("[a-zA-Z]").fillna(False)]
    except AttributeError:
    	pass
    
    data = data.apply(pd.to_numeric)
    
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
    
    df1['a'] = pd.Series(a)
    df2['a'] = pd.Series(a)
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

def import_slice(slice_file,sep="SEPARATOR"):
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
    while l1!='':        
        if l1!=sep:
            raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
        if l2==sep or l3==sep:
            raise ValueError("Expected either metric or field mesh in lines 2 and 3. Instead found the separator.")
        a = float(l2)
        mesh = list(map(float, l3.split()))
        mesh = np.array(mesh).reshape(res,res).T
        
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
    cbar = fig.colorbar(pl3d, ax=ax3)
    
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

def import_energy(energy_file,sep="SEPARATOR"):
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

    while l1!='' and l2!='' and l3!='' and l4!='' and l5!='':        
        if l1!=sep:
            raise ValueError("Expected separator in second line. Instead found '"+l1+"'.")    
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

    