# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:38:40 2020

@author: James
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
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
    return df1,df2


def import_fields(fields_file):
    file = open(fields_file,'r')
    field_list = []

    l1 = file.readline()
    while True:
        l2 = file.readline()
        if l1=='' or l2=='':
            break
        la = file.readline()
        lb = file.readline()
        while la!="SEPARATOR" and lb!="SEPARATOR":
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
    return field_df

def plot_fields(field_df,cond=['x', 1]):
    plane1 = field_df[field_df[cond[0]]==cond[1]]
    fvals = np.array(plane1.zs.values.tolist())
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

    ax.contourf(fvals)
    pl3d = ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
    cbar = fig.colorbar(pl3d, ax=ax3)
    
    axcolor = 'lightgoldenrodyellow'
    slid_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    x_slider = Slider(slid_ax, 'X', 1, plane1[anticond].max(), valinit=1, valstep=1)
    print(plane1[anticond].max())
    def update(val):
        xv = x_slider.val
        plane1 = field_df[field_df[cond[0]]==xv]
        print(plane1)
        fvals = np.array(plane1.zs.values.tolist())
        
        ax.cla()
        ax3.cla()
        ax.contourf(fvals)
        ax3.plot_surface(X=X, Y=Y, Z=fvals, cmap='YlGnBu_r')
        fig.canvas.draw_idle()
        
    x_slider.on_changed(update)
    
    
    
