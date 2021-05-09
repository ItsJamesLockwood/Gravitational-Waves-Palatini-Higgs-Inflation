# -*- coding: utf-8 -*-
"""
Created on Sun May  9 01:43:29 2021

@author: James
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
import logging
from plotUtils import *
from physUtils import *


    
def metric_mission_control(data,ns,pw_data1,rows=[],error=True,save_panel=False,save_plots=False,path='',truncate=0,ts=False,xi=3.8*10**6*50**2*10**-4):    
    strs, vals = sim_settings(path)
    RESOLUTION = vals[6] 
    # Try first to import gravitational wave files
    GW_files_found = True
    try:
        gw1, gw2 = import_GW(trim_name(path) + '_GW.log')
    except FileNotFoundError:
        print("No GW files found. Proceeding without...")
        GW_files_found = False
        
    png_name = trim_file_name(path)
    if truncate!=0:
        png_name += '_trunc' + str(truncate)
        truncate = np.searchsorted(data['a'],truncate)
    else:
        truncate = data.shape[0]-1
    #Default colors    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    prop_colors = prop_cycle.by_key()['color']
    
    fig, ax = plt.subplots(2,3)
    fig.set_figwidth(13)
    fig.set_figheight(9.5)
    fig.suptitle(png_name)
    plt.subplots_adjust(top=.93)
    #Subplt 0,0
    ax[0,0].set_title("Inflaton field evolution")
    ax[0,0].plot(data['a'][:truncate], data['mean1'][:truncate])
    ax[0,0].axhline(1.05/np.sqrt(xi),linestyle='dashed',color='grey')
    ax[0,0].axhline(-1.05/np.sqrt(xi),linestyle='dashed',color='grey')
    
    if ts:
        ax[0,0].plot(data['a'][:truncate], data['mean1'][:truncate],'r.')
    c2 = np.flip(cm.magma(np.linspace(0,1,len(ns.columns)-1)),axis=0)
    ax_twin =  ax[0,0].twinx()
    for j in range(len(ns.columns)-1):
        ax_twin.plot(ns['a'][:truncate],np.log(pw_data1[:truncate].iloc[:,j]/pw_data1[:truncate].iloc[0,j]),
                     color=c2[j],
                     alpha=1)
    
        
    #Subplot 0,1
    if (not ts) and GW_files_found:
        global freqs
        freqs = gw1.drop('a',axis=1)
        gw_intensity = gw2.drop('a',axis=1)
        a_list = gw1['a']
        #Find those frequencies which vary little over the course of the simulation
        tolerance = 0.99
        i=0 
        
        # Find i such that the difference between the start and end is within the tolerance compared to the next band.
        while i<freqs.shape[1]-1 and (freqs.iloc[0,i]-freqs.iloc[-1,i]<0.5*(freqs.iloc[-1,i+1]-freqs.iloc[-1,i])):
            i += 1
            print(i)
        print("Tolerance %i %% -> truncation: "%(tolerance*100),i)
        #Ensure there are not too many labels
        if i<14:
            step = 1
        else: 
            step = ceil(freqs.shape[1]/14)
        colors = np.flip(cm.viridis(np.linspace(0,1,i)),axis=0)
        ax_twin =  ax[0,1].twinx()
        ax[0,1].patch.set_visible(False)
        ax_twin.set_zorder(0)
        ax[0,1].set_zorder(1)
    
        for j in range(len(freqs.iloc[:,:i].columns)):
            if j%step==0:
                ax[0,1].plot(a_list,gw_intensity.iloc[:,j],color=colors[j],label=freqs.iloc[:,j].median())
            else:
                ax[0,1].plot(a_list,gw_intensity.iloc[:,j],color=colors[j])
        c3 = np.flip(cm.magma(np.linspace(0,1,len(ns.columns)-1)),axis=0)
        for j in range(len(ns.columns)-1):
            ax_twin.plot(ns['a'][:truncate],np.log(pw_data1[:truncate].iloc[:,j]/pw_data1[:truncate].iloc[0,j]),
                     color=c3[j],
                     alpha=.6)        
        
        ax[0,1].set_yscale('log')
        ax[0,1].set_title("Evolution of the gravity waves modes (tolerance: %i%%)"%(100*tolerance))
        ax[0,1].legend()
        
    else:
        diffs = data['a'].diff()[1:]
        ax[0,1].set_title("Increment of a for each step j")
        ax[0,1].plot(data.index[1:][:truncate], diffs[:truncate])
        #ax[0,1].plot(data.a[1:][:truncate], diffs[:truncate])    
        ax[0,1].plot(data.index[1:][:truncate], diffs[:truncate],'r.')
        ewidth = (data.kratio-0.5*(4*data.pratio+2*data.gratio)).max() - (data.kratio-0.5*(4*data.pratio+2*data.gratio)).min()
        emid= (data.kratio-0.5*(4*data.pratio + 2*data.gratio)).min() + ewidth/2
        eavg = (data.kratio-0.5*(4*data.pratio + 2*data.gratio)).mean()
        dwidth, davg = diffs.max()-diffs.min(), diffs.mean()
        dmid = diffs.min() + dwidth/2
        print(dwidth,ewidth)
        print(davg,eavg)
        print(dmid,emid)
        energies = (data.kratio-0.5*(4*data.pratio+2*data.gratio)-emid)*dwidth/ewidth + dmid
        ax[0,1].plot(data.index[1:][:truncate], energies[1:][:truncate],alpha=0.5)
        ax[0,1].axhline(dmid-emid*dwidth/ewidth,linestyle='dashed',color='grey')
        #ax[0,1].set_yscale('log')
        #ax[0,1].set_xscale('log')
        xpower = np.linspace(data.index.min(),data.index.max(),100)
        ypower = xpower **(2/3) * 10**-5
        #ax[0,1].plot(xpower,ypower)
        
    #Subplot 1,0
    if not ts:
        if rows==[]:
            rows.append(int(ns.shape[0])/2)
            colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
        else:
            if len(rows)>14: 
                step = int(np.floor(len(rows)/10))
                rows = rows[::step]
            colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
        for j in range(len(rows)):
            xs = np.array(ns.columns[:-1]) / (2 * np.pi)
            ys = ns.iloc[rows[j],:-1] * (2 * xs**4)
            ax[1,0].plot(xs,ys,label=ns['a'][rows[j]], color=colors[j])
        ax[1,0].set_yscale('log')
        #plt.xscale('log')    
        if len(rows)<14:
            ax[1,0].legend(title='Spectrum at $a=$',loc='lower right')
        ax[1,0].set_title("Occupation number $n_k$")
        ax[1,0].set_xlabel(r"$k\Delta / 2 \pi$")
        ax[1,0].set_ylabel(r"$ k^4\; n_k /\; 2 \pi^2\; \rho}$")
      
    #Alternative Subplot 1,0
    else:
        dmetric = data['a'].diff()[1:]
        dfield = data['mean1'].diff()[1:]/dmetric
        ax[1,0].plot(data['a'][1:], dfield,label="Derivative from the field")
        ax[1,0].set_title("Derivative $d\phi/da$ of the inflaton")
        Etot = data['ef'] * data['eratio']
        #ax[1,0].set_yscale('log')
        #ax[1,0].set_xscale('log')
        
        fp = np.sqrt(2*RESOLUTION**3 * data['a']**0 * data['kratio']*Etot)
        ax[1,0].plot(data['a'],fp,label="$f_p$ from $E_{tot}$")
        ax[1,0].legend()
        
    
    #Subplot 1,1
    etot = data['ef'][:truncate]
    ax[1,1].set_yscale('log')
    ax[1,1].set_title('Reproduction of Fig.1: ratios of energies')
    ax[1,1].set_xlabel('a')
    ax[1,1].set_ylabel('$\log_{10}(|E|/E_{tot})$')
    ax[1,1].plot(data['a'][:truncate],data['pratio'][:truncate],linestyle='dashed',label='Potential energy')
    ax[1,1].plot(data['a'][:truncate],data['kratio'][:truncate],linestyle='dashed',label='Kinetic energy')
    ax[1,1].plot(data['a'][:truncate],data['gratio'][:truncate],'b',label='Gradient of field energy')
    #ax[1,1].plot(data['a'][:truncate],data['gratio'][:truncate]*etot,color='orange',label='Gradient energy (net)')
    ax[1,1].plot(data['a'][:truncate],etot/etot[0],color='purple',label='Total energy $E^{tot}/E_0^{tot}$')
    ax[1,1].legend(loc=1)
    
    c2 = np.flip(cm.magma(np.linspace(0,1,len(ns.columns)-1)),axis=0)
    for j in range(len(ns.columns)-1):
        #pass
        ax[1,1].plot(ns['a'][:truncate],ns[:truncate].iloc[:,j]/ns[:truncate].drop('a',axis=1).min().min(),color=c2[j])
    if error==True:
        ax[1,1].plot(data['a'][:truncate],abs(3*data['h'][:truncate]**2/(data['ef'][:truncate]*data['eratio'][:truncate])-1),linestyle='dashed',label='Fractional energy noises')
        ax[1,1].legend()
    
    #Subplot 0,2
    if GW_files_found:
        ax[0,2].set_yscale('log')    
        ax[0,2].set_xscale('log')
        
        frequencies = gw1.drop(['a'],axis=1)
        grav_intensity = gw2.drop(['a'],axis=1)
        a_list = gw1['a']
        if rows==[]:
            #rows.append(int(gw1.shape[0]/2))
            #colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
            rows.append(int(gw1.shape[0]/2))
            colors = np.flip(cm.magma(np.linspace(0,1,gw1.shape[0])),axis=0)
        else:
            #colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
            colors = np.flip(cm.magma(np.linspace(0,1,gw1.shape[0])),axis=0)
    
        print(rows)
        # Plot all of the spectra, but only label those in the selected rows
        for j in range(gw1.shape[0]):
            #Retrieve k values from the column headers, discarding the 'a' in the last column
            xs= frequencies.iloc[j,:]
            ys = grav_intensity.iloc[j,:]
            if j in rows:
                ax[0,2].plot(xs,ys,color=colors[j], label=a_list[j])
            else:
                ax[0,2].plot(xs,ys,color=colors[j])
                
        ax[0,2].legend(title='Spectrum at $a=$',loc='lower right')
        ax[0,2].set_title("Gravitational waves spectrum at different scale factors $a$ ")
        ax[0,2].set_xlabel(r"Frequency spectrum in Hz")
        ax[0,2].set_ylabel("(Omega_{gw}h^2)")
    else:
        ax[0,2].set_title("FILES NOT FOUND: Gravitational waves spectrum")
        
    #Subplot 1,2
    ax_twin = ax[1,2].twinx()
    p1, = ax[1,2].plot(data['a'], np.sqrt(data['mean1']**2 + data['rms1']**2),label=r"$\langle \phi^{2} \rangle $")
    p2, = ax_twin.plot(data['a'],np.abs(data['mean1']),color=prop_colors[1],label=r"$\langle \phi \rangle$")
    p3, = ax_twin.plot(data['a'],np.sqrt(data['mean1']**2 + data['rms1']**2),color=prop_colors[2],label=r"$\sqrt{\langle \phi^2 \rangle}$")
    
    lns = [p1,p2,p3]
    ax[1,2].legend(handles=lns)
    ax[1,2].set_title("Comparison of the average field and squared field values")
    ax[1,2].axhline(1.05/np.sqrt(xi),linestyle='dashed',color='grey')
    
    
    #End of function: save figures if requested
    if save_panel==True:
        print("Saving panel...")
        fig.savefig(png_name + '_panel.png')
    if save_plots==True:
        print("Saving individual plots...")
        #extent = ax[0,0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        extents = [ax[i,j].get_tightbbox(fig.canvas.renderer).transformed(fig.dpi_scale_trans.inverted()) for i in range(ax.shape[0]) for j in range(ax.shape[1])]
        
        # Pad the saved area by 10% in the x-direction and 20% in the y-direction
        padx = 1.05
        pady = 1.1
        fig.savefig(png_name + '_field.png', bbox_inches=extents[0].expanded(padx, pady))
        fig.savefig(png_name + '_increments.png', bbox_inches=extents[1].expanded(padx, pady))
        fig.savefig(png_name + '_n_k.png', bbox_inches=extents[2].expanded(padx, pady))
        fig.savefig(png_name + '_energies.png', bbox_inches=extents[3].expanded(padx, pady))
        if GW_files_found:
            fig.savefig(png_name + '_gw.png', bbox_inches=extents[4].expanded(padx, pady))
