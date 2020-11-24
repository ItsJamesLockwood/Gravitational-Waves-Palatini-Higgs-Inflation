# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:02:21 2020

@author: James
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from plotUtils import *
from physUtils import *
#TODO: obtain these values from the data files
#TODO: check conformal time used correctly


#%% Define variables 
global lam, Mpl, g
lam = 10**-2
Mpl = 1024
Nstar= 50
xi = 3.8*10**6 * lam * Nstar**2
g = np.sqrt(3 * lam) #TODO: ?

#wd = os.getcwd()
wd_path = "D:/Physics/MPhys Project/gw-local-repo/HLatticeV2.0/"
os.chdir(wd_path)
print("Current working dir: ",os.getcwd())

def  trim_name(file):
    ind = file.index('_screen')
    return file[5:ind]

#%% File management
file_name = "data/run_with_GW_17_11_screen.log"
GW_file = "data/run_with_GW_17_11_GW.log"
eliot_file = "data/Dataset-Friday_screen.log"
higgs_file = "data/higgs-vev-run1_screen.log"
tanhfile = "data/higgs-tanh4-run1_screen.log"

filefile = eliot_file
pw_field_number = 2 #Choose which field spectrum to plot (start: 1)
my_img_name =  trim_name(filefile) + '_img'



#%% Function definitions


def plot_pw_t(df1,save_img=True,img_name='df1',trim=10):
    dfs = df1.iloc[:,:trim]
    a = df1['a']
    fig,ax = plt.subplots()
    for c in dfs.columns:
        ax.plot(a,dfs[c])
    plt.yscale('log')
    if save_img==True:
        fig.savefig("pw_img_"+img_name)
    fig.show()

def plot_pw_k(df1,save_img=True,img_name='df1',trim=10):
    if trim>=df1.shape[1]: trim = -1
    dfs = df1.T.iloc[:-1,:trim]
    fig,ax = plt.subplots()
    for c in dfs.columns:
        ax.plot(dfs[c])
    plt.yscale('log')
    if save_img==True:
        fig.savefig("pw_img_"+img_name)
    fig.show()

def plot_fig1(data,error=True,save_img=True,img_name="img",truncate=0):
    plt.yscale('log')
    plt.title('Reproduction of Fig.1: ratios of energies')
    plt.xlabel('a')
    plt.ylabel('$\log_{10}(|E|/E_{tot})$')
    plt.plot(data['a'],data['pratio'],linestyle='dashed',label='Potential energy')
    plt.plot(data['a'],data['kratio'],linestyle='dashed',label='Kinetic energy')
    plt.plot(data['a'],data['gratio'],'b',label='Gradient of field energy')
    plt.legend()
    
    if save_img==True:
        plt.savefig(img_name)
    if error==True:
        plt.plot(data['a'][truncate:],abs(1/(data['omega'][truncate:]+1)-1),linestyle='dashed',label='Fractional energy noises')
        plt.legend()
        if save_img==True: plt.savefig(img_name + "_with_error")
    
    plt.show()
def plot_fig2_tchakev(data):
    plt.yscale('log')
    plt.xscale('log')
    plt.title("Reproduction of Tchakev's Fig2")
    #fftchi = np.fft.fft(data['mean2'])
    #plt.plot(fftchi)
    m_chi2 = 1**2
    m = 1

    a2 = data['a']**2
    print(a2.size)
    k = np.logspace(-1,2,a2.size)

    eps2 = k**2/a2 + m_chi2
    phi0 = data['a']**(-3/2)
    g = 1 #1.095 * 10**(-6)

    beta2 = np.exp(np.pi * eps2 / g / phi0 /m)
    NMAX = 30
    print(beta2[:NMAX])
    plt.plot(k[:NMAX],beta2[:NMAX])

    plt.show()
 
def k_list(pw_data1,L=64):
    ki = 2 * np.pi / L
    #Note: We only want to create ks for the spectrum values and not for the a column
    #Save as np.array for extra functionality
    return np.array([ki * i for i in range(1,pw_data1.shape[1])]) 

def n_k(pw_data1,pw_data2,L=64):
    pw_a = pw_data1.drop(['a'],axis=1)
    pw_b = pw_data2.drop(['a'],axis=1)
    a_list = pw_data1['a'] # By construction, this is the same as pw_data2['a']
    
    ks = k_list(pw_data1,L=L)
    
    #Retrieve actual field eignemode values
    fk_df = pw_a.multiply(a_list**2,axis=0) / ks**5
    fkdot_df = pw_b / ks**3
    #fk_df = pw_a
    #fkdot_df = pw_b
    
    #Some pretty cool element-wise multiplication between vector ks and dataframes fk_df and fkdot_df!
    ns = 1/(2*ks) * (ks**2 * fk_df + fkdot_df)
    #Store ks in the column headers
    ns.columns = ks
    #Add a values back for plotting purposes
    ns['a'] = a_list
    return ns

def plot_n_k(ns,rows=[]):
    if rows==[]:
        rows.append(int(ns.shape[0]))
                    
    for r in rows:
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns.columns[:-1]) 
        ys = ns.iloc[r,:-1] * xs**4
        plt.plot(xs,ys,label=ns['a'][r])
    #plt.yscale('log')
    plt.legend(title='Spectrum at $a=$',loc='lower right')
    plt.show()
    
def plot_fig6(ns,rows=[],vlines=False):
    if rows==[]:
        rows.append(int(ns.shape[0])/2)
    colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
    for j in range(len(rows)):
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns.columns[:-1]) / (2 * np.pi)
        ys = ns.iloc[rows[j],:-1] * (2 * xs**4)
        plt.plot(xs,ys,label=ns['a'][rows[j]], color=colors[j])
    plt.yscale('log')
    
    if vlines: 
        vert_lines = np.array([1.25, 1.8, 3.35])/ (2*np.pi) 
        for vl in vert_lines:
            plt.axvline(x = vl)
    
    plt.legend(title='Spectrum at $a=$',loc='lower right')
    plt.show()

#%% Main

data = import_screen(filefile)
pw_data1, pw_data2 = import_pw('data/'+trim_name(filefile) + '_pw_%i.log'%pw_field_number)

n_df = n_k(pw_data1,pw_data2)
plot_fig6(n_df, rows=[1,10,20,30,40,50,60,70,80,90])
#plot_gw(pw_data1,trim=2,save_img=False)



#plt.plot(pw_data1.iloc[4,:-1]**2)
#plot_pw_k(pw_data1,save_img=True,trim=10)
#plt.plot(np.abs(data['a']*data['mean1']))
#plt.yscale('log')
#plt.plot(data['a'][truncate:],1/(data['omega'][truncate:]+1)-1)
#plt.plot(np.fft.fft(np.cos(np.linspace(0,100,1000))))
#plt.plot(sp.ellipj(np.linspace(0,10,100),1/2**0.5)[1])
#plt.show()
#plot_fig1(data,error=True,img_name=my_img_name)
#plot_fig2_tchakev(data)

#df1,df2 = import_GW(GW_file)
#plot_gw(df1,img_name='df1')
#plot_gw(df2,img_name='df2')

