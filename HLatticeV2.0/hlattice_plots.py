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
from itertools import chain

#%% Define variables 
global lam, Mpl, g
lam = 10**-2
Mpl = 1024
Nstar= 50
xi = 3.8*10**6 * lam * Nstar**2
g = np.sqrt(3 * lam) #TODO: ?
L = 64

#wd = os.getcwd()
wd_path = "D:/Physics/MPhys Project/gw-local-repo/HLatticeV2.0/"
os.chdir(wd_path)
print("Current working dir: ",os.getcwd())

def  trim_name(file):
    ind = file.index('_screen')
    return file[:ind]

#%% File management
file_name = "data/run_with_GW_17_11_screen.log"
GW_file = "data/run_with_GW_17_11_GW.log"
eliot_file = "data/lf4-std-run1_screen.log"
higgs_file = "data/higgs-vev-run1_screen.log"
tanhfile = "data/higgs-tanh4-run1_screen.log"
lf4_tkachev1 = "data/lf4-tkachev-coupling1_screen.log"
lf4_tkachev2 = "data/lf4-tkachev-coupling2_screen.log"
lf4_tkachev3 = "data/lf4-tkachev-coupling3_screen.log"
file128 = "data/lf4-std-h2-128-run1_screen.log"
fref = r"D:\Physics\MPhys Project\DatasetArcive\lf4-128-H2-M1-LH30_screen.log"
unstable1 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r1_screen.log"
unstable2 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r2_screen.log"
unstable3 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r3_screen.log"
unstable4 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r4_screen.log"
unstable5 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r5_screen.log"
unstable6 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\test-exit-r6_screen.log"

filefile = unstable4
pw_field_number =1 #Choose which field spectrum to plot (start: 1)
form = 'log'
rows=[1,10,20,30,40,50,60,70,80,90]
rsmall = [1,5,10,15,20,25,30,35,40,45]
tk1 = [1,20,50,75,100] + list(range(2,19,3))
tk2 = range(0,22,3)
rmid = [1,10,90] + list(range(30,60,3))
my_rows = rows
tk_rows = tk2
my_rows = sorted(my_rows)
my_img_name =  trim_name(filefile) + '_img'
save_img = False



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
    fk_df = pw_a / ks**5
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

def n_k_red(pw_data1,pw_data2,L=64):
    pw_a = pw_data1.drop(['a'],axis=1)
    pw_b = pw_data2.drop(['a'],axis=1)
    a_list = pw_data1['a'] # By construction, this is the same as pw_data2['a']

    ks = k_list(pw_data1,L=L)
    #Retrieve actual field eignemode values
    fk_df = pw_a.multiply(a_list**2,axis=0) 
    fkdot_df = pw_b

    #Some pretty cool element-wise multiplication between vector ks and dataframes fk_df and fkdot_df!
    ns = 1/2 * (fk_df + fkdot_df)
    #Store ks in the column headers
    ns.columns = ks
    #Add a values back for plotting purposes
    ns['a'] = a_list
    return ns

def plot_n_t(ns,cols=[],save_img=False,img_name='Fig3_rubio',data=[]):
    if cols==[]:
        cols.append(int(ns.shape[1]/2))
    colors = np.flip(cm.magma(np.linspace(0,1,len(cols))),axis=0)
 
    for c in range(len(cols)):
        print('C:',c)
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns['a']) 
        ks = np.array(ns.columns[:-1])
        ys = ns.iloc[:,cols[c]] * ks[cols[c]]**4
        plt.plot(xs,ys,label=ks[cols[c]],color=colors[c])
    if type(data)==pd.DataFrame:
        plt.plot(data['a'],10**3*np.exp(data['mean1']),label='Field oscillation',linestyle='dashed')
    plt.yscale('log')
    plt.xlabel('$a(t)$')
    plt.ylabel('$n_k(t)$')
    plt.title("Occupation number $n_k(t)$ for different eigenmodes $k$ (Fig.3, Rubio)\n"+trim_name(filefile)+'_field_%i'%pw_field_number)
    plt.legend(title='Spectrum at $k=$',loc='lower right')
    if save_img:
        plt.savefig(img_name + '_fig3_Rubio_f%i'%pw_field_number)
        
    plt.show()
    


def plot_n_k(ns,rows=[]):
    if rows==[]:
        rows.append(int(ns.shape[0]/2))
                    
    for r in rows:
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns.columns[:-1]) 
        ys = ns.iloc[r,:-1] * xs**4
        plt.plot(xs,ys,label=ns['a'][r])
    #plt.yscale('log')
    plt.legend(title='Spectrum at $a=$',loc='lower right')
    plt.show()
 
def plot_n_k_red(ns,rows=[]):
    if rows==[]:
        rows.append(int(ns.shape[0]/2))                   
    for r in rows:
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns.columns[:-1]) 
        ys = ns.iloc[r,:-1]
        plt.plot(xs,ys,label=ns['a'][r])
    plt.yscale('log')
    plt.legend(title='Spectrum at $a=$',loc='lower right')
    plt.show()
     
def plot_tkachev2(ns,rows=[],save_img=False,img_name="Tkachev2"):
    if rows==[]:
        rows.append(int(ns.shape[0])/2)
        colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
    else:
        colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
    for j in range(len(rows)):
        #Retrieve k values from the column headers, discarding the 'a' in the last column
        xs = np.array(ns.columns[:-1]) / (2 * np.pi)
        ys = ns.iloc[rows[j],:-1] * 2
        plt.plot(xs,ys,label="$a=$%f"%ns['a'][rows[j]], color=colors[j])
    plt.yscale('log')
    plt.xscale('log')
    xl1 = np.linspace(min(xs),max(xs),1000)

    plt.plot(xl1, xl1**(-3/2)*10**2,linestyle='dashed',label='~$k^{-3/2}$')
    plt.plot(xl1, xl1**(-1)*10**4,linestyle='dashed',label='~$k^{-1}$')
    plt.legend(loc='lower right')
    plt.title("Occuptation number $n_k$ (Fig.2, Tkachev)")
    plt.xlabel('k')
    plt.ylabel('$n_k$')
    if save_img==True:
        plt.savefig(img_name + '_fig2_tkachev_f%i'%pw_field_number)

    plt.show()

def plot_fig6(ns,rows=[],vlines=False,save_img=False,img_name="Fig 6"):
    if rows==[]:
        rows.append(int(ns.shape[0])/2)
        colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
    else:
        colors = np.flip(cm.magma(np.linspace(0,1,len(rows))),axis=0)
        #colors = np.flip(cm.magma(np.array(rows)/max(rows)),axis=0)
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
    plt.title("Occupation number $n_k$ at different scale factors $a$ (Fig.6, Z.Huang)\n"+trim_name(filefile)+'_field_%i'%pw_field_number)
    plt.xlabel(r"$k\Delta / 2 \pi$")
    plt.ylabel(r"$ k^4\; n_k /\; 2 \pi^2\; \rho}$")
    if save_img==True:
        plt.savefig(img_name + '_fig6_field_%i'%pw_field_number)

    plt.show()

#%% Main

data = import_screen(filefile)
pw_data1, pw_data2 = import_pw(trim_name(filefile) + '_pw_%i.%s'%(pw_field_number,form))

n_df = n_k(pw_data1,pw_data2,L=L)
nred_df  = n_k_red(pw_data1,pw_data2,L=L)
print(data['a'])
#plot_n_t(n_df,cols=[1,4,5,20,30],save_img=save_img,img_name=my_img_name,data=data)
my_rows = np.searchsorted(n_df['a'],my_rows)
my_rows = range(1,n_df.shape[0],2)
plot_fig6(n_df, rows=my_rows, save_img=save_img,img_name=my_img_name)
tk_rows = sorted(tk_rows)
#plot_tkachev2(n_df,rows=tk_rows,save_img=save_img,img_name=my_img_name)
#plot_gw(pw_data1,trim=2,save_img=False)
loc_min = 0
loc_max = None

#plt.plot(data['a'][loc_min:loc_max]**.5,data['mean1'][loc_min:loc_max])
#plt.plot(data.diff()['a'][loc_min:loc_max])

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

#%% Temp test of omega
'''templ = [3.0818286343172027E-004,
         3.2058427750028849E-004,
         3.4025032220321951E-004,
         3.6601186544975902E-004,
         3.9668311112553732E-004,
         4.3121769069237337E-004,
         4.6876252409028659E-004,
         5.0865146449070083E-004,
         5.5037507372443560E-004,
         5.9354657132985465E-004,
         6.3787204360222561E-004,
         6.8312689655197897E-004,
         7.2913810124543754E-004,
         7.7577109104069513E-004,
         8.2292016733869411E-004,
         8.7050147561197339E-004,
         9.1844784140641727E-004,
         9.6670494818077574E-004]
tl = np.array(templ)
xs = np.arange(1,len(tl)+1)
plt.yscale('log')
plt.xscale('log')
plt.plot(xs,tl,'ro')
'''