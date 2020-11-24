# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
wd = os.getcwd()
wd_path = r"C:\Users\eliot\OneDrive\Documents\physics 4th year\Mphys project\HLatticeV2.0.tar\HLatticeV2.0\data"
os.chdir(wd_path)
print("Current working dir: ",os.getcwd())
GW_file = r"Dataset-Friday_pw_1.log"
my_img_name = "fig2-Thermalization-Paper"
L=10**8
m_chi=0
global psdf


def import_GW(GW_file):
    file = open(GW_file,'r')
    global psf
    global a
    global psdf
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
    #print(psf[1])
                                  
    return df1,df2 
df1,df2 = import_GW(GW_file)
ks = [2 * np.pi / L * (i+1) for i in range(len(psf[1]))]
omega_k=[ks[n]**2+m_chi**2 for n in range(len(ks))]
N_q_all=[]#list of list: each list represent N_k at different time k

for i in range(len(a)):
    psfi=psf[i]
    psdfi=psdf[i]
    phi_k_i=[2*psfi[l]*a[i]**2/ks[l]**5 for l in range(len(ks))]
    dot_phi_k_i=[2*psdfi[p]/ks[p]**3 for p in range(len(ks))]
    N_qi=[1/2/ks[k]*(dot_phi_k_i[k]/2/omega_k[k]+omega_k[k]**0.5*phi_k_i[k]/2) for k in range(len(ks))]
    N_q_all.append(N_qi)
#print(N_q_all[1])
plt.yscale('log')
for i in range(len(N_q_all)):
    plt.plot(ks,np.array(ks)**4 * np.array(N_q_all[i]))

plt.show()

#np.array(ks)**4 * np.array(N_q_all[i])

  
   

