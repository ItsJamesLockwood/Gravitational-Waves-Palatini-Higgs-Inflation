# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:08:31 2021

@author: James
"""

from physUtils import *
from plotUtils import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import math as maths
import numpy as np

#%% Define variables


#%% Define files
r_math = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\rtanh-math-test%i%s_screen.log"
t_16 = r"D:\Physics\MPhys Project\gw-local-repo\HLatticeV2.0\data\tanh-16-lh1-run%i_screen.log"

tv = 5
t_16f = t_16 % tv

t_math_v = 10
t_math_s = ''
r_math_f = r_math% (t_math_v, t_math_s)

eliot1 = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\noMpl_tach_1_1024_slices_screen.log"
eliot2 = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\long_sim_1_64_GW_screen.log"
eliot3 = r"D:\Physics\MPhys Project\DatasetArcive\long_sim_1_64_GW_screen.log"
eliot4 = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\long_1_64_GW_80_screen.log"
eliot5 = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\1_128_80_long_slice_real_screen.log"

superlong = r"C:\Users\James\Downloads\super_long_lf4_screen.log"


analysis_file = eliot4
#analysis_file = t_16f
pw_field_number = 1
form = 'log'
#%% Get sim settings variables
try:
    strs, vals = sim_settings(analysis_file)
    RESOLUTION = vals[6] 
    SKIP = int(vals[11]/vals[10] )
    BOXSIZE_H = vals[1] 
    Mpl = vals[14]
    CONF = vals[7]
    sim_vals = {'skip':SKIP, 'lh': BOXSIZE_H, 'mpl': Mpl, 'resolution':RESOLUTION,'conformal':CONF}
    for k,v in sim_vals.items():
        print(k,': ',v)
except FileNotFoundError:
    print("Did not update settings...")
    BOXSIZE_H = 15
    SKIP = 5 
    CONF = 1
    RESOLUTION = 64
    sim_vals = {}
#%% Placholder analysis function
'''
!Plot 1: Energies 
!Plot 2: phi^2 vs phi
Plot 3: potential animation
!Plot 4: field oscillation w/ perturbations
Plot 5: perturbations integrated
Plot 6: loglog plot of \eta
    
'''

def placeholder(data,pw1,pw2,ns,rows=[],sim_vals={}, save_panel=False,save_plots=False,path=analysis_file,truncate=0,ts=False):
    # Try first to import gravitational wave files
    GW_files_found = True
    try:
        gw1, gw2 = import_GW(trim_name(analysis_file) + '_GW.log')
    except FileNotFoundError:
        print("No GW files found. Proceeding without...")
        GW_files_found = False
        
    png_name = trim_file_name(path)
    if truncate!=0:
        png_name += '_trunc' + str(truncate)
        truncate = np.searchsorted(data['a'],truncate)
    else:
        truncate = data.shape[0]-1
    
    #Default sim values
    if len(sim_vals)==0:
        SKIP = 10
        LH = 15
        RESOLUTION = 64
        Mpl = 1024
        conformal = 1
    else:
        SKIP = sim_vals['skip']
        LH = sim_vals['lh']
        RESOLUTION = sim_vals['resolution']
        Mpl = sim_vals['mpl']
        conf = sim_vals['conformal']
    #Default colors    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    prop_colors = prop_cycle.by_key()['color']
    
    fig, ax = plt.subplots(2,3)
    fig.set_figwidth(13)
    fig.set_figheight(9.5)
    fig.suptitle(png_name)
    plt.subplots_adjust(top=.93)
    
    #Subplot 0,0
    ax[0,0].set_title("Inflaton field evolution")
    ax[0,0].plot(data['a'][:truncate], data['mean1'][:truncate])
    if ts:
        ax[0,0].plot(data['a'][:truncate], data['mean1'][:truncate],'r.')
    c2 = np.flip(cm.magma(np.linspace(0,1,len(ns.columns)-1)),axis=0)
    ax_twin =  ax[0,0].twinx()
    for j in range(len(ns.columns)-1):
        ax_twin.plot(ns['a'][:truncate],np.log(pw_data1[:truncate].iloc[:,j]/pw_data1[:truncate].iloc[0,j]),
                     color=c2[j],
                     alpha=1)
    
    #Subplot 1,0
    etot = 3* data['h'][:truncate]**2 * Mpl**2 * (data['omega'][:truncate]+1)
    ax[1,0].set_yscale('log')
    ax[1,0].set_title('Reproduction of Fig.1: ratios of energies')
    ax[1,0].set_xlabel('a')
    ax[1,0].set_ylabel('$\log_{10}(|E|/E_{tot})$')
    ax[1,0].plot(data['a'][:truncate],data['pratio'][:truncate],linestyle='dashed',label='Potential energy')
    ax[1,0].plot(data['a'][:truncate],data['kratio'][:truncate],linestyle='dashed',label='Kinetic energy')
    ax[1,0].plot(data['a'][:truncate],data['gratio'][:truncate],'b',label='Gradient of field energy')
    #ax[1,0].plot(data['a'][:truncate],data['gratio'][:truncate]*etot,color='orange',label='Gradient energy (net)')
    ax[1,0].plot(data['a'][:truncate],etot/etot[0],color='purple',label='Total energy $E^{tot}/E_0^{tot}$')
    
    c2 = np.flip(cm.magma(np.linspace(0,1,len(ns.columns)-1)),axis=0)
    for j in range(len(ns.columns)-1):
        #pass
        ax[1,0].plot(ns['a'][:truncate],ns[:truncate].iloc[:,j]/ns[:truncate].drop('a',axis=1).min().min(),color=c2[j])
    ax[1,0].plot(data['a'][:truncate],abs(1/(data['omega'][:truncate]+1)-1),linestyle='dashed',label='Fractional energy noises')
    ax[1,0].legend()
        
    #Subplot 0,1
    n_pal = n_palatini(pw1,pw2, data, L=RESOLUTION, LH=LH, skip=SKIP)
    int_perts = integrate_perturbation(n_pal,Mpl=Mpl)
    ax[0,1].plot(n_pal['a'],int_perts)
    ax[0,1].set_yscale('log')
    ax[0,1].set_title(r"Integrated perturbations $\rho_{pert}/\rho_{BG}$")
    ax[0,1].set_xlabel("$a$")
    ax[0,1].set_ylabel(r"$\log_{10} (\rho_{pert}/\rho_{BG} )$")

    
    #Subplot 1,1
    ax_twin = ax[1,1].twinx()
    #p1, = ax[1,1].plot(data['a'], (data['mean1']**2 + data['rms1']**2),label=r"$\langle \phi^{2} \rangle $")
    p2, = ax_twin.plot(data['a'],np.abs(data['mean1']),color=prop_colors[1],label=r"$\langle \phi \rangle$")
    p3, = ax_twin.plot(data['a'],np.sqrt(data['mean1']**2 + data['rms1']**2),color=prop_colors[2],label=r"$\sqrt{\langle \phi^2 \rangle}$")
    
    lns = [p2,p3]
    ax[1,1].legend(handles=lns)
    ax[1,1].set_title("Comparison of the average field and squared field values")
        
    
    #Subplot 0,2 : animated potential (stolen from the plotUtils function)
    anim_running = True
    
    #Work out the interval and fps to optimise the animation
    t_max= 25
    fps =3
    interval = maths.ceil(1000/fps)  #Time per frame: 1 second / frames per second
    step = ceil(data.shape[0]/t_max / fps) #How large a step between rows to ensure fps and t_max are met
    print("Animation: fps: ",fps,", milliseconds: ",interval,", step: ",step,", data.size: ",data[::step].shape[0])
    
    phi_max = np.abs(data['mean1']).max() * 1.02 #Provide 2% extra margin
    xs = np.linspace(0,phi_max,1000)
    ys = palatiniV(xs)
    v_phi_interp = scipy.interpolate.PchipInterpolator(xs, ys, extrapolate=False)
    
    def onClick(event):
        nonlocal anim_running
        if anim_running:
            #anim.pause()
            anim.event_source.stop()
            anim_running = False
        else:
            #anim.resume()
            anim.event_source.start()
            anim_running = True
    
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
    
    p0, = ax[0,2].plot(xs,ys)
    p1, = ax[0,2].plot(xp1,yp1,'ro')
    p2, = ax[0,2].plot(xp2,yp2,'go')
    time_template = "Metric $a$= %.7f \nRow: %i "
    time_text = ax[0,2].text(0.05,0.87,'', transform=ax[0,2].transAxes)
    time_text.set_text(time_template%(data['a'][0],0))
    ax[0,2].set_title(r"Potential for $\langle \phi \rangle(t)$ and $\sqrt{\langle \phi^2 \rangle(t)}$")
    ax[0,2].set_xlabel("Field value (units of $M_{P}=$%d)"%Mpl)
    ax[0,2].set_ylabel("Potential ($M_{P}=$%d)"%Mpl)
    print(anim_running)
    anim = animation.FuncAnimation(fig,animate,anim_data,interval=interval,repeat=True,save_count=50)
    
    #fig.canvas.mpl_connect('button_press_event', onClick)    
    #Subplot 1,2
    ax[1,2].plot(data.index,data['a'])
    ax[1,2].set_yscale('log')
    ax[1,2].set_xscale('log')          
    if conf==1:
        ax[1,2].set_xlabel(r'Conformal time $\eta$')
    else:
        ax[1,2].set_xlabel(r"Physical time $t$")
    ax[1,2].set_ylabel("Metric $a$")
    ax[1,2].set_title("Log-log plot of the metric v. %s time"%('conformal' if conf else 'physical'))
    
    fig.tight_layout()
    
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
    return anim

def plot_loglog_time(data,conf = CONF):
    fig, ax = plt.subplots()
    ax.plot(data.index,data['a'])
    ax.set_yscale('log')
    ax.set_xscale('log')
    if conf==1:
        ax.set_xlabel(r'Conformal time $\eta$')
    else:
        ax.set_xlabel(r"Physical time $t$")
    ax.set_ylabel("Metric $a$")
    ax.set_title("Log-log plot of the metric v. %s time"%('conformal' if conf else 'physical'))
    
    fig.tight_layout()

#%% import relevant files
data = import_screen(analysis_file)
pw_data1, pw_data2 = import_pw(trim_name(analysis_file) + '_pw_%i.%s'%(pw_field_number,form))

n_df = n_k(pw_data1,pw_data2,L=RESOLUTION)

my_rows = range(1,n_df.shape[0],2)

#%% Plot the panel
anim = placeholder(data,pw_data1,pw_data2,n_df,sim_vals=sim_vals)


