# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:42:02 2020

@author: James
"""
import numpy as np
import pandas as pd
import scipy.special as sp
from scipy import interpolate, signal,misc


def phi_amplitude(data):
    '''Deprecated: was designed for a time-dependent inflaton amplitude. 
    However, a comoving amplitude has ~ constant amplitude
    '''
    phis = np.abs(data['mean1'])
    ts = np.linspace(data.index.min(),data.index.max(),data.shape[0]+1)
    pks = signal.find_peaks(phis)[0]
    
    tp = ts.iloc[pks]
    phip = phis.iloc[pks]
    spl = interpolate.interp1d(tp,phip,kind='cubic')
    
    within_range = np.where((ts>=np.min(tp)) & (ts<=np.max(tp)))[0]
    tinterp = ts.iloc[within_range]
    phi_amp = spl(xinterp)
    
    return tinterp,phi_amp,spl

def cmvphi_amp(data):
    comoving_phis = np.abs(data['mean1'] * data['a'])
    #ts = np.linspace(data.index.min(),data.index.max(),data.shape[0]+1)
    pks = signal.find_peaks(comoving_phis)[0]
    cmv_phi_pks = comoving_phis.iloc[pks]
    
    return cmv_phi_pks.mean() * 2
    

def xvar(data,conformal=True,amp=False):
    if conformal:
        if type(amp)==bool: raise UserWarning("What the fuck are you doing with the amplitude?")
        xvar = np.sqrt(lam) * amp * np.array(data.index)
    else:   
        ts = np.linspace(data.index.min(),data.index.max(),data.shape[0]+1)
        const = (6 * lam * Mpl**2 / np.pi)**0.25
        xvar = const * np.sqrt(ts)
    return xvar

def init_momenta(v,dvdf,field0, use_PM=False):
    result = -dvdf(field0,*args) / np.sqrt( 3 * v(field0, *args))
    if use_PM:
        result *= (1/np.sqrt(8*np.pi)) **2

def palatiniV(field0,l=10**-4,xi=-1,nstar=50,use_PM=False):
    if xi==-1:
        xi = 3.8*10**6 * nstar**2 * l
    a = l / (4 * xi**2)
    b = np.sqrt(xi)
    eta = np.sqrt(8*np.pi)
    
    if use_PM:
        a /= eta**4
        b *= eta
    
    potential = a * np.tanh(b*field0)**4
    return potential

def palatiniDV(field0,l=10**-4,xi=-1,nstar=50,use_PM=False):
    if xi==-1:
        xi = 3.8*10**6 * nstar**2 * l
    a = l / (4 * xi**2)
    b = np.sqrt(xi)
    eta = np.sqrt(8*np.pi)
    
    if use_PM:
        a /= eta**4
        b *= eta

    dvdf = 4*a*b * np.tanh(b*field0)**3 / np.cosh(b*field0)**2
    return dvdf
    
    
    
    
    



