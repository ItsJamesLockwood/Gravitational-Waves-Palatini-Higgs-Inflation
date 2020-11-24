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