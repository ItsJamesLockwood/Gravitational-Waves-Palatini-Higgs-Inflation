# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:42:02 2020

@author: James
"""
import numpy as np
import pandas as pd
import scipy.special as sp
from scipy import interpolate, signal,misc, integrate

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

def init_momenta(v,dvdf,field0, use_PM=True, **kwargs):
    result = -dvdf(field0,**kwargs,use_PM=use_PM) / np.sqrt( 3 * v(field0, **kwargs,use_PM=use_PM))
    if use_PM:
        result *= (1/np.sqrt(8*np.pi)) **2
    return result

def palatiniV(field0,use_PM=True,**kwargs):
    default = {'l':10**-4,'xi':-1,'nstar':50}
    ks = default
    for k in default.keys():
        if k in kwargs.keys():
            ks[k] = kwargs[k]
    #print(ks)

    if ks['xi']==-1:
        ks['xi'] = 3.8*10**6 * ks['nstar']**2 * ks['l']
    a = ks['l'] / (4 * ks['xi']**2)
    b = np.sqrt(ks['xi'])
    eta = np.sqrt(8*np.pi)

    
    if use_PM:
        a /= eta**4
        b *= eta
    
    potential = a * np.tanh(b*field0)**4
    return potential

def palatiniDV(field0,use_PM=True, **kwargs):
    default = {'l':10**-4,'xi':-1,'nstar':50}
    ks = default
    for k in default.keys():
        if k in kwargs.keys():
            ks[k] = kwargs[k]
    #print(ks)

    if ks['xi']==-1:
        ks['xi'] = 3.8*10**6 * ks['nstar']**2 * ks['l']
    a = ks['l'] / (4 * ks['xi']**2)
    b = np.sqrt(ks['xi'])
    eta = np.sqrt(8*np.pi)
    
    if use_PM:
        a /= eta**4
        b *= eta

    dvdf = 4*a*b * np.tanh(b*field0)**3 / np.cosh(b*field0)**2
    return dvdf

def palatiniD2V(field0,use_PM=True, **kwargs):
    default = {'l':10**-4,'xi':-1,'nstar':50}
    ks = default
    for k in default.keys():
        if k in kwargs.keys():
            ks[k] = kwargs[k]
    #print(ks)

    if ks['xi']==-1:
        ks['xi'] = 3.8*10**6 * ks['nstar']**2 * ks['l']
    a = ks['l'] / (4 * ks['xi']**2)
    b = np.sqrt(ks['xi'])
    eta = np.sqrt(8*np.pi)
    
    if use_PM:
        a /= eta**4
        b *= eta

    dvdf = 4*a*(b**2) * (3*np.tanh(b*field0)**2 / np.cosh(b*field0)**4 - 2*np.tanh(b*field0)**4 / np.cosh(b*field0)**2)
    return dvdf

    
def lf4V(field0,use_PM=True, **kwargs):
    default = {'l':10**-4}
    ks = default
    for k in default.keys():
        if k in kwargs.keys():
            ks[k] = kwargs[k]
    eta = np.sqrt(8*np.pi)

    potential = ks['l'] * field0**4 /4
    if use_PM:
        potential /= eta**4
    return potential

def lf4DV(field0,use_PM=True, **kwargs):
    default = {'l':10**-4}
    ks = default
    for k in default.keys():
        if k in kwargs.keys():
            ks[k] = kwargs[k]
    eta = np.sqrt(8*np.pi)

    dvdf = ks['l'] * field0**3 
    if use_PM:
        dvdf /= eta**3
    return dvdf
    
    
def slow_roll(v,dv):
    return (dv/v)**2 / (16 * np.pi)


def k_list(pw_data1,L=64):
    ki = 2 * np.pi / L
    #Note: We only want to create ks for the spectrum values and not for the a column
    #Save as np.array for extra functionality
    return np.array([ki * i for i in range(1,pw_data1.shape[1])]) 

def phi_k(pw_data1,L=64):
    pw1 = pw_data1.drop('a',axis=1)
    ks = k_list(pw_data1,L=L)
    pw1 /= ks**5
    return pw1

def phi_dot_k(pw_data2,L=64):
    pw2 = pw_data2.drop('a',axis=1)
    ks = k_list(pw_data2,L=L)
    pw2 /= ks**5
    return pw2
    
def k_list(pw_data1,L=64):
    ki = 2 * np.pi / L
    #Note: We only want to create ks for the spectrum values and not for the a column
    #Save as np.array for extra functionality
    return np.array([ki * i for i in range(1,pw_data1.shape[1])]) 
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

def omegas(pw1,data,L=64,skip=50,LH=0.8,h0=-1):
    data_size = data[::skip].shape[0]
    pw_size = pw1.shape[0]
    min_size,max_size = min(data_size,pw_size),max(data_size,pw_size)
    if (int(min_size*1.1) > max_size):
        data = data[:min_size*skip]
        pw1 = pw1[:min_size]
        
    else:
        raise ValueError("data and pw1/pw2 have two large a discrepancy in shapes: ",data.shape,pw1.shape)
    
    ks = k_list(pw1,L=L)
    a_list = pw1.a
    if h0==-1:
        h0 = data.h[0]
    ks *= (L*h0)/LH
    
    ks_a_matrix = np.einsum('a,b->ba',ks**2, 1/a_list**2)
    # Column containing the U''(chi) terms for the different values of 'a'.
    d2v_column = np.array(palatiniD2V(data[::skip]['mean1'],use_PM=False))
    
    # Assemble into an omega matrix (which will multiply the second term in the square brackets in Rubio, eq.4.20).
    omegas = ks_a_matrix + d2v_column.reshape(-1,1)
    # For information: print information about tachyonic modes. 
    print("Tachyonic modes present: ", (omegas<0).any())
    if (omegas<0).any(): print("Number of tachyonic modes: ", np.count_nonzero(omegas<0))
    return ks_a_matrix, d2v_column.reshape(-1,1)
    
def n_palatini(pw1,pw2,data, L=64, skip=2,h0=-1,LH=6,use_ks=True):
    """
    Generate the 'number density' object found in the integral of the perturbation energy density.
    For a derivation, see Rubio eq.4.20.

    Parameters
    ----------
    pw1 : pw_data1, pd.DataFrame
        Field modes data. Dataframe (n x (m+1): n rows for 'a' evolution; m+1 columns for 'k' evolution and one for 'a' values.)
    pw2 : pw_data2, pd.DataFrame
        Momentum modes data. Dataframe (n x (m+1): n rows for 'a' evolution; m+1 columns for 'k' evolution and one for 'a' values.)
    data : screen datafile, pd.DataFrame
        Contains all the generic time-evolution behaviours and variables.
    L : int, optional
        Resolution of simulation: number of points in the side of the lattice (i.e. given by 'n' in HLATTICE). The default is 64.
    skip : int, optional
        Interval between saves: How often the checkpoint files are updated. Required to know how to extract corresponding values from 'data'. The default is 2.
    h0 : float, optional
        Initial hubble value. If not specified (i.e. h0==-1), obtained from the first line of 'data'. The default is -1.
    LH : float, optional
        Initial_boxsize_times_h. Value from the HLATTICE simulation. The default is 6.
    use_ks : bool, optional
        Option to extract the mode values from the power spectra (default) or to use the raw values. The default is True.

    Returns
    -------
    df : pd.DataFrame
        Contains the full time-evolution of the perturbation for all modes.

    """    
    # Due to mismatches due to downloads during runtime of simulation: trim to correct size
    data_size = data[::skip].shape[0]
    pw_size = pw1.shape[0]
    min_size,max_size = min(data_size,pw_size),max(data_size,pw_size)
    if (int(min_size*1.1) > max_size):
        data = data[:min_size*skip]
        pw1 = pw1[:min_size]
        pw2 = pw2[:min_size]
        print("Trimmed data, pw1 and pw2: ",data.shape,pw1.shape,pw2.shape)
    else:
        raise ValueError("data and pw1/pw2 have two large a discrepancy in shapes: ",data.shape,pw1.shape)

    pw_a = pw1.drop('a',axis=1)
    pw_b = pw2.drop('a',axis=1)
    a_list = pw1['a']
    
    # Ensure that the ks have the correct dimensionality/scale.
    ks = k_list(pw1,L=L)
    if h0==-1:
        h0 = data.h[0]
    ks *= (L*h0)/LH
    
    if use_ks:
        phis = phi_k(pw1,L=L)
        phi_dots = phi_dot_k(pw2,L=L)
    else:
        phis= pw_a
        phi_dots = pw_b
    
    # Mode + time-evolution: corresponds to k^2/a^2 column repeated multiple times to form a matrix of size: (ks.size,a_list.size)
    # Note: may not be the only way/optimal way to achieve this, but maths (probably) check out.
    ks_a_matrix = np.einsum('a,b->ba',ks**2, 1/a_list**2)
    # Column containing the U''(chi) terms for the different values of 'a'.
    d2v_column = np.array(palatiniD2V(data[::skip]['mean1'],use_PM=False))
    
    
    # Assemble into an omega matrix (which will multiply the second term in the square brackets in Rubio, eq.4.20).
    omegas = ks_a_matrix + d2v_column.reshape(-1,1)
    # For information: print information about tachyonic modes. 
    print("Tachyonic modes present: ", (omegas<0).any())
    if (omegas<0).any(): print("Number of tachyonic modes: ", np.count_nonzero(omegas<0))
    
    # Add the terms together to form the number density.
    term1 = phi_dots
    term2 = phis * omegas
    df = 1/2 * (term1 + term2)
    # Add the information about ks, metric, and hubble to the dataframe (useful during integration calculations).
    df.columns = ks
    df['a'] = a_list
    df['h'] = data[::skip]['h'].to_numpy()
    return df

def integrate_perturbation(n_pal,Mpl=1024):
    """
    WORK IN PROGRESS: integrate over perturbations to get the energy density growth. 
    Currently not producing the expected results.    
    
    Parameters
    ----------
    n_pal : n x m Matrix
        m columns for each k mode. n rows for the evolution of the metric 'a'.
    Mpl : TYPE, optional
        Value of the Reduced Planck Mass in the HLATTICE simulation. The default is 1024.

    Returns
    -------
    perturbations : one-dimensional array
        Time-evolution of the integrated energy density.

    """
    # Configure the dataframe by first storing the columns 'a' and 'h' before removing them.
    a_list = n_pal['a']
    hubbles = n_pal['h']
    ns = n_pal.drop(['a','h'],axis=1)
    # Ensure the k-modes array has a (n, 1) shape, instead of (n,), which allows it to be included in calculations.
    ks = ns.columns.to_numpy().reshape(-1,1)

    # List of set size which will contain the results of the integration.
    integrated_column = [0] * ns.shape[0]
    for i in range(ns.shape[0]):
        # Evaluate the y values to be integrated over. Reshape statements ensure the objects can be used in calculations.
        ys = ns.iloc[i,:].to_numpy().reshape(-1,1) * ks**2
        # Store integral for a given metric value 'a' as a temp float (useful while troubleshooting).
        temp = integrate.simps( ys.reshape(-1),x=ks.reshape(-1))
        # Update the corresponding value in the list.
        integrated_column[i] = temp
    
    # Convert list to numpy array.
    integrated_arr = np.array(integrated_column)
    print("Shape: ",integrated_arr.shape)
    # Perform the remaining non-integral calculations on the perturbation energy densities.
    perturbations =  integrated_arr / (2*np.pi**2) / (3*hubbles**2 * Mpl**2)
    return perturbations
    




if __name__=="__main__":
    print('foo')
    conds = {'l':10**-4}
    p0 = init_momenta(palatiniV, palatiniDV, 5*10**-4, **conds)    
    print(p0)
    xs = np.linspace(0,10**-5,1000)
    import matplotlib.pyplot as plt
    v = palatiniV(xs)
    dv = palatiniDV(xs)
    #plt.plot(xs,lf4DV(xs))
    
    
    #plt.plot(xs,lf4V(xs)*10**5)
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.plot(xs, slow_roll(v, dv))


