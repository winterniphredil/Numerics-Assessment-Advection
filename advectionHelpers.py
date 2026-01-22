# Helper functions for testing advection schemes.
# Authors: sws02hs

import numpy as np

def cosBell(x, a, b):
    '''
    Returns a cosine bell shape between parameters a and b
    Author: sws02hs
    Parameters:
        x (1darray): The x-locations for calculating the cosine bell.
        a (float): The start of the non-zero bell
        b (float): The end of the non-zero bell
    Returns:
        phi : 1darray. A cosine bell shape between x=a and x=b
    '''
    return np.where((x > a) & (x < b), 0.5*(1 - np.cos(2*np.pi*(x - a)/(b-a))), 0.)

def squareWave(x, a, b):
    '''
    Returns a square wave that is non-zero between parameters a and b
    Author: Hannah
    Parameters:
        x (1darray): The x-locations for calculating the square wave.
        a (float): The start of the wave
        b (float): The end of the wave
    Returns:
        phi : 1darray. A square wave between x=a and x=b
    '''
    return np.where((x > a) & (x < b), 1., 0.)

def l2error(phi, phiA):
    '''
    Returns the l2 error norm of phi in comparison to exact solution phiA.
    Author: sws02hs
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution
    Returns
        l2 (float): The l2 error norm
    '''
    phiError = phi - phiA
    l2 = np.sqrt(sum(phiError**2)/sum(phiA**2))
    return l2

def total_mass(phi):
    '''
    Returns the total mass.
    Author: Hannah
    Parameters:
        phi (1darray): The numerical solution
    Returns
        total mass
    '''
    return sum(phi)

def total_variation(phi):
    '''
    Returns the total variation of phi.
    Author: Hannah
    Parameters:
        phi (1darray): The numerical solution
    Returns
        tv: measurement of strength of oscillations
    '''
    nx = len(phi)
    tv = sum(abs(phi[(i+1) % nx] - phi[i]) for i in range(nx))
    return tv