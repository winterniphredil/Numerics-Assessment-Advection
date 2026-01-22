# Numerical schemes for solving the advection equation
# Authors: sws02hs

import numpy as np

def FTBS(phi, c):
    '''
    Author: sws02hs
    Solves the linear advection equation for one time step using FTBS.
    phi is updated. Nothing is returned
    Parameters:
        phi ([1darray]): A list length one containing a 1darray which is the
                         latest value of phi, to be updated.
        c (float): The Cournat number
    Returns:
        None
    '''
    # add another time level to the list, copying from the most recent
    phi.append(phi[-1].copy())
    # Labels for phi at time levels n and n+1
    phiOld = phi[0]
    phiNew = phi[1]
    
    # Update phi for all points
    nx = len(phiOld)
    for i in range(nx):
        phiNew[i] = phiOld[i] - c*(phiOld[i] - phiOld[(i-1)%nx])
    
    # Remove the old value of phi
    phi.pop(0)
    return None


def FTCS(phi, c):
    '''
    Author: sws02hs
    Solves the linear advection equation for one time step using FTCS.
    phi is updated. Nothing is returned
    Parameters:
        phi ([1darray]): A list length one containing a 1darray which is the
                         latest value of phi, to be updated.
        c (float): The Cournat number
    Returns:
        None
    '''
    # add another time level to the list, copying from the most recent
    phi.append(phi[-1].copy())
    # Labels for phi at time levels n and n+1
    phiOld = phi[0]
    phiNew = phi[1]
    
    # Update phi for all points
    nx = len(phiOld)
    for i in range(nx):
        phiNew[i] = phiOld[i] - 0.5*c*(phiOld[(i+1)%nx] - phiOld[(i-1)%nx])
    
    # Remove the old value of phi
    phi.pop(0)
    return None


def CTCS(phi, c):
    '''
    Author: sws02hs
    Solves the linear advection equation for one time step using CTCS.
    phi is updated. Nothing is returned.
    Parameters:
        phi ([1darray,1darray]): A list of two consecutive values of phi. If
                    only one is present, the next is created using FTCS.
        c (float): The Cournat number
    Returns:
        None
    '''
    # If phi[1] is not set, calculate it using FTCS
    if len(phi) == 1:
        phiOld = phi[0].copy()
        FTCS(phi, c)
        phi.insert(0, phiOld)
        return None
    
    # add another time level, copying from the oldest
    phi.append(phi[0].copy())

    # Labels for phi at time levels n-1, n and n+1
    phiOld = phi[0]
    phiMid = phi[1]
    phiNew = phi[2]
    
    # Update phiNew using and phiOld and phiMid
    nx = len(phiOld)
    for i in range(nx):
        phiNew[i] = phiOld[i] - c*(phiMid[(i+1)%nx] - phiMid[(i-1)%nx])
    
    # Remove the old value of phi
    phi.pop(0)
    return None