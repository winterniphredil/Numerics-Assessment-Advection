# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 09:34:46 2026

@author: hanna
"""

import numpy as np
import math

def semi_lagrangian(phi, c):
    '''
    Author: Hannah
    Solves the linear advection equation for one time step using semi-Lagrangian advection with cubic interpolation
    phi is updated. Nothing is returned
    Parameters:
        phi ([1darray]): A list length one containing a 1darray which is the
                         latest value of phi, to be updated.
        c (float): The Courant number
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
        k = math.floor(i-c)
        beta = i - k - c
        phiNew[i] = -1/6. * beta * (1-beta) * (2-beta) * phiOld[(k-1) % nx] + 1/2. * (1+beta) * (1-beta) * (2-beta) * phiOld[k % nx] + 1/2. * (1+beta) * beta * (2-beta) * phiOld[(k+1) % nx]  - 1/6. * (1+beta) * beta * (1-beta) * phiOld[(k+2) % nx]
        
    
    # Remove the old value of phi
    phi.pop(0)
    return None


def tvd(phi, c):
    '''
    Author: Hannah
    Solves the linear advection equation for one time step using a TVD scheme wth Van Leer limiter and FV method
    phi is updated. Nothing is returned
    Parameters:
        phi ([1darray]): A list length one containing a 1darray which is the
                         latest value of phi, to be updated.
        c (float): The Courant number
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
    phi_half = [0]*nx  # phi_half[i] corresponds to phi[i+1/2]
    
    for i in range(nx):
        denom = (phiOld[(i+1) % nx] - phiOld[i])
        num = (phiOld[i] - phiOld[(i-1) % nx])
        r = ( num / denom if denom != 0 else (1 if num == 0 else 0))
        phi_H = 1/2. * (1+c) * phiOld[i] + 1/2. * (1-c) * phiOld[(i+1) % nx]
        phi_L = (phiOld[i] if c>=0 else phiOld[(i+1) % nx])
        limiter = (r + abs(r)) / (1+ abs(r))
        
        phi_half[i] = limiter * phi_H + (1-limiter) * phi_L
        
    for i in range(nx):
        phiNew[i] = phiOld[i] - abs(c) * (phi_half[i] - phi_half[(i-1) % nx])
    
    # Remove the old value of phi
    phi.pop(0)
    return None


def A_matrix(nx, c):
    A = np.diag(np.array([1+c*c/2.]*nx))
    for i in range(nx):
        A[i,(i-1)%(nx)] = - c/4. * (1+c)
        A[i,(i+1)%(nx)] = c/4. * (1-c)
    return A

def B_matrix(nx, c):
    B = np.diag(np.array([1-c*c/2.]*nx))
    for i in range(nx):
        B[i,(i-1)%(nx)] = c/4. * (c+1)
        B[i,(i+1)%(nx)] = c/4. * (c-1)
    return B


def lax_wendroff_cn(phi, c):
    '''
    Author: Hannah
    Solves the linear advection equation for one time step using a Lax-Wendorff scheme with Crank-Nicolson in time 
    phi is updated. Nothing is returned
    Parameters:
        phi ([1darray]): A list length one containing a 1darray which is the
                         latest value of phi, to be updated.
        c (float): The Courant number
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
    
    A = A_matrix(nx, c)
    B = B_matrix(nx, c)
    
    phiNew = np.linalg.solve(A, B @ phiOld)
    
    phi[1] = phiNew
    # Remove the old value of phi
    phi.pop(0)
    return None