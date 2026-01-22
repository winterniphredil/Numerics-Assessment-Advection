###############################################################################
# Date Revised: 19 Sept 2025
# Author: sws02hs
# Description: Create a graph of errors growing over time for advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import *
from schemes_sws02hs import *
from schemes import *

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

   

def plot_mass_error(init, schemes, outFile,
                  nx = 40, nt = 50, u = 1, endTime = 0.4, 
                  xLim = [0,1], yLim = None):
    '''
    Author: sws02hs, Hannah
    Creates a graph of mass difference over time for various advection schemes
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The maximum and minimun spacial dimensions
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    '''
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions for all schemes
    phi0 = init(x)
    
    # Store an array of all times
    times = np.linspace(0, endTime, nt+1)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    
    # Loop over all the schemes and plot
    for i in range(len(schemes)):
        # the mass for each time step
        mass_diff = np.zeros(nt+1)
        mass_0 = total_mass(phi0)
        # The solution for this scheme
        phi = [phi0.copy()]
        # Loop through all time steps
        for n in range(nt):
            # Update phi
            schemes[i](phi, c)
            # calculate mass
            mass_diff[n+1] = abs(total_mass(phi[-1]) - mass_0)

        # Plot the error for this scheme
        plt.semilogy(times, mass_diff, cols[i], label = schemes[i].__name__)
    
    plt.xlim([0, endTime])
    plt.ylim(yLim)
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel("mass difference")
    plt.savefig("plots/"+outFile)


def plot_tvd(init, schemes, outFile,
                  nx = 40, nt = 50, u = 1, endTime = 0.4, 
                  xLim = [0,1], yLim = None):
    '''
    Author: sws02hs, Hannah
    Creates a graph of TVD over time for various advection schemes
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The maximum and minimum spacial dimensions
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    '''
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions for all schemes
    phi0 = init(x)
    
    # Store an array of all times
    times = np.linspace(0, endTime, nt+1)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    
    # Loop over all the schemes and plot
    for i in range(len(schemes)):
        # the TV for each time step
        tv = np.zeros(nt+1)
        tv[0] = total_variation(phi0)
        # The solution for this scheme
        phi = [phi0.copy()]
        # Loop through all time steps
        for n in range(nt):
            # Update phi
            schemes[i](phi, c)
            # calcuate tv
            tv[n+1] = np.array(total_variation(phi[-1]))

        # Plot the error for this scheme
        plt.semilogy(times, tv, cols[i], label = schemes[i].__name__)
    
    plt.xlim([0, endTime])
    plt.ylim(yLim)
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel("Total variation")
    plt.savefig("plots/"+outFile)
    

def plot_inc_c(init, schemes, outFile,
                  start_nx = 40, nt = 50, u = 1, endTime = 0.4, 
                  xLim = [0,1], yLim = None):
    '''
    Author: sws02hs
    Creates a graph of error growth over time for various advection schemes
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The maximum and minimun spacial dimensions
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    '''   
    
    # Derived parameters
    dx = (xLim[1] - xLim[0])/start_nx
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    
    # Set the initial conditions for all schemes
    phi0 = init(x)
    
    # Store an array of all times
    times = np.linspace(0, endTime, nt+1)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    
    # Loop over all the schemes and plot
    for i in range(len(schemes)):
        # clear plot to start new scheme
        plt.cla()
        for j in range(5):
            # parameters
            nx = start_nx + j * (nt/4)
            dx = (xLim[1] - xLim[0])/nx
            x = np.arange(xLim[0], xLim[1], dx)
            dt = endTime/nt
            phi0 = init(x)
            c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
            print('Courant number =', c)
            # The l2 error for each time step
            l2errors = np.zeros(nt+1)
            # The solution for this scheme
            phi = [phi0.copy()]
            # Loop through all time steps
            for n in range(nt):
                # Update phi
                schemes[i](phi, c)
                # Calculate the analytic solution
                phiA = init((x - u*times[n+1])%(xLim[1] - xLim[0]))
                # Calculate the l2 error
                l2errors[n+1] = l2error(phi[-1], phiA)

            # Plot the error for this scheme
            plt.semilogy(times, l2errors, cols[j], label = " c = "+str(c))
        
        plt.xlim([0, endTime])
        plt.ylim(yLim)
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel(r'$\ell_2$ error')
        plt.savefig("plots/"+schemes[i].__name__+outFile)
        plt.pause(0.05)
            
           


def cosBell_0_05(x):
    "cosBell initial conditions with a=0, b=0.5"
    return cosBell(x, 0, 0.5)

def sqWave_0_05(x):
    "squareWave initial conditions with a=0, b=0.5"
    return squareWave(x, 0, 0.5)

def mass_test(x):
    "sum of square waves"
    return squareWave(x, 0, 0.5) + 0.5 * squareWave(x, 0.3, 0.6) + 0.75 * squareWave(x, 0.4, 0.8) + 0.25 * squareWave(x, 0.9, 1)

def tv_test(x):
    "sum of both square wave and cosine bell"
    return squareWave(x, 0, 0.5) + cosBell(x, 0.25, 0.75)
    
schemes = [FTBS, CTCS, FTCS]

schemes_new = [semi_lagrangian, tvd, lax_wendroff_cn]

#growingErrors(cosBell_0_05, schemes_new, "basicSchemesCosMass.jpg", nx = 80, nt = 100, u = 1, endTime = 1)

#plot_tvd(tv_test, schemes_new, "tv_nx_80_nt_160.jpg", nx = 80, nt = 160, u = 1, endTime = 1)

plot_inc_c(cosBell_0_05, schemes_new, "_courant.jpg", start_nx = 60, nt = 100, endTime = 1)