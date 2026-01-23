###############################################################################
# Date Revised: 17 Sept 2025
# Author: sws02hs
# Description: Basic tests of advection schemes
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

def compareSchemeSnapshots(init, schemes, outFile,
                           nx = 40, nt = 10, u = 1, endTime = 0.4,
                           xLim = [0,1], yLim = [None,None]):
    '''
    Author: sws02hs
    Compares snapshots of results of advection schemes using periodic boundaries
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The minimim (xLim[0]) and maximim (xLim[1]) x-values of
                      the domain. xLim[1] is not part of the domain when periodic
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If [None,None] is specified then default limits
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
    
    # Set the initial conditions
    phi0 = init(x)
    plt.cla()
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    
    # Find and plot the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]))
    plt.plot(x, phiA, 'k', label = "Analytic", linewidth=2)
    
    # Calculate the results of all of the schemes and plot
    for i in range(len(schemes)):
        # The solution, phi, is stored as a list of arrays different
        # time steps. Initially just the initial conditions. The latest
        # time level array is always phi[-1]
        phi = [phi0.copy()]
        for n in range(nt):
            schemes[i](phi, c)
    
        plt.plot(x, phi[-1], cols[i], label = schemes[i].__name__)

    
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.xlabel('x')
    plt.ylabel(r'$\phi$')
    plt.legend()
    plt.savefig("plots/"+outFile)

def cosBell_0_05(x):
    "cosBell initial conditions with a=0, x=0.5"
    return cosBell(x, 0, 0.5)

def sqWave_0_05(x):
    "squareWave initial conditions with a=0, x=0.5"
    return squareWave(x, 0, 0.5)

def mass_test(x):
    "sum of square waves"
    return squareWave(x, 0, 0.5) + 0.5 * squareWave(x, 0.3, 0.6) + 0.75 * squareWave(x, 0.4, 0.8) + 0.25 * squareWave(x, 0.9, 1)

def tv_test(x):
    "sum of both square wave and cosine bell"
    return squareWave(x, 0, 0.5) + cosBell(x, 0.25, 0.75)

# List of functions to do the advection
schemes = [FTBS, CTCS, FTCS]
schemes_new = [semi_lagrangian, tvd, implicit]

#compareSchemeSnapshots(cosBell_0_05, schemes_new, "newSchemesCosBell.jpg", nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])
#compareSchemeSnapshots(sqWave_0_05, schemes_new, "newSchemesSqWave.jpg", nx = 40, nt = 20, u = 1, endTime = 0.2, xLim = [0,1])
#compareSchemeSnapshots(mass_test, schemes_new, "mass_test_nx_80_nt_167_profiles.jpg", nx = 80, nt = 167, u = 1, endTime = 1)
#compareSchemeSnapshots(sqWave_0_05, schemes_new, "tv_nx_80_nt_160_profiles.jpg", nx = 80, nt = 160, u = 1, endTime = 1)