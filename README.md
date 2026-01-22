MTMFMD advection assignment. Original code from Hilary Weller <h.weller@reading.ac.uk>

Code to solve the advection equation.

advection.py          : snapshots of results for each scheme
                        
growingErrors.py      : l2 errors, total variation and mass conservation 

advectionHelpers.py   : functions for initialising and calculating properties

schemes_sws02hs.py    : FTCS, CTCS and FTBS advection schemes
schemes.py            : semi-Lagrangian, TVD, and implicit schemes
