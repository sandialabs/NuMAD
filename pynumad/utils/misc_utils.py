########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import numpy as np

# SED-like substitution
def copy_and_replace(fin, fout, replacements):
    inf = open(fin, 'r')
    outf = open(fout, 'w')
    for line in inf:
        for src, target in replacements.items():
            line = line.replace(src, target)
        outf.write(line)
    inf.close()
    outf.close()

def LARCetaT(alp0):
    #TODO complete docstring
    """Compute the coefficient of transverse influence required for Larc failure criteria.

    "In the absence of biaxial test data, ?L can be estimated from the longitudinal
    and transverse shear strengths." Failure Criteria for FRP Laminates in Plane Stress
    Carlos G. Dávila, Pedro P. Camanho

    Parameters
    ----------
    alp0
        Material fracture angle, degrees
    
    Returns
    -------
    etaT 
    """
    num = -1
    denom = np.tan(np.deg2rad(2*alp0))
    with np.errstate(divide='ignore', invalid='ignore'):
        etaT = num / denom
    return etaT



def LARCetaL(SL,YC,alp0):
    #TODO complete docstring
    """Compute the coefficient of longitudinal influence required for Larc failure criteria.

    "In the absence of biaxial test data, ?L can be estimated from the longitudinal
    and transverse shear strengths." Failure Criteria for FRP Laminates in Plane Stress
    Carlos G. Dávila, Pedro P. Camanho

    Parameters
    ----------
    SL
        Lateral shear strength 
    YC
        Transverse compressive strength 
    alp0
        Material fracture angle, degrees
    Returns
    -------
    etaL
    """
    if alp0:
        num = -SL * np.cos(np.deg2rad(2*alp0))
        denom = (YC*np.cos(np.deg2rad(alp0))**2)
        with np.errstate(divide='ignore', invalid='ignore'):
            etaL = num / denom
    else:
        etaL = None
    return etaL

def _parse_data(data):
    """Helper function for parsing data from blade yaml files.
    
    Parameters
    ----------
    data 
        a number or list of numbers where numbers can be floats or strings
        e.g. 3.0 or '1.2e2'
    
    Returns
    -------
    parsed_data
        single float or array of floats
    """
    try:
        # detect whether data is list
        data+[]
    except TypeError: # case for single data point
        parsed_data = float(data)
    else:
        parsed_data = np.array([float(val) for val in data]) # case for list of data points
    finally:
        return parsed_data