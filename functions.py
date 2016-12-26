##############################
# Import necessary libraries #
##############################

import numpy as np
from scipy.optimize import fsolve

##################################
# Define various math functions. #
##################################

def norm(v):
    return np.sqrt(np.dot(v,v))
    
def S(z):
    return ( np.sqrt(z) - np.sin(np.sqrt(z)) ) / np.sqrt(z**3)

def C(z):
    return ( 1 - np.cos(np.sqrt(z)) ) / z
