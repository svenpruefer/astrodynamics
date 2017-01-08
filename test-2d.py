##########################################
# Import necessary classes and libraries #
##########################################

from kepler import *
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


############################
# Define Figure parameters #
############################

fig = plt.figure()
#ax = fig.gca#(projection='3d')

#####################
# Define Parameters #
#####################

mu = 1.0 # Gravitational Parameter, equal to product of gravitational constant and central mass in restricted 2-body approximation

###################
# Define Objects #
###################

position_1 = np.array([3.0 / 4 * np.sqrt(3), 3.0 /4, 0],float)
velocity_1 = np.array([-1.0/(2*np.sqrt(2)), np.sqrt(3) / (1.6 * np.sqrt(2)), 0],float)
test_body_1 = celestial_body.create_object(position_1,velocity_1)

position_2 = np.array([0.4, 1, 0],float)
velocity_2 = np.array([1.2, 0.7, 0],float)
test_body_2 = celestial_body.create_object(position_2,velocity_2)

#############
# Plotting #
#############

def plot(body):
    orbit = body.export_orbit(100)
    pos,vec = body.export_position_velocity()

    # Sun
    plt.scatter(np.array([0]),np.array([0]),color='yellow')

    # Celestial body
    plt.scatter([pos[0]],[pos[1]],color='red')

    # Velocity vector
    plt.quiver(pos[0],pos[1],vec[0],vec[1],color='red',scale=1,units='xy',scale_units='xy',width=0.02)

    #Orbit
    plt.plot(orbit[:,0],orbit[:,1])
    
plot(test_body_2)
test_body_2.advance_in_true_anomaly(0.3)
plot(test_body_2)

plt.xlim(-3,3)
plt.ylim(-3,3)
plt.show()
