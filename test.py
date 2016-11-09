##########################################
# Import necessary classes and libraries #
##########################################

from celestial_object import *
from kepler import *
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

############################
# Define Figure parameters #
############################

fig = plt.figure()
ax = fig.gca(projection='3d')

#####################
# Define Parameters #
#####################

mu = 1.0 # Gravitational Parameter, equal to product of gravitational constant and central mass in restricted 2-body approximation

###################
# Define Objects #
###################

position_1 = np.array([3.0 / 4 * np.sqrt(3), 3.0 /4, -0.1],float)
velocity_1 = np.array([-1.0/(2*np.sqrt(2)), np.sqrt(3) / (1.6 * np.sqrt(2)), 1.0 / np.sqrt(2)],float)
mass_1 = 2
#test_body_1 = celestial_body.from_position_velocity(1,mu,position_1,velocity_1)

position_2 = np.array([2.5 / 4 * np.sqrt(3), 2.9 /4, 0.3],float)
velocity_2 = np.array([0.7/(2*np.sqrt(2)), -np.sqrt(3) / (2.4 * np.sqrt(2)), -1.0 / np.sqrt(2)],float)
mass_2 = 1

mu = 1 * (mass_1 + mass_2)
#test_body_2 = celestial_body.from_position_velocity(1,mu,position_2,velocity_2)

#body_1, body_2 = kepler_problem(1,position_1,velocity_1,1,position_2,velocity_2)

total_body = kepler_problem(mass_1,position_1,velocity_1,mass_2,position_2,velocity_2)

#############
# Plotting #
#############


#position1 = body_1.export_orbit(100)
#position2 = body_2.export_orbit(100)

position1 = mu / mass_1 * total_body.export_orbit(100)
position2 = - mu / mass_2 * total_body.export_orbit(100)

#lines = np.vstack( (position1[20,:],position2[20,:]) )

#print(position1)
#print(lines)

ax.scatter(np.array([0]),np.array([0]),np.array([0]))
ax.plot(position1[:,0],position1[:,1],position1[:,2])
ax.plot(position2[:,0],position2[:,1],position2[:,2])
#ax.plot(lines[:,0],lines[:,1],lines[:,2])

ax.set_zlim(-5.0,5.0)
plt.show()
