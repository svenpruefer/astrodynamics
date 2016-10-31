##############################
# Import necessary libraries #
##############################

import numpy as np

#####################
# Define parameters #
#####################

mu = 1.0 # Gravitational Parameter, equal to product of gravitational constant and central mass in restricted 2-body approximation

##################################
# Define various math functions. #
##################################
# Are they implemented somewhere else?

# Cross product is implemented as np.cross(a,b)
def cross_product(v,w):
    product = np.array([v[1]*w[2]-w[1]*v[2],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]],float)
    return product

def norm(v):
    return np.sqrt(np.dot(v,v))

######################################
# Define class for celestial bodies. #
######################################

class celestial_body:
    # This class assumes a reference coordinate system such that a large mass is situated at the origin. It might actually assume some more things.
    
    def __init__(self,mass,semi_major_axis,eccentricity,inclination,longitude_ascending_node,argument_periapsis,true_anomaly_epoch):
        # Initialization of class using classical orbital elements a, e, i, Omega, omega, nu_0
        self.semi_major_axis = semi_major_axis
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.longitude_ascending_node = longitude_ascending_node
        self.argument_periapsis = argument_periapsis
        self.true_anomaly_epoch = true_anomaly_epoch
        self.mass = mass
    
    @classmethod
    def from_position_velocity(self,mass,position,velocity):
        # Initialization of class using position and momentum
        # For this purpose we need to calculate various intermediate objects. Should we save them for later? Is it more clever to just use position and momentum all the time?
        
        h = np.cross(position,velocity) # Calculate angular momentum h
        n = np.cross(np.array([0,0,1],float),h) # Calculate node vector
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu
        
        # Is it better to just save the cosine of the angles?
        semi_major_axis = p / (1-np.dot(e,e))
        eccentricity = norm(e)
        inclination = np.arccos(h[2] / norm(h))
        longitude_ascending_node = np.arccos(n[0] / norm(n))
        argument_periapsis = np.arccos(np.dot(n,e) / (norm(n) * norm(e)))
        true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
        
        body = celestial_body(mass,semi_major_axis,eccentricity,inclination,longitude_ascending_node,argument_periapsis,true_anomaly_epoch)
        return body

position = np.array([3.0 / 4 * np.sqrt(3), 3.0 /4, 0],float)
velocity = np.array([-1.0/(2*np.sqrt(2)), np.sqrt(3) / (2 * np.sqrt(2)), 1.0 / np.sqrt(2)],float) 
test_body = celestial_body.from_position_velocity(1,position,velocity)
print(test_body.inclination)













