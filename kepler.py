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
        self.parameter = semi_major_axis * (1 - eccentricity**2)
    
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
    
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?
        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))

        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.
        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [ np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
        ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity

#############
# Test Code #
#############

position = np.array([3.0 / 4 * np.sqrt(3), 3.0 /4, 0],float)
velocity = np.array([-1.0/(2*np.sqrt(2)), np.sqrt(3) / (2 * np.sqrt(2)), 1.0 / np.sqrt(2)],float) 
test_body = celestial_body.from_position_velocity(1,position,velocity)
print(test_body.export_position_velocity())













