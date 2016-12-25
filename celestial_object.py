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

######################################
# Define class for celestial bodies. #
######################################

# This works at the moment only for elliptical (generic) orbits. Fix this!

class celestial_body:
    # This class assumes a reference coordinate system such that a large mass is situated at the origin. It might actually assume some more things.
    
    ####### Init #######
    
    def __init__(self,mass,mu,semi_major_axis,eccentricity,inclination,longitude_ascending_node,argument_periapsis,true_anomaly_epoch):
        # Initialization of class using classical orbital elements a, e, i, Omega, omega, nu_0
        self.semi_major_axis = semi_major_axis # a
        self.energy = - mu / ( 2.0 * self.semi_major_axis ) # E
        self.eccentricity = eccentricity # e
        if self.energy < 0:
            if self.eccentricity == 0:
                self.type = "circular"
            else:
                self.type = "elliptical"
        elif self.energy == 0:
            self.type = "parabolic"
        else:
            self.type = "hyperbolic"
        self.inclination = inclination # i
        if inclination == 0:
            self.planar == True
        else:
            self.planar == False
        
        if self.planar == False:
            self.longitude_ascending_node = longitude_ascending_node # Omega
            self.argument_periapsis = argument_periapsis # omega
        else:
            self.longitude_ascending_node = 0
            self.argument_periapsis = 0
        self.true_anomaly_epoch = true_anomaly_epoch # nu
        self.mass = mass # m
        self.parameter = semi_major_axis * (1 - eccentricity**2) # p
        if ( 0 <= self.true_anomaly_epoch ) and ( self.true_anomaly_epoch <= np.pi):
            self.eccentric_anomaly = np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E, at the moment the cases dont't cover everything.
        else:
            self.eccentric_anomaly = 2 * np.pi - np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E
        self.mean_anomaly = self.eccentric_anomaly - self.eccentricity * np.sin(self.eccentric_anomaly) # M
        self.mean_motion = np.sqrt(mu / self.semi_major_axis**3 ) # n
        self.period = 2 * np.pi / np.sqrt(mu) * np.sqrt(self.semi_major_axis**3) # T
        self.mu = mu # mu
        self.X = 0 # X for universal formulation of time of flight
    
    @classmethod
    def from_position_velocity(self,mass,mu,position,velocity):
        # Initialization of class using position and momentum
        # For this purpose we need to calculate various intermediate objects. Should we save them for later? Is it more clever to just use position and momentum all the time?
        
        h = np.cross(position,velocity) # Calculate angular momentum h
        if h != [0,0,0]:
            n = np.cross(np.array([0,0,1],float),h) # Calculate node vector
            e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
            p = np.dot(h,h) / mu
            
            # Is it better to just save the cosine of the angles?
            semi_major_axis = p / (1-np.dot(e,e))
            eccentricity = norm(e)
            inclination = np.arccos(h[2] / norm(h))
            if position[1] >= 0:
                longitude_ascending_node = np.arccos(n[0] / norm(n))
            else:
                longitude_ascending_node = 2 * np.pi - np.arccos(n[0] / norm(n))
            if e[2] >= 0:
                argument_periapsis = np.arccos(np.dot(n,e) / (norm(n) * norm(e)))
            else:
                argument_periapsis = 2 * np.pi - np.arccos(np.dot(n,e) / (norm(n) * norm(e)))
            if np.dot(position,velocity) >= 0:
                true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
            else:
                true_anomaly_epoch = 2 * np.pi - np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
            
            body = celestial_body(mass,mu,semi_major_axis,eccentricity,inclination,longitude_ascending_node,argument_periapsis,true_anomaly_epoch)
            return body
        else:
            return celestial_object.initialize_collision_orbit(mass,mu,position,velocity)
    
    @classmethod
    def initialize_collision_orbit(self,mass,mu,position,velocity):
        pass
    
    ####### Export #######
    
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?
        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))

        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.
        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [ np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
        ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points):
        # Returns a list of three dimensional coordinates for the orbit.
        position = np.zeros( (number_points,3) )
        interval = 2 * np.pi / number_points
        for i in range(number_points):
            position[i,:] = self.calculate_advance_in_true_anomaly(i * interval)[0]
        return np.vstack( (position,position[0,:]) )
        
    ###### Advance along orbit #######
    
    def advance_in_time(self,delta_t):
        # This method advances the object on its course by delta t in time. This means that it needs to translate the time difference into changes in the true anomaly at epoch and then add this number to the existing value.
        # delta_t should be small enough such that the body does not evolve more than one period. Is this necessary?
        
        # Update mean anomaly. Ignore full rotations.
        new_mean_anomaly = self.mean_motion * delta_t + self.mean_anomaly
        
        # Solve E-e*sin(E)=M numerically
        new_eccentric_anomaly = fsolve(lambda E : E - self.eccentricity * np.sin(E) -new_mean_anomaly,new_mean_anomaly) 
        
        # Calculate new true anomaly at epoch
        if new_eccentric_anomaly <= np.pi:
            new_true_anomaly_epoch = np.arccos( ( np.cos(new_eccentric_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cos(new_eccentric_anomaly)))
        else:
            new_true_anomaly_epoch = 2 * np.pi - np.arccos( ( np.cos(new_eccentric_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cos(new_eccentric_anomaly)))
        
        # Update values of true anomaly at epoch and eccentric anomaly and mean anomaly
        self.true_anomaly_epoch = new_true_anomaly_epoch
        self.mean_anomaly = new_mean_anomaly
        self.eccentric_anomaly = new_eccentric_anomaly
    
    def t_in_dep_of_X(self, X):
        r_0, v_0 = self.export_postion_velocity()
        return 1 / np.sqrt(self.mu) * ( np.dot(r_0,v_0) /np.sqrt(self.mu) * X**2 * C(X) + ( 1 - norm(r_0) / self.semi_major_axis ) * X**3 * S(X) + norm(r_0) * X )
        
    def advance_in_time_universal(self,delta_t):
        # This method advances the object on its course by delta t in time using the universal time of fligt formulation. This means it should be usable for all kinds of orbits.
        
        # Solve for new X
        new_X = fsolve(lambda X : self.t_in_dep_of_X(X) - delta_t,delta_t)
    
    def advance_in_true_anomaly(self,delta_nu):
        # This method increases the true anomaly by a given input. It can be used to find equi-distant-angle points on the orbit for visualization purposes. It also updates eccentric anomaly and mean anomaly.
        self.true_anomaly_epoch = self.true_anomaly_epoch + delta_nu
        if self.true_anomaly_epoch <= np.pi:
            self.eccentric_anomaly = np.arccos( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))
        else:
            self.eccentric_anomaly = 2 * np.pi - np.arccos( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))
        self.mean_anomaly = self.eccentric_anomaly - self.eccentricity * np.sin( self.eccentric_anomaly )
    
    def calculate_advance_in_true_anomaly(self,delta_nu):
        # This method advances the object on its course by delta nu in true anomaly and returns the new position. It is useful for calculating points on the orbit without actually advancing the object itself.
        new_true_anomaly_epoch = self.true_anomaly_epoch + delta_nu
        
        r = self.parameter  / ( 1 + self.eccentricity * np.cos(new_true_anomaly_epoch))

        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.
        position_perifocal_system = np.array([r * np.cos(new_true_anomaly_epoch),r * np.sin(new_true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(new_true_anomaly_epoch),self.eccentricity + np.cos(new_true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [ np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
        ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
        
    
        











