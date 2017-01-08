
##############################
# Import necessary libraries #
##############################

from functions import *
from celestial_body_base import *
import numpy as np
from scipy.optimize import fsolve

#################################################
# Define class for non-planar celestial bodies. #
#################################################

class celestial_body_on_nonplanar_parabolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "parabolic"
        self.planar = False

        h = np.cross(position,velocity) # Calculate angular momentum h
        n = np.cross(np.array([0,0,1],float),h) # Calculate node n
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu
        self.parameter = p
        
        self.eccentricity = norm(e)
        self.inclination = np.arccos(h[2] / norm(h))
        if position[1] >= 0:
            self.longitude_ascending_node = np.arccos(n[0] / norm(n))
        else:
            self.longitude_ascending_node = 2 * np.pi - np.arccos(n[0] / norm(n))
        if e[2] >= 0:
            self.argument_periapsis = np.arccos(np.dot(n,e) / (norm(n) * norm(e)))
        else:
            self.argument_periapsis = 2 * np.pi - np.arccos(np.dot(n,e) / (norm(n) * norm(e)))
        if np.dot(position,velocity) >= 0:
            self.true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
        else:
            self.true_anomaly_epoch = 2 * np.pi - np.arccos(np.dot(e,position) / (norm(e) * norm(position)))

        self.energy = 0 # E

        # Parabolic Anomaly
        self.parabolic_anomaly = np.sqrt(self.parameter) * np.tan(self.true_anomaly / 2.0) # D
                
        # Mean anomaly
        self.mean_anomaly = self.parameter * self.parabolic_anomaly + 1.0 / 3 * self.parabolic_anomaly**3
        self.mean_motion = 2 * np.sqrt(mu) # n

        #Universal time of flight
        self.X = 0 # X

    ######################
    ####### Export #######
    ######################
        
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?
        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))
            
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.
        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
            ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points=60):
        # Returns a list of three dimensional coordinates for the orbit. Remove point opposite perigree.
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
        new_parabolic_anomaly = fsolve(lambda D : self.parameter * D + 1.0 / 3 * D**3 - new_mean_anomaly, new_mean_anomaly) 
        
        # Calculate new true anomaly at epoch
        new_true_anomaly_epoch = 2 * np.arctan( new_parabolic_anomaly / np.sqrt(self.parameter) )
            
        # Update values of true anomaly at epoch and eccentric anomaly and mean anomaly
        self.true_anomaly_epoch = new_true_anomaly_epoch
        self.mean_anomaly = new_mean_anomaly
        self.parabolic_anomaly = new_parabolic_anomaly
        
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
        self.parabolic_anomaly = np.sqrt(self.parameter) * np.tan( self.true_anomaly / 2.0 )
        self.mean_anomaly = self.parameter * self.parabolic_anomaly + 1.0 / 3 * self.parabolic_anomaly**3
        
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
   

######################################################
# Define class for planar elliptic celestial bodies. #
#####################################################

class celestial_body_on_planar_parabolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "parabolic"
        self.planar = True

        h = np.cross(position,velocity) # Calculate angular momentum h
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu
        self.parameter = p
        self.eccentricity = norm(e)
        
        self.argument_periapsis = np.arccos(np.dot(np.array([1,0,0],float),e) / (norm(e)))
        
        if np.dot(position,velocity) >= 0:
            self.true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
        else:
            self.true_anomaly_epoch = 2 * np.pi - np.arccos(np.dot(e,position) / (norm(e) * norm(position)))

        self.energy = 0 # E

        # Parabolic Anomaly
        self.parabolic_anomaly = np.sqrt(self.parameter) * np.tan(self.true_anomaly / 2.0) # D
                
        # Mean anomaly
        self.mean_anomaly = self.parameter * self.parabolic_anomaly + 1.0 / 3 * self.parabolic_anomaly**3
        self.mean_motion = 2 * np.sqrt(mu) # n

        #Universal time of flight
        self.X = 0 # X
    
        ####### Export #######
        
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?

        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))
            
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.

        rotation_matrix = np.array([[np.cos(self.argument_periapsis), - np.sin(self.argument_periapsis), 0],\
        [np.sin(self.argument_periapsis), - np.cos(self.argument_periapsis), 0],\
        [0, 0, 1]],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points=60):
        # Returns a list of three dimensional coordinates for the orbit. Remove point opposite of perigree.
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
        new_parabolic_anomaly = fsolve(lambda D : self.parameter * D + 1.0 / 3 * D**3 - new_mean_anomaly, new_mean_anomaly) 
        
        # Calculate new true anomaly at epoch
        new_true_anomaly_epoch = 2 * np.arctan( new_parabolic_anomaly / np.sqrt(self.parameter) )
            
        # Update values of true anomaly at epoch and eccentric anomaly and mean anomaly
        self.true_anomaly_epoch = new_true_anomaly_epoch
        self.mean_anomaly = new_mean_anomaly
        self.parabolic_anomaly = new_parabolic_anomaly
        
    def t_in_dep_of_X(self, X):
        r_0, v_0 = self.export_postion_velocity()
        return 1 / np.sqrt(self.mu) * ( np.dot(r_0,v_0) / np.sqrt(self.mu) * X**2 * C(X) + ( 1 - norm(r_0) / self.semi_major_axis ) * X**3 * S(X) + norm(r_0) * X )
    
    def advance_in_time_universal(self,delta_t):
        # This method advances the object on its course by delta t in time using the universal time of fligt formulation. This means it should be usable for all kinds of orbits.
        
        # Solve for new X
        new_X = fsolve(lambda X : self.t_in_dep_of_X(X) - delta_t,delta_t)
        
    def advance_in_true_anomaly(self,delta_nu):
        # This method increases the true anomaly by a given input. It can be used to find equi-distant-angle points on the orbit for visualization purposes. It also updates parabolic anomaly and mean anomaly.

        self.true_anomaly_epoch = self.true_anomaly_epoch + delta_nu
        self.parabolic_anomaly = np.sqrt(self.parameter) * np.tan( self.true_anomaly / 2.0 )
        self.mean_anomaly = self.parameter * self.parabolic_anomaly + 1.0 / 3 * self.parabolic_anomaly**3
        
    def calculate_advance_in_true_anomaly(self,delta_nu):
        # This method advances the object on its course by delta nu in true anomaly and returns the new position. It is useful for calculating points on the orbit without actually advancing the object itself.

        new_true_anomaly_epoch = self.true_anomaly_epoch + delta_nu
        
        r = self.parameter / ( 1 + self.eccentricity * np.cos(new_true_anomaly_epoch))
        
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(new_true_anomaly_epoch),r * np.sin(new_true_anomaly_epoch),0],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(new_true_anomaly_epoch),self.eccentricity + np.cos(new_true_anomaly_epoch),0],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.

        rotation_matrix = np.array([[np.cos(self.argument_periapsis), - np.sin(self.argument_periapsis), 0],\
        [np.sin(self.argument_periapsis), - np.cos(self.argument_periapsis), 0],\
        [0, 0, 1]],float)

        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
   
