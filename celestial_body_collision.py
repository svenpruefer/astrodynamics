##############################
# Import necessary libraries #
##############################

from functions import *
from celestial_body_base import *
import numpy as np
from scipy.optimize import fsolve
 
######################################
# Define class for colission orbits. #
######################################

class celestial_body_on_collision_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type. Should this class know wether it is planar?
        self.typ = "collision"
        
        self.inclination = np.arccos(np.dot(np.array([0,0,1],float),position) / norm(position)) # i
        self.energy = mass * np.dot(velocity,velocity) / 2.0 - mu / norm(position) # E        
        self.position = position
        self.velocity = velocity
        
    ####### Export #######
        
    def export_position_velocity(self):
        # Exports position and velocity of celestial body.
        return self.position, self.velocity
    
    def export_orbit(self):
        # Returns a list of three dimensional coordinates for the orbit. Perhaps plus most distant point?

        if self.energy >= 0:
            if np.dot(self.position,self.velocity) >= 0:
                position = np.zeros( (2,3) )
                position[0,:] = self.position
                position[1,:] = 5 * self.position
            else:
                position = np.zeros( (2,3) )
                position[0,:] = self.position
                position[1,:] = [0,0,0]
        else:
            if np.dot(self.position,self.velocity) < 0:
                position = np.zeros( (2,3) )
                position[0,:] = self.position
                position[1,:] = [0,0,0]
            else:
                maximum_distance = - self.energy / mu
                position = np.zeros( (2,3) )
                position[0,:] = self.position / norm(self.position) * maximum_distance
                position[1,:] = [0,0,0]
        
        return position # Orbit is non-periodic.
        
    ###### Advance along orbit #######
    
    def advance_in_time(self,delta_t):
        # This method advances the object on its course by delta t in time. It needs to integrate the equation of motion or solve it. How? Temporarily it will just use naive Euler integration.
        
        new_position = self.position + delta_t * self.velocity
        new_velocity = self.velocity + delta_t * self.mu / self.position**2
        self.position = new_position
        self.velocity = new_velocity
        
    def t_in_dep_of_X(self, X):
        # Is this still applicable for collision  orbits?
        r_0, v_0 = self.export_postion_velocity()
        return 1 / np.sqrt(self.mu) * ( np.dot(r_0,v_0) /np.sqrt(self.mu) * X**2 * C(X) + ( 1 - norm(r_0) / self.semi_major_axis ) * X**3 * S(X) + norm(r_0) * X )
    
    def advance_in_time_universal(self,delta_t):
        # This method advances the object on its course by delta t in time using the universal time of fligt formulation. Is this applicable to collision orbits?
        
        # Solve for new X
        new_X = fsolve(lambda X : self.t_in_dep_of_X(X) - delta_t,delta_t)
        
    def advance_in_true_anomaly(self,delta_nu):
        # What is this supposed to do for collision orbits?
        pass
    
    def calculate_advance_in_true_anomaly(self,delta_nu):
        # This method advances the object on its course by delta nu in true anomaly and returns the new position. It is useful for calculating points on the orbit without actually advancing the object itself.
        # How should this work?
        
        
        return position, velocity
