##############################
# Import necessary libraries #
##############################

from functions import *
import numpy as np
from scipy.optimize import fsolve

######################################
# Define class for celestial bodies. #
######################################
x
class celestial_body:
    # This class assumes a reference coordinate system such that a large mass is situated at the origin. It might actually assume some more things.
    
    ####### Init #######
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Technical properties
        self.mass = mass
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction

        # Position and Velocity
        self.position = position
        self.velocity = velocity

    def export_position_velocity(self):
        return self.position, self.velocity

    def advance_in_true_anomaly(self,delta_nu):
        pass

    def advance_in_time(self,delta_t):
        pass

    def export_orbit(self,number_points=60):
        pass
    
    def t_in_dep_of_X(self, X):
        pass

    def advance_in_time_universal(self,delta_t):
        pass

    def calculate_advance_in_true_anomaly(self,delta_nu):
        pass
    
    @classmethod
    def create_object(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determine type of orbit, planarity and initialize correct class of orbit
        h = np.cross(position,velocity) # Calculate angular momentum h
        energy = mass * np.dot(velocity,velocity) / 2.0 - mu / norm(position)
        eccentricity = norm(1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity))

        if norm(h) != 0:
            ### Non-collision ###
            if energy < 0:
                ### Elliptic/Cirvular ###
                if eccentricity > 0:
                    ### Elliptic ###
                    if h[2] == norm(h):
                        ### Planar ###
                        return celestial_body_on_planar_elliptic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                    else:
                        ### Non-planar ###
                        return celestial_body_on_nonplanar_elliptic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Circular ###
                    if h[2] == norm(h):
                        ### Planar ###
                        return celestial_body_on_planar_circular_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                    else:
                        ### Non-planar ###
                        return celestial_body_on_nonplanar_circular_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
            elif energy == 0:
                ### Parabolic ###
                if h[2] == norm(h):
                    ### Planar ###
                    return celestial_body_on_planar_parabolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Non-planar ###
                    return celestial_body_on_nonplanar_parabolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
            else:
                ### Hyperbolic ###
                if h[2] == norm(h):
                    ### Planar ###
                    return celestial_body_on_planar_hyperbolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Non-planar ###
                    return celestial_body_on_nonplanar_hyperbolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
        else:
            return celestial_body_on_collision_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
