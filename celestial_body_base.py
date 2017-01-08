##############################
# Import necessary libraries #
##############################

from functions import *
import numpy as np
from scipy.optimize import fsolve

##############################################
# Define general class for celestial bodies. #
##############################################

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
    def create_object(self,position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
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
                    if h[0] == 0 and h[1] == 0:
                        ### Planar ###
                        return celestial_body_on_planar_elliptic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                    else:
                        ### Non-planar ###
                        return celestial_body_on_nonplanar_elliptic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Circular ###
                    if h[0] == 0 and h[1] == 0:
                        ### Planar ###
                        return celestial_body_on_planar_circular_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                    else:
                        ### Non-planar ###
                        return celestial_body_on_nonplanar_circular_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
            elif energy == 0:
                ### Parabolic ###
                if h[0] == 0 and h[1] == 0:
                    ### Planar ###
                    return celestial_body_on_planar_parabolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Non-planar ###
                    return celestial_body_on_nonplanar_parabolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
            else:
                ### Hyperbolic ###
                if h[0] == 0 and h[1] == 0:
                    ### Planar ###
                    return celestial_body_on_planar_hyperbolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
                else:
                    ### Non-planar ###
                    return celestial_body_on_nonplanar_hyperbolic_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)
        else:
            return celestial_body_on_collision_orbit(position,velocity,mass,mu,ejection_speed,fuel_fraction)

    ###################
    # Applying thrust #
    ###################

    def apply_thrust(self,delta_v):
        # This function applies an acceleration of delta_v to the body. This happens instanteneous without actually moving the object. It actually returns the new object as the type or subclass might change.

        position, old_velocity = self.export_position_velocity()
        velocity = old_velocity + delta_v

        old_mass = self.mass
        mu = self.mu
        ejection_speed = self.ejection_speed
        fuel_fraction = self.fuel_fraction

        # Calculate new fuel_fraction and mass
        
        return celestial_body.create_object(position,velocity,new_mass,mu,ejection_speed,new_fuel_fraction)


############################################################
# Define class for non-planar elliptical celestial bodies. #
############################################################

class celestial_body_on_nonplanar_elliptic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "elliptical"
        self.planar = False

        h = np.cross(position,velocity) # Calculate angular momentum h
        n = np.cross(np.array([0,0,1],float),h) # Calculate node n
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass

        # Orbit definitions
        self.semi_major_axis = p / (1-np.dot(e,e))
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

        self.energy = - mu / ( 2.0 * self.semi_major_axis ) # E
        self.parameter = self.semi_major_axis * (1 - self.eccentricity**2) # p
        
        # Eccentric Anomaly
        if ( 0 <= self.true_anomaly_epoch ) and ( self.true_anomaly_epoch <= np.pi):
            self.eccentric_anomaly = np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E
        else:
            self.eccentric_anomaly = 2 * np.pi - np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E
                
        # Mean anomaly
        self.mean_anomaly = self.eccentric_anomaly - self.eccentricity * np.sin(self.eccentric_anomaly) # M
        self.mean_motion = np.sqrt(mu / self.semi_major_axis**3 ) # n

        # Period
        self.period = 2 * np.pi / np.sqrt(mu) * np.sqrt(self.semi_major_axis**3) # T
        
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
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
            ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points=60):
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
   

######################################################
# Define class for planar elliptic celestial bodies. #
#####################################################

class celestial_body_on_planar_elliptic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "elliptical"
        self.planar = True

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
        h = np.cross(position,velocity) # Calculate angular momentum h
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu
        
        self.semi_major_axis = p / (1-np.dot(e,e))
        self.eccentricity = norm(e)
                
        self.argument_periapsis = np.arccos(np.dot(np.array([1,0,0],float),e) / (norm(e)))
        
        if np.dot(position,velocity) >= 0:
            self.true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
        else:
            self.true_anomaly_epoch = 2 * np.pi - np.arccos(np.dot(e,position) / (norm(e) * norm(position)))

        self.energy = - mu / ( 2.0 * self.semi_major_axis ) # E
        self.parameter = self.semi_major_axis * (1 - self.eccentricity**2) # p
        
        # Eccentric Anomaly
        if ( 0 <= self.true_anomaly_epoch ) and ( self.true_anomaly_epoch <= np.pi):
            self.eccentric_anomaly = np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E
        else:
            self.eccentric_anomaly = 2 * np.pi - np.arccos((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # E
                
        # Mean anomaly
        self.mean_anomaly = self.eccentric_anomaly - self.eccentricity * np.sin(self.eccentric_anomaly) # M
        self.mean_motion = np.sqrt(mu / self.semi_major_axis**3 ) # n

        # Period
        self.period = 2 * np.pi / np.sqrt(mu) * np.sqrt(self.semi_major_axis**3) # T
        
        #Universal time of flight
        self.X = 0 # X
    
        ####### Export #######
        
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?

        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))
            
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch)],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch)],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.

        rotation_matrix = np.array([[np.cos(self.argument_periapsis), - np.sin(self.argument_periapsis)],\
        [np.sin(self.argument_periapsis), np.cos(self.argument_periapsis)]],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points=60):
        # Returns a list of three dimensional coordinates for the orbit.
        position = np.zeros( (number_points,2) )
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
        new_eccentric_anomaly = fsolve(lambda E : E - self.eccentricity * np.sin(E) - new_mean_anomaly,new_mean_anomaly) 
        
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
        return 1 / np.sqrt(self.mu) * ( np.dot(r_0,v_0) / np.sqrt(self.mu) * X**2 * C(X) + ( 1 - norm(r_0) / self.semi_major_axis ) * X**3 * S(X) + norm(r_0) * X )
    
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
        
        r = self.parameter / ( 1 + self.eccentricity * np.cos(new_true_anomaly_epoch))
        
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(new_true_anomaly_epoch),r * np.sin(new_true_anomaly_epoch)],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(new_true_anomaly_epoch),self.eccentricity + np.cos(new_true_anomaly_epoch)],float)
        
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.

        rotation_matrix = np.array([[np.cos(self.argument_periapsis), - np.sin(self.argument_periapsis)],\
        [np.sin(self.argument_periapsis),np.cos(self.argument_periapsis)]],float)

        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
   

############################################################
# Define class for non-planar hyperbolic celestial bodies. #
############################################################

class celestial_body_on_nonplanar_hyperbolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "hyperbolic"
        self.planar = False

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
        h = np.cross(position,velocity) # Calculate angular momentum h
        n = np.cross(np.array([0,0,1],float),h) # Calculate node n
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu
        
        self.semi_major_axis = p / (1-np.dot(e,e))
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

        self.energy = - mu / ( 2.0 * self.semi_major_axis ) # E
        self.parameter = self.semi_major_axis * (1 - self.eccentricity**2) # p

        # Hyperbolic Anomaly
        if ( 0 <= self.true_anomaly_epoch ) and ( self.true_anomaly_epoch <= np.pi):
            self.hyperbolic_anomaly = np.arccosh((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # F
        else:
            self.hyperbolic_anomaly = - np.arccosh((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # F
                
        # Hyperbolic mean anomaly
        self.mean_anomaly = - self.hyperbolic_anomaly + self.eccentricity * np.sinh(self.hyperbolic_anomaly) # M
        self.mean_motion = np.sqrt(mu / ( (-self.semi_major_axis)**3 ) ) # n

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
        rotation_matrix = np.array([[np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) - np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , np.sin(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.longitude_ascending_node) * np.cos(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.sin(self.argument_periapsis) * np.cos(self.inclination) , - np.sin(self.longitude_ascending_node) * np.sin(self.argument_periapsis) + np.cos(self.longitude_ascending_node) * np.cos(self.argument_periapsis) * np.cos(self.inclination) , - np.cos(self.longitude_ascending_node) * np.sin(self.inclination)],\
        [np.sin(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.argument_periapsis) * np.sin(self.inclination) , np.cos(self.inclination)]\
            ],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)
        
        return position, velocity
    
    def export_orbit(self,number_points=60):
        # Returns a list of three dimensional coordinates for the orbit.

        interval = 2 * np.pi / number_points
        points = [2 * np.pi / number_points * i for i in range(- number_points // 2, number_points // 2) if np.cos(i * interval) > -1 / self.eccentricity]
        position = np.array([self.calculate_advance_in_true_anomaly(point)[0] for point in points],float)
        return position
        
    ###### Advance along orbit #######
    
    def advance_in_time(self,delta_t):
        # This method advances the object on its course by delta t in time. This means that it needs to translate the time difference into changes in the true anomaly at epoch and then add this number to the existing value.
            # delta_t should be small enough such that the body does not evolve more than one period. Is this necessary?
        
        # Update mean anomaly. Ignore full rotations.
        new_mean_anomaly = self.mean_motion * delta_t + self.mean_anomaly
        
        # Solve e*sinh(F)-F=M numerically
        new_hyperbolic_anomaly = fsolve(lambda F : self.eccentricity * np.sinh(F) - F - new_mean_anomaly, new_mean_anomaly) 
        
        # Calculate new true anomaly at epoch
        if new_hyperbolic_anomaly >= 0:
            new_true_anomaly_epoch = np.arccos( ( np.cosh(new_hyperbolic_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cosh(new_hyperbolic_anomaly)))
        else:
            new_true_anomaly_epoch = 2 * np.pi - np.arccos( ( np.cosh(new_hyperbolic_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cosh(new_eccentric_anomaly)))
            
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
            self.hyperbolic_anomaly = np.arccosh( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))
        else:
            self.hyperbolic_anomaly = 2 * np.pi - np.arccosh( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))

        self.mean_anomaly = - self.hyperbolic_anomaly + self.eccentricity * np.cosh( self.hyperbolic_anomaly )
        
    def calculate_advance_in_true_anomaly(self,nu):
        # This method advances the object on its course by delta nu in true anomaly and returns the new position. It is useful for calculating points on the orbit without actually advancing the object itself. In contrast to the elliptic case, the argument is not added to the true anomaly but it is in fact the new value. This is because hyperbolic orbits are not periodic. 
        
        new_true_anomaly_epoch = nu
        
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
# Define class for planar hyperbolic celestial bodies. #
#####################################################

class celestial_body_on_planar_hyperbolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "hyperbolic"
        self.planar = True

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
        h = np.cross(position,velocity) # Calculate angular momentum h
        e = 1.0 / mu * ((np.dot(velocity,velocity) - mu / norm(position)) * position - np.dot(position,velocity) * velocity) # Calculate eccentricity vector pointing in direction of perihelion
        p = np.dot(h,h) / mu

        if h[2] > 0:
            self.direction = 'positive'
        else:
            self.direction = 'negative'
            
        self.semi_major_axis = p / (1-np.dot(e,e))
        self.eccentricity = norm(e)
                
        self.argument_periapsis = np.arccos(np.dot(np.array([1,0,0],float),e) / (norm(e)))

        if np.dot(position,velocity) >= 0:
            self.true_anomaly_epoch = np.arccos(np.dot(e,position) / (norm(e) * norm(position)))
        else:
            self.true_anomaly_epoch = 2 * np.pi - np.arccos(np.dot(e,position) / (norm(e) * norm(position)))

        self.energy = - mu / ( 2.0 * self.semi_major_axis ) # E
        self.parameter = self.semi_major_axis * (1 - self.eccentricity**2) # p

        
        # Hyperbolic Anomaly
        if ( 0 <= self.true_anomaly_epoch ) and ( self.true_anomaly_epoch <= np.pi):
            self.hyperbolic_anomaly = np.arccosh((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # F
        else:
            self.hyperbolic_anomaly = - np.arccosh((self.eccentricity + np.cos(self.true_anomaly_epoch)) / (1 + self.eccentricity * np.cos(self.true_anomaly_epoch))) # F
                
        # Hyperbolic mean anomaly
        self.mean_anomaly = - self.hyperbolic_anomaly + self.eccentricity * np.sinh(self.hyperbolic_anomaly) # M
        self.mean_motion = np.sqrt(mu / ( (-self.semi_major_axis)**3 ) ) # n
        
        #Universal time of flight
        self.X = 0 # X
    
        ####### Export #######
        
    def export_position_velocity(self):
        # Exports position and velocity of celestial body. How should time dependence be incorparated? Should it be a parameter for this function?

        r = self.parameter  / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch))
            
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(self.true_anomaly_epoch),r * np.sin(self.true_anomaly_epoch)],float)
        if self.direction == 'positive':
            velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch)],float)
        else:
            velocity_perifocal_system = - np.sqrt(self.mu / self.parameter) * np.array([-np.sin(self.true_anomaly_epoch),self.eccentricity + np.cos(self.true_anomaly_epoch)],float)
                           
        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.

        rotation_matrix = np.array([[np.cos(self.argument_periapsis), - np.sin(self.argument_periapsis)],\
        [np.sin(self.argument_periapsis), np.cos(self.argument_periapsis)]],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)

        return position, velocity

    
    def export_orbit(self,number_points=60):
        # Returns a list of three dimensional coordinates for the orbit.

        interval = 2 * np.pi / number_points
        points = [2 * np.pi / number_points * i for i in range(- number_points // 2, number_points // 2) if np.cos(i * interval) > -1 / self.eccentricity]
        position = np.array([self.calculate_advance_in_true_anomaly(point)[0] for point in points],float)
        return position

    ###### Advance along orbit #######

    
    def advance_in_time(self,delta_t):
        # This method advances the object on its course by delta t in time. This means that it needs to translate the time difference into changes in the true anomaly at epoch and then add this number to the existing value.
            # delta_t should be small enough such that the body does not evolve more than one period. Is this necessary?
        
        # Update mean anomaly. Ignore full rotations.
        new_mean_anomaly = self.mean_motion * delta_t + self.mean_anomaly
        
        # Solve e*sinh(F)-F=M numerically
        new_hyperbolic_anomaly = fsolve(lambda F : self.eccentricity * np.sinh(F) - F - new_mean_anomaly, new_mean_anomaly) 
        
        # Calculate new true anomaly at epoch
        if new_hyperbolic_anomaly >= 0:
            new_true_anomaly_epoch = np.arccos( ( np.cosh(new_hyperbolic_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cosh(new_hyperbolic_anomaly)))
        else:
            new_true_anomaly_epoch = 2 * np.pi - np.arccos( ( np.cosh(new_hyperbolic_anomaly) - self.eccentricity ) / ( 1 - self.eccentricity * np.cosh(new_eccentric_anomaly)))
            
        # Update values of true anomaly at epoch and eccentric anomaly and mean anomaly
        self.true_anomaly_epoch = new_true_anomaly_epoch
        self.mean_anomaly = new_mean_anomaly
        self.eccentric_anomaly = new_eccentric_anomaly
        
    def t_in_dep_of_X(self, X):
        r_0, v_0 = self.export_postion_velocity()
        return 1 / np.sqrt(self.mu) * ( np.dot(r_0,v_0) / np.sqrt(self.mu) * X**2 * C(X) + ( 1 - norm(r_0) / self.semi_major_axis ) * X**3 * S(X) + norm(r_0) * X )
    
    def advance_in_time_universal(self,delta_t):
        # This method advances the object on its course by delta t in time using the universal time of fligt formulation. This means it should be usable for all kinds of orbits.
        
        # Solve for new X
        new_X = fsolve(lambda X : self.t_in_dep_of_X(X) - delta_t,delta_t)

    def advance_in_true_anomaly(self,delta_nu):
        # This method increases the true anomaly by a given input. It can be used to find equi-distant-angle points on the orbit for visualization purposes. It also updates eccentric anomaly and mean anomaly.
        self.true_anomaly_epoch = self.true_anomaly_epoch + delta_nu
        if self.true_anomaly_epoch <= np.pi:
            self.hyperbolic_anomaly = np.arccosh( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))
        else:
            self.hyperbolic_anomaly = 2 * np.pi - np.arccosh( ( np.cos(self.true_anomaly_epoch) + self.eccentricity ) / ( 1 + self.eccentricity * np.cos(self.true_anomaly_epoch)))

        self.mean_anomaly = - self.hyperbolic_anomaly + self.eccentricity * np.cosh( self.hyperbolic_anomaly )
    
    def calculate_advance_in_true_anomaly(self,nu):
        # This method advances the object on its course by delta nu in true anomaly and returns the new position. It is useful for calculating points on the orbit without actually advancing the object itself. In contrast to the elliptic case, the argument is not added to the true anomaly but it is in fact the new value. This is because hyperbolic orbits are not periodic. 
        
        new_true_anomaly_epoch = nu
        
        r = self.parameter / ( 1 + self.eccentricity * np.cos(new_true_anomaly_epoch))
        
        # The perifocal coordinate system uses coordinate axes P, Q, W in this order, where P points in the direction of the periapsis and Q is perpendicular in positive direction in the plane of the orbit.

        position_perifocal_system = np.array([r * np.cos(new_true_anomaly_epoch),r * np.sin(new_true_anomaly_epoch)],float)
        velocity_perifocal_system = np.sqrt(self.mu / self.parameter) * np.array([-np.sin(new_true_anomaly_epoch),self.eccentricity + np.cos(new_true_anomaly_epoch)],float)

        # Calculate the rotation matrix from perifocal to fixed frame. Bate says, one should avoid this technique.
        
        rotation_matrix = np.array([[np.cos(self.argument_periapsis),  - np.sin(self.argument_periapsis)],\
                                    [np.sin(self.argument_periapsis), np.cos(self.argument_periapsis)]],float)
        
        position = np.dot(rotation_matrix,position_perifocal_system)
        velocity = np.dot(rotation_matrix,velocity_perifocal_system)

        return position, velocity
   

###########################################################
# Define class for non-planar parabolic celestial bodies. #
###########################################################

class celestial_body_on_nonplanar_parabolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "parabolic"
        self.planar = False

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
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
    
    def export_orbit(self,number_points=61):
        # Returns a list of three dimensional coordinates for the orbit. Remove point opposite perigree.

        if number_points % 2 == 1:
            number_points += 1
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
   

#######################################################
# Define class for planar parabolic celestial bodies. #
#######################################################

class celestial_body_on_planar_parabolic_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type
        self.typ = "parabolic"
        self.planar = True

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
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
    
    def export_orbit(self,number_points=61):
        # Returns a list of three dimensional coordinates for the orbit. Removes point opposite of perigree by using odd number of points.

        if number_points % 2 == 0:
            number_of_points += 1
        
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
   
######################################
# Define class for collision orbits. #
######################################

class celestial_body_on_collision_orbit(celestial_body):
    def __init__(self, position, velocity, mass = 1, mu = 1, ejection_speed = 1, fuel_fraction = 1):
        # Determination of orbit type. Should this class know wether it is planar?
        self.typ = "collision"

        # Technical definitions
        self.mu = mu
        self.ejection_speed = ejection_speed
        self.fuel_fraction = fuel_fraction
        self.mass = mass
        
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
