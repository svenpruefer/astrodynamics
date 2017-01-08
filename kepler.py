##########################################
# Import necessary classes and libraries #
##########################################

from celestial_body_base import *
#from celestial_body_elliptic import *
#from celestial_body_hyperbolic import *
#from celestial_body_parabolic import *
#from celestial_body_collision import *
#from celestial_body_create import *
import numpy as np
from scipy.optimize import fsolve


##################################
# Define Function Kepler Problem #
##################################

def kepler_problem(mass_1,position_1,velocity_1,mass_2,position_2,velocity_2):
    # Function takes initial data and parameters for two celestial objects, calculates their common center of mass and then returns two instances of the class celestia_object corresponding to those two objects in the intertial system translated to the common center of mass.
    
    center_of_mass = ( mass_1 * position_1 + mass_2 * position_2 ) / ( mass_1 + mass_2 )
    center_of_mass_velocity = ( mass_1 * velocity_1 + mass_2 * velocity_2 ) / ( mass_1 + mass_2 )
    mu = 1.0 * ( mass_1 + mass_2 ) # Gravitation constant = 1
    
    relative_position = position_1 - position_2
    relative_velocity = velocity_1 - velocity_2
    
    total_body = celestial_body.create_object(relative_position,relative_velocity)
    
    #position_1_in_center_of_mass_system = mu / mass_1 * total_body.export_position_velocity()[0]
    #position_2_in_center_of_mass_system = - mu / mass_2 * total_body.export_position_velocity()[0]
    
    #velocity_1_in_center_of_mass_system = mu / mass_1 * total_body.export_position_velocity()[1]
    #velocity_2_in_center_of_mass_system = - mu / mass_2 * total_body.export_position_velocity()[1]
    
    #position_1_in_center_of_mass_system = position_1 - center_of_mass #mu / mass_1 * relative_position
    #position_2_in_center_of_mass_system = position_2 - center_of_mass #- mu / mass_2 * relative_position
    #velocity_1_in_center_of_mass_system = velocity_1 - center_of_mass_velocity
    #velocity_2_in_center_of_mass_system = velocity_2 - center_of_mass_velocity
    
    #print(position_1_in_center_of_mass_system)
    #print(velocity_1_in_center_of_mass_system)
    #print(position_2_in_center_of_mass_system)
    #print(velocity_2_in_center_of_mass_system)
    
    #body_1 = celestial_body.from_position_velocity(mass_1,mu,position_1_in_center_of_mass_system,velocity_1_in_center_of_mass_system)
    #body_2 = celestial_body.from_position_velocity(mass_2,mu,position_2_in_center_of_mass_system,velocity_2_in_center_of_mass_system)
    
    return total_body#body_1, body_2  
    

