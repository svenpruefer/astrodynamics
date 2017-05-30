##############################
# Import necessary libraries #
##############################

from functions import *
from celestial_body_base import *
import numpy as np
from scipy.optimize import fsolve

##############################
# Change system of reference #
##############################

def copy_system_of_reference(base_body,old_center_body):
    # This function takes two celestial bodies, assuming that base_body's orbital elements are defined with respect to old_center_body which in turn has orbital elements defined with respect to some new_center_body. Furthermore, it assumes that the reference orthonormal coordinate systems are parallel.
    # It then calculates the orbital elements of base_body in terms of those of new_center_body. Also the case new_center_body = [0,0,0,0,0,0] needs to be dealt with separately.

    # Fixme: Origin ist only positional origin at the moment. Is there a celestial object "origin" ?
    if old_center_body.export_position_velocity()[0] != [0,0,0]:
        new_base_body = celestial_body.create_object(base_body.export_position_velocity()[0]+old_center_body.export_position_velocity()[0],base_body.export_position_velocity()[1]+old_center_body.export_position_velocity()[1])
    else:
        new_base_body = base_body

    return new_base_body

def change_system_of_reference(base_body,new_center_body):
    # This function takes two celestial bodies, assuming that both base_body's and new_center_body's orbital elements are defined with respect to the same old_center_body. Furthermore, it assumes that the reference orthonormal coordinate systems are parallel.
    # It then calculates the orbital elements of base_body in terms of those of new_center_body. Also the case new_center_body = [0,0,0,0,0,0] needs to be dealt with separately as this is ust the identity transformation.

    # Fixme: Origin ist only positional origin at the moment. Is there a celestial object "origin" ?
    if new_center_body.export_position_velocity()[0:3] != [0,0,0]:
        new_base_body = celestial_body.create_object(base_body.export_position_velocity()[0]-old_center_body.export_position_velocity()[0],base_body.export_position_velocity()[1]+old_center_body.export_position_velocity()[1])
    else:
        new_base_body = base_body

    return new_base_body


