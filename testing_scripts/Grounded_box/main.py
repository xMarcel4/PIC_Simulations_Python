# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:38:32 2024

@author: marce
"""

from vec3 import Vec3
from fields_2 import Field
from fields_2 import *
from world import World
from world import *
#print("\033[H\033[J")


if __name__ == "__main__":
    # Initialize the World object with dimensions 4x4x4
    world = World(4, 4, 4)
    
    # Set the extents of the world
    world.set_extents([-1, -1, 0], [2, 2, 3])

    # Print the fields and properties
    print("rho\n",field_output(world.rho), "\n")
    print("phi\n",field_output(world.phi), "\n")
    print("ef\n",field_output(world.ef), "\n")
    print("node_vol\n",field_output(world.node_vol), "\n")

    print("x0\t", world.get_x0(), "\n")
    print("xc\t", world.get_xc(), "\n")
    print("xm\t", world.get_xm(), "\n")
    print("dh\t", world.get_dh(), "\n")
    print("nn\t", world.nn, "\n")
