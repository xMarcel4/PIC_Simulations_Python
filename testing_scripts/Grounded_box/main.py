# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:38:32 2024

@author: marce
"""
from fields import *
from output import *
from vec3 import *
from potential_solver import *
from world import *
from species import *
print("\033[H\033[J")

if __name__ == "__main__":
    print("Initializing the world...")
    # Initialize the domain
    world = World(21, 21, 21)  # mesh size
    world.set_extents([-0.1, -0.1, 0.0], [0.1, 0.1, 0.2])  # start and end of world
    world.set_time(2e-10, 1000)  # time step size and number of time steps

    # Set up particle species
    species = [
        Species("O+", 16 * Const.AMU, Const.QE, world),
        Species("e-", Const.ME, -1 * Const.QE, world)
    ]

    # Initialize potential solver and solve initial potential
    solver = PotentialSolver(world, 10000, 1e-6)
    solver.solve()
    solver.compute_ef()

    # Create particles
    np_ions = 800000 # number of simulation ions
    np_eles = 100000  # number of simulation electrons
    species[0].load_particles_box_qs(world.get_x0(), world.get_xm(), 1e11, np_ions)  # ions
    species[1].load_particles_box_qs(world.get_x0(), world.get_xc(), 1e11, np_eles)  # electrons
    print('particles loaded')
    
    # Update fields
for index, sp in enumerate(species):
    print(f"Computing number density for species {index + 1}/{len(species)}: {sp.name}")
    sp.compute_number_density()


world.compute_charge_density(species)
solver.solve()
solver.compute_ef()

# Output results
filename = 'Particles_QS'
Output.fields(world, species,filename)
