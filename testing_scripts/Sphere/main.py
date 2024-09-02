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
import time
#from numba import jit
print("\033[H\033[J")

if __name__ == "__main__":
    
    # Define the project directory
    project_dir = "Sphere_EF_with_inlet2"
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
        
    print("Initializing the world...")
    # Initialize the domain
    world = World(21, 21, 41)  # mesh size
    world.set_extents([-0.1, -0.1, 0.0], [0.1, 0.1, 0.4])  # start and end of world
    world.set_time(1e-7, 401)  # time step size and number of time steps

    phi_sphere = -100
    world.add_sphere((0,0,0.15),0.05,phi_sphere)
    world.add_inlet()
    
    
    # # Set up particle species
    # species = [
    #     Species("O+", 16 * Const.AMU, Const.QE, world),
    #     Species("e-", Const.ME, -1 * Const.QE, world)
    # ]

   # Initialize the solver
    solver = PotentialSolver(world, 10000, 1e-4)
    
    # Set reference values
    solver.set_reference_values(0, 1.5, 1e12)
    
    # Solve the potential
    solver.solve()
    
    solver.compute_ef()
    # # Define grid sizes for ions and electrons
    # np_ions_grid = [41, 41, 41]  # Equivalent to int3 np_ions_grid = { 21,21,21 }
    # np_eles_grid = [21, 21, 21]  # Equivalent to int3 np_eles_grid = { 41,41,41 }
    
    # # Load particles using the quiet start method
    # species[0].load_particles_box_qs(world.get_x0(), world.get_xm(), 1e11, np_ions_grid)  # Ions
    # species[1].load_particles_box_qs(world.get_x0(), world.get_xc(), 1e11, np_eles_grid)  # Electrons
    # print('particles loaded')
    filename = "fields"
    Output.fields_without_species(world, filename, project_dir)

    # # Main loop
    # while world.advance_time():
    #     # print(f"Time step: {world.get_ts()}, Time: {world.get_time()}")
    #    # print("")
    #     # Update particle velocities and positions
    #     for sp in species:
    #         sp.advance()
    #         sp.compute_number_density()
    #     # Update charge density
    #     # print("Computing charge density...")
        
    #     world.compute_charge_density(species)
    #     # Assuming world.rho is a 3D NumPy array or a similar structure

    #     solver.solve()
        
    #     solver.compute_ef()
        
    #     Output.screen_output(world, species)
    #     # Save data
    #     # Save data
    #     if world.get_ts() % 10 == 0 or world.is_last_time_step():
    #         # Construct the filename with the time step and timestamp
    #         ts = world.get_ts()
            
    #         filename = f"Particles_QS_ts{ts:05d}"
            
    #         # Call the Output.fields method with the constructed filename
    #         Output.fields(world, species, filename, project_dir)
        
    #         # Call the diag_output method to save diagnostics data
    #         Output.diag_output(world, species, project_dir)
            
    # print(f"\nSimulation took {world.get_wall_time()} seconds.")
    print("\nProgram finished! :)")




