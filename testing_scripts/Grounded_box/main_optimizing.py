# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:07:10 2024

@author: marce
"""

from fields import *
from output import *
from vec3 import *
from potential_solver import *
from world import *
from species import *
from numba import njit
import time
import os

print("\033[H\033[J")

if __name__ == "__main__":
    
    # Define the project directory
    project_dir = "Grounded_Box_Final_Project_optimizing"
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
            
    # Timing the initialization
    start_time = time.time()
    # print("Initializing the world...")
    
    # Initialize the domain
    world = World(21, 21, 21)  # mesh size
    world.set_extents([-0.1, -0.1, 0.0], [0.1, 0.1, 0.2])  # start and end of world
    world.set_time(2e-9, 5)  # time step size and number of time steps
    
    # Measure time for world initialization
    init_world_time = time.time() - start_time
    print(f"World initialization: \t\t\t {init_world_time:.6f} seconds.")
    
    # Set up particle species
    species = [
        Species("O+", 16 * Const.AMU, Const.QE, world),
        Species("e-", Const.ME, -1 * Const.QE, world)
    ]
    
    # Timing the potential solver initialization
    start_time = time.time()
    # print("Initializing potential solver and solving initial potential...")
    
    # Initialize potential solver and solve initial potential
    solver = PotentialSolver(world, 10000, 1e-4)
    solver.solve()
    
    # Measure time for potential solver initialization and solving
    init_solver_time = time.time() - start_time
    print(f"Potential solver: \t\t\t\t {init_solver_time:.6f} seconds.")
    
    # Timing the electric field computation
    start_time = time.time()
    # print("Computing electric field...")
    
    # Compute electric field
    solver.compute_ef()
    
    # Measure time for electric field computation
    compute_ef_time = time.time() - start_time
    print(f"Electric field computation: \t {compute_ef_time:.6f} seconds.")
    
    # Define grid sizes for ions and electrons
    np_ions_grid = [21, 21, 21]  # Equivalent to int3 np_ions_grid = { 21,21,21 }
    np_eles_grid = [41, 41, 41]  # Equivalent to int3 np_eles_grid = { 41,41,41 }
    
    # Timing the particle loading
    start_time = time.time()
    # print("Loading particles using the quiet start method...")
    
    # Load particles using the quiet start method
    species[0].load_particles_box_qs(world.get_x0(), world.get_xm(), 1e11, np_ions_grid)  # Ions
    species[1].load_particles_box_qs(world.get_x0(), world.get_xc(), 1e11, np_eles_grid)  # Electrons
    
    # Measure time for particle loading
    load_particles_time = time.time() - start_time
    print(f"Particles loaded: \t\t\t\t {load_particles_time:.6f} seconds.")

    # Main loop
    while world.advance_time():
        print(f"\nLoop Number {world.get_ts()}:")

        loop_start_time = time.time()  # Start timing the loop iteration

        # Update particle velocities and positions
        update_particles_start = time.time()
        for sp in species:
            sp.advance()
            sp.compute_number_density()
        update_particles_time = time.time() - update_particles_start
        print(f"Update particles: \t\t\t\t {update_particles_time:.6f} seconds.")

        # Update charge density
        compute_charge_density_start = time.time()
        world.compute_charge_density(species)
        compute_charge_density_time = time.time() - compute_charge_density_start
        print(f"Compute charge density: \t\t {compute_charge_density_time:.6f} seconds.")
    
        # Solve the potential
        solve_potential_start = time.time()
        solver.solve()
        solve_potential_time = time.time() - solve_potential_start
        print(f"Solve potential: \t\t\t\t {solve_potential_time:.6f} seconds.")
    
        # Compute electric field
        compute_ef_start = time.time()
        solver.compute_ef()
        compute_ef_time = time.time() - compute_ef_start
        print(f"Compute electric field: \t\t {compute_ef_time:.6f} seconds.")
        
        #Output.screen_output(world, species)
        
        # Save data
        if world.get_ts() % 10 == 0 or world.is_last_time_step():
            save_data_start = time.time()
            ts = world.get_ts()
            filename = f"Particles_QS_ts{ts:05d}"
            Output.fields(world, species, filename, project_dir)
            Output.diag_output(world, species, project_dir)
            save_data_time = time.time() - save_data_start
            print(f"Save data: \t\t\t\t\t {save_data_time:.6f} seconds.")
        
        loop_time = time.time() - loop_start_time  # Measure the total loop iteration time
        print(f"Total loop iteration: \t\t\t {loop_time:.6f} seconds.")

    print(f"\nSimulation took {world.get_wall_time()} seconds.")
