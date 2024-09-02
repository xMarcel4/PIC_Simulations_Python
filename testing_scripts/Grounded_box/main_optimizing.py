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

    
    # Initialize an array to store the results
    results = []
    
    # Define the range for w
    w_values = np.arange(1.3, 1.7, 0.1)  # 1.2, 1.3, ..., 1.6
    
    # Loop over each w value
    for w in w_values:
        print(f"Running simulation with w = {w:.1f}")
        
        solve_times = []
        
        # Initialize the solver with the current w value
        solver = PotentialSolver(world, max_it=1000, tol=1e-5, w=w)
        
        while world.advance_time():
            for sp in species:
                sp.advance()
                sp.compute_number_density()
            
            world.compute_charge_density(species)
    
            # Measure the time taken by solver.solve()
            start_time = time.time()
            solver.solve()
            end_time = time.time()
            solve_duration = end_time - start_time
            solve_times.append(solve_duration)
            
            solver.compute_ef()
            Output.screen_output(world, species)
            
            # Save data periodically
            if world.get_ts() % 10 == 0 or world.is_last_time_step():
                ts = world.get_ts()
                filename = f"Particles_QS_ts{ts:05d}"
                Output.fields(world, species, filename, project_dir)
                Output.diag_output(world, species, project_dir)
        
        # Calculate the average solve time for this value of w
        avg_solve_time = np.mean(solve_times)
        
        # Store the result in the results array
        results.append((w, avg_solve_time))
        
        # Reset or reinitialize the world, species, and other variables as necessary
        # (This depends on your specific setup, e.g., reloading initial conditions, resetting time, etc.)
    
    # Print out the results
    print("\nResults:")
    for w, avg_time in results:
        print(f"w = {w:.1f}, Average solve time = {avg_time:.6f} seconds")
    
    # Optionally, you can save the results array to a file
    np.savetxt("solve_times.txt", results, header="w value, Average solve time")
