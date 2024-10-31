# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 03:00:32 2024

@author: marce
"""
from fields import *
from output import *
from vec2 import *
from potential_solver import *
from species import *
import time
from Source import *

from fields_2d import *
from output_2d import *
from world_2d import *
from potential_solver_2d import *
from species_2d import *
from Source_2d import *

import os
import csv
import time

def frange(start, stop, step):
    while start < stop:
        yield start
        start += step

if __name__ == "__main__":
    
    # Define the project directory
    project_dir = "TEST_overrelaxation_factor_2"
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    
    # Define the range of w values (relaxation factor)
    w_values = [round(x, 2) for x in list(frange(1.95, 1.99, 0.05))]  # Generate values from 1 to 1.95 in steps of 0.05
    
    for w in w_values:
        print(f"Running simulation for w = {w}...")

        # Create a unique CSV file for this relaxation value
        csv_filename = os.path.join(project_dir, f'iteration_times_w_{w:.2f}.csv')
        
        # Create a CSV file to store results for this relaxation value
        with open(csv_filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['w', 'Main Loop Time (s)', 'Iterations for Solve'])

            # Initialize the world
            world = World2D(10, 31)
            world.set_extents([x_start, y_start], [x_end, y_end])
            world.set_time(1e-8, 2)

            # Set up planes and inlets as before
            scale = 1e3
            world.add_plane_with_hole(y_position=1.5e-3, thickness=0.5e-3, phi_plane=1200/scale, hole_half_width=0.96e-3)
            world.add_plane_with_hole(y_position=2.75e-3, thickness=1e-3, phi_plane=-100/scale, hole_half_width=0.6e-3)
            world.add_plane_with_hole(y_position=4.51e-3, thickness=1e-3, phi_plane=175/scale, hole_half_width=0.96e-3)
            world.add_inlet()
            
            # Set up species and sources
            species_factor = 1e2
            species_list = [Species2D(name="Xe+", mass=131.293 * Const.AMU/species_factor, charge=Const.QE/species_factor, mpw0=1e2, world=world)]
            sources = [ColdBeamSource2D(species=species_list[0], world=world, v_drift=5000, den=1e13)]

            # Initialize the solver with current w value
            solver = PotentialSolver2D(world, 100000, 1e-4, w)
            solver.set_reference_values(0, 1.5, 1e10)
            
            # Record the start time

            # Solve the initial potential
            solver_iterations = solver.solve()
            solver.compute_ef()

            # Run the main simulation loop
            while world.advance_time():
                start_time = time.time()

                # Inject particles
                for source in sources:
                    source.sample()
                
                # Move particles and compute number density
                for sp in species_list:
                    sp.advance()
                    sp.compute_number_density()

                # Compute charge density
                world.compute_charge_density(species_list)

                # Update potential and electric field
                solver_iterations = solver.solve()
                solver.compute_ef()
            
                # Record the end time
                end_time = time.time()
                main_loop_time = end_time - start_time

                # Write the results to the CSV file
                writer.writerow([w, main_loop_time, solver_iterations])

    print("Simulation completed for all w values.")
