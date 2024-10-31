# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:38:32 2024

@author: marce
"""
from fields import *
from output import *
from vec2 import *
from potential_solver import *
# from world import World2D
from species import *
import time
from Source import *

from fields_2d import *
from output_2d import *
from world_2d import *
from potential_solver_2d import *
from species_2d import *
from Source_2d import *

#from numba import jit
print("\033[H\033[J")

if __name__ == "__main__":
    
    # Define the project directory
    project_dir = "C:/Users/Marcel/OneDrive/Dokumente/GitHub_repositories/GIT_particle_trajectory_simulation/testing_scripts/GIT_simulation_2D/Potential_Layout"
    
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
        
    print("Initializing the world...")
     
    # Convert the dimensions to meters (since SI units are typically used in physics simulations)
    x_length = 2.64e-3  # 2.64 mm in meters (X-dimension)
    y_length = 5.7e-3  # 7.92 mm in meters (Y-dimension for 2D)
    
    # Define the start and end points for the world extents
    x_start = -x_length / 2  # Centering the world around the origin in X
    x_end = x_length / 2
    y_start = 0.0  # Start from 0 for the positive Y-axis
    y_end = y_length  # Total height in the Y-axis

    
    # Create the 2D world object with the specified mesh size (scaled mesh density)
    world = World2D(10*4, 22*4)  # 2D mesh size with 'ni' and 'nj'
    # world = World2D(10, 30)  # 2D mesh size with 'ni' and 'nj'

    
    # Set the world extents using the calculated start and end points
    world.set_extents([x_start, y_start], [x_end, y_end])

    world.set_time(7.1e-9, 500)  # time step size and number of time steps
    
    # phi_circle = -100
    # world.add_circle(center=(0, 0.15), radius=0.05, phi_circle=phi_circle)
    
    # world.add_plane(y_position=1.5e-3, thickness=0.5e-3, phi_plane=1200/scale)
    # world.add_plane(y_position=2.75e-3, thickness=1e-3, phi_plane=-100/scale)
    # world.add_plane(y_position=4.51e-3, thickness=1e-3, phi_plane=175/scale)

    scale = 1e0

    world.add_plane_with_hole(y_position=1.5e-3,thickness=0.5e-3,phi_plane=1200/scale,hole_half_width=0.96e-3)
    
    # Plane 2 with a hole
    world.add_plane_with_hole(y_position=2.75e-3,thickness=1e-3,phi_plane=-100/scale,hole_half_width=0.6e-3)
    
    # Plane 3 with a hole
    world.add_plane_with_hole(y_position=4.51e-3,thickness=1e-3, phi_plane=10/scale,hole_half_width=0.96e-3)  

    # Add the inlet
    world.add_inlet()
    
    species_factor = 1.5
    w = 1.9
    species_list = []
    #species_list.append(Species2D(name="O+", mass=1 * Const.AMU, charge=Const.QE, mpw0=1e2,world=world))
    species_list.append(Species2D(name="Xe+", mass=131.293 * Const.AMU/species_factor, charge=Const.QE/species_factor, mpw0=1e2, world=world))

    
    # Set up injection sources
    sources = []
    sources.append(ColdBeamSource2D(species=species_list[0], world=world, v_drift=7000, den=1e12))
    
    solver = PotentialSolver2D(world, 100000, 1e-4,w)
    
    # # Set reference values
    solver.set_reference_values(0, 1.5, 1e10)
    
    # print("Calculating initial fields...")
    solver.solve(project_dir,save_data=True)
    solver.compute_ef()
       
    print("Initial fields are set. Proceeding with the simulation...")
   # Add the rest of your simulation loop here
    
    # # Main loop
    while world.advance_time():
        # Inject particles
        for source in sources:
            source.sample()
            
        # Move particles and compute number density
        for sp in species_list:
            sp.advance()
            sp.compute_number_density()
        
        # Compute charge density
        world.compute_charge_density(species_list) 
    
        # Update potential
        solver.solve()
    
        # Obtain electric field
        solver.compute_ef()
    
        # Screen and file output
        Output2D.screen_output(world, species_list)
        Output2D.diag_output(world, species_list,project_dir)
    
        # Periodically write out results
        if world.get_ts() % 20 == 0 or world.is_last_time_step():
            ts = world.get_ts()
            filename = f"Particles_QS_ts{ts:05d}"
            Output2D.fields(world, species_list,filename,project_dir)
            
    
    # filename = "fields"
    Output2D().fields_without_species_2d(world,filename,project_dir)
    Output2D().fields(world,species_list,filename,project_dir)
    print(f"\nSimulation took {world.get_wall_time()} seconds.")
    print("\nProgram finished! :)")



