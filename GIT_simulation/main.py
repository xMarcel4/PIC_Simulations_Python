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
from Source import *
#from numba import jit
print("\033[H\033[J")

if __name__ == "__main__":
    
    # Define the project directory
    project_dir = "Test_without_phi_rho_2"
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
        
    print("Initializing the world...")
     
    
    # Convert the dimensions to meters (since SI units are typically used in physics simulations)
    x_length = 2.64e-3  # 2.6 mm in meters
    y_length = 2.64e-3  # 2.6 mm in meters
    z_length = 7.92e-3  # 6.0 mm in meters
    
    # Define the start and end points for the world extents
    x_start = -x_length / 2  # Centering the world around the origin
    x_end = x_length / 2
    y_start = -y_length / 2  # Centering the world around the origin
    y_end = y_length / 2
    z_start = 0.0  # Starting from 0 to positive z-axis
    z_end = z_length
    
    f = 3
    # Create the world object with the specified mesh size
    world = World(21*f, 21*f, 61*f)  # mesh size
    
    # Set the world extents using the calculated start and end points
    world.set_extents([x_start, y_start, z_start], [x_end, y_end, z_end])

    # Initialize the domain
    # world = World(21, 21, 41)  # mesh size
    # world.set_extents([-0.1, -0.1, 0.0], [0.1, 0.1, 0.4])  # start and end of world
    world.set_time(1e-7, 401)  # time step size and number of time steps

    # Add the inlet
    world.add_inlet()
    
    #plane worked perfectly
    # world.add_plane(z_start=3e-3, thickness=1e-3, phi_plane=-50)
    # world.add_plane_with_hole(z_start=3e-3, thickness=1e-3, phi_plane=-50, hole_center=(0, 0), hole_radius=4e-4)

    
    # Add three planes with holes
    world.add_plane_with_hole(z_start=1.5e-3, thickness=0.5e-3, phi_plane=1210, hole_center=(0, 0), hole_radius=0.6e-3)
    
    world.add_plane_with_hole(z_start=2.75e-3, thickness=1e-3, phi_plane=-100, hole_center=(0, 0), hole_radius=0.6e-3)
    
    world.add_plane_with_hole(z_start=4.51e-3, thickness=1e-3, phi_plane=0, hole_center=(0, 0), hole_radius=0.4e-3)
    
    # Apply all the objects to the fields
    # world.apply_objects_to_fields()
   # Initialize the solver
    solver = PotentialSolver(world, 10000, 1e-4)
    
    # Set reference values
    solver.set_reference_values(0, 1.5, 1e10)
    
    # Check if initial fields are already saved
    print("Calculating initial fields...")
    # Solve the potential and compute electric field
    #solver.solve()
    #solver.compute_ef()
       
    print("Initial fields are set. Proceeding with the simulation...")
   # Add the rest of your simulation loop here
    
    # # # Main loop
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
        Output.screen_output(world, species_list)
        Output.diag_output(world, species,project_dir)
    
        # Periodically write out results
        if world.get_ts() % 20 == 0 or world.is_last_time_step():
            ts = world.get_ts()
            filename = f"Particles_QS_ts{ts:05d}"
            Output.fields(world, species,filename,project_dir)

    filename = "fields"
    Output.fields_without_species(world,filename,project_dir)
    print(f"\nSimulation took {world.get_wall_time()} seconds.")
    print("\nProgram finished! :)")



