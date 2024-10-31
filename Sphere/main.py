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
    project_dir = "Sphere_particle_flow_2"
    
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
    
    species_list = []
    species_list.append(Species(name="O+", mass=16 * Const.AMU, charge=Const.QE, mpw0=1e2,world=world))
    
    # Set up injection sources
    sources = []
    sources.append(ColdBeamSource(species=species_list[0], world=world, v_drift=7000, den=1e10))

    # # Set up particle species
    # species = [
    #     Species("O+", 16 * Const.AMU, Const.QE, world),
    #     Species("e-", Const.ME, -1 * Const.QE, world)
    # ]

   # Initialize the solver
    solver = PotentialSolver(world, 10000, 1e-4)
    
    # Set reference values
    solver.set_reference_values(0, 1.5, 1e10)
    
    # Check if initial fields are already saved
    if os.path.exists(os.path.join(project_dir, 'phi.npy')) and os.path.exists(os.path.join(project_dir, 'ef.npy')):
       print("Loading initial fields from files...")
       world.load_phi(project_dir)
       world.load_ef(project_dir)
    else:
       print("Calculating initial fields...")
       # Solve the potential and compute electric field
       solver.solve()
       solver.compute_ef()
       
       # Save the calculated fields
       world.save_phi(project_dir)
       world.save_ef(project_dir)

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
        Output.screen_output(world, species_list)
        Output.diag_output(world, species,project_dir)
    
        # Periodically write out results
        if world.get_ts() % 20 == 0 or world.is_last_time_step():
            ts = world.get_ts()
            filename = f"Particles_QS_ts{ts:05d}"
            Output.fields(world, species,filename,project_dir)

    print(f"\nSimulation took {world.get_wall_time()} seconds.")
    print("\nProgram finished! :)")




