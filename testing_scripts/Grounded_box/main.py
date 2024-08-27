# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:38:32 2024

@author: marce
"""
import functions as f
import time
print("\033[H\033[J")
 

if __name__ == "__main__":
    start_time = time.time()

    # Initialize the domain
    world = f.World(21, 21, 21)
    world.set_extents([-0.1, -0.1, 0.0], [0.1, 0.1, 0.2])
    
    # Set phi[i=0] = 1
    for j in range(world.nj):
        for k in range(world.nk):
            world.phi.w_at(0, j, k, 1)

    # Set phi[k=0] = 2
    for i in range(world.ni):
        for j in range(world.nj):
            world.phi.w_at(i, j, 0, 2)
    
    # Initialize and solve potential
    solver = f.PotentialSolver(world, max_it=1000, tol=1e-6)
    solver.solve()
    solver.compute_ef()
    
    # Output results
    f.Output.create_results_dir()
    f.Output.fields(world)
    
    
    elapsed_time =  time.time() - start_time
    
    # Print the elapsed time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")


    # Print the fields and properties
    # print("rho",field_output(world.rho), "\n")
    # print("phi",field_output(world.phi), "\n")
    # print("ef",field_output(world.ef), "\n")
    # print("node_vol",field_output(world.node_vol), "\n")

    # print("x0\t", world.get_x0(), "\n")
    # print("xc\t", world.get_xc(), "\n")
    # print("xm\t", world.get_xm(), "\n")
    # print("dh\t", world.get_dh(), "\n")
    # print("nn\t", world.get_nn(), "\n")
