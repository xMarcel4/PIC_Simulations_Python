# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:44:56 2024

@author: marce
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:00:43 2024

@author: marce (adapted for 2D)
"""
from fields_2d import *
from world_2d import *
import numpy as np
from math import sqrt, exp
import time

import os

import decimal
from decimal import Decimal

class PotentialSolver2D:
    def __init__(self, world, max_solver_it, tolerance, w):
        self.world = world
        self.max_solver_it = max_solver_it
        self.tolerance = tolerance
        self.w = w  # SOR relaxation factor

        # Set precision for decimal calculations
        decimal.getcontext().prec = 200  # Adjust this if necessary

        # Initialize reference values as Decimal
        self.phi0 = Decimal(0.0)
        self.Te0 = Decimal(0.0)
        self.n0 = Decimal(0.0)
        
        self.max_exp_value = Decimal(1000)  # This value can be adjusted based on your domain


    def set_reference_values(self, phi0, Te0, n0):
        self.phi0 = Decimal(phi0)
        self.Te0 = Decimal(Te0)
        self.n0 = Decimal(n0)
    
    def solve(self, project_dir, save_data=False):
        """
        Solves the potential using a relaxation method and optionally saves phi and L2 values.
        
        Args:
            project_dir (str): The directory where data will be saved.
            save_data (bool): If True, saves phi for every iteration and L2 values. Defaults to False.
        """
        if save_data:
            # Ensure the project directory exists
            if not os.path.exists(project_dir):
                os.makedirs(project_dir)
    
            # Create a subdirectory for phi iterations within the project directory
            phi_dir = os.path.join(project_dir, 'phi_iterations')
            if not os.path.exists(phi_dir):
                os.makedirs(phi_dir)
    
            # Open a file to store L2 values in the project directory
            l2_file_path = os.path.join(project_dir, "L2_values.txt")
            l2_file = open(l2_file_path, "w")
    
        # Convert numpy arrays to Decimal where necessary
        phi = np.array([[Decimal(p) for p in row] for row in self.world.phi.data], dtype=object)
        rho = np.array([[Decimal(r) for r in row] for row in self.world.rho.data], dtype=object)
    
        # Precompute 1/(dx^2) and 1/(dy^2) as Decimal
        dh = self.world.get_dh()
        idx2 = Decimal(1.0) / Decimal(dh[0]) ** 2
        idy2 = Decimal(1.0) / Decimal(dh[1]) ** 2
    
        L2 = Decimal(0)  # norm
        converged = False
    
        # Save the initial guess for phi before entering the loop
        if save_data:
            initial_phi_file_path = os.path.join(phi_dir, "phi_initial.txt")
            with open(initial_phi_file_path, "w") as phi_file:
                for row in phi:
                    phi_file.write(' '.join([str(p) for p in row]) + "\n")
        
        # Solve potential
        for it in range(self.max_solver_it):
            for i in range(self.world.ni):
                for j in range(self.world.nj):
                    # Skip Dirichlet boundaries (solid, fixed nodes)
                    if self.world.objectid.data[i, j] > 0:
                        continue
    
                    if i == 0:
                        phi[i, j] = phi[i + 1, j]
                    elif i == self.world.ni - 1:
                        phi[i, j] = phi[i - 1, j]
                    elif j == 0:
                        phi[i, j] = phi[i, j + 1]
                    elif j == self.world.nj - 1:
                        phi[i, j] = phi[i, j - 1]
                    else:
                        # Convert phi[i, j] to Decimal for computation
                        phi_val = Decimal(phi[i, j])
    
                        # Evaluate electron density using the Boltzmann relation
                        ne = self.n0 * Decimal.exp((phi_val - self.phi0) / self.Te0)
    
                        phi_new = ((rho[i, j] - Decimal(Const.QE) * ne) / Decimal(Const.EPS_0) +
                                   idx2 * (phi[i - 1, j] + phi[i + 1, j]) +
                                   idy2 * (phi[i, j - 1] + phi[i, j + 1])) / (2 * idx2 + 2 * idy2)
    
                        # Successive Over-Relaxation (SOR)
                        phi[i, j] = phi_val + Decimal(self.w) * (phi_new - phi_val)
    
            # Optionally save the current phi field to a file at each iteration
            if save_data:
                phi_file_path = os.path.join(phi_dir, f"phi_iter_{it}.txt")
                with open(phi_file_path, "w") as phi_file:
                    for row in phi:
                        phi_file.write(' '.join([str(p) for p in row]) + "\n")
    
            # Check for convergence every 25 iterations
            if it % 25 == 0:
                sum_sq = Decimal(0)
                for i in range(self.world.ni):
                    for j in range(self.world.nj):
                        if self.world.objectid.data[i, j] > 0:
                            continue
    
                        R = Decimal(0)
                        if i == 0:
                            R = phi[i, j] - phi[i + 1, j]
                        elif i == self.world.ni - 1:
                            R = phi[i, j] - phi[i - 1, j]
                        elif j == 0:
                            R = phi[i, j] - phi[i, j + 1]
                        elif j == self.world.nj - 1:
                            R = phi[i, j] - phi[i, j - 1]
                        else:
                            phi_val = Decimal(phi[i, j])
                            ne = self.n0 * Decimal.exp((phi_val - self.phi0) / self.Te0)
                            R = (-phi_val * (2 * idx2 + 2 * idy2) +
                                 (rho[i, j] - Decimal(Const.QE) * ne) / Decimal(Const.EPS_0) +
                                 idx2 * (phi[i - 1, j] + phi[i + 1, j]) +
                                 idy2 * (phi[i, j - 1] + phi[i, j + 1]))
    
                        sum_sq += R * R
    
                L2 = Decimal.sqrt(sum_sq / (self.world.ni * self.world.nj))
                print(L2)
    
                # Optionally write the L2 value to the L2 file
                if save_data:
                    l2_file.write(f"Iteration: {it}\tL2: {L2}\n")
    
                # Check for convergence
                if L2 < Decimal(self.tolerance):
                    if save_data:
                        l2_file.close()
                    return it
    
        if not converged:
            print(f"GS failed to converge, L2={L2}")
    
        if save_data:
            l2_file.close()
    
        return converged
    


    # Gauss-Seidel for grounded box example in 2D
    def solve_GS(self):
        phi = self.world.phi
        rho = self.world.rho

        # precompute 1/(dx^2) and 1/(dy^2)
        dh = self.world.get_dh()
        idx2 = 1.0 / (dh[0] * dh[0])
        idy2 = 1.0 / (dh[1] * dh[1])

        L2 = 0
        converged = False

        # solve potential
        for it in range(self.max_solver_it):
            for i in range(1, self.world.ni - 1):
                for j in range(1, self.world.nj - 1):
                    phi_new = (rho[i, j] / Const.EPS_0 +
                               idx2 * (phi[i - 1, j] + phi[i + 1, j]) +
                               idy2 * (phi[i, j - 1] + phi[i, j + 1])) / (2 * idx2 + 2 * idy2)

                    # SOR
                    phi.w_at(i, j, phi[i, j] + self.w * (phi_new - phi[i, j]))

            # check for convergence every 25 iterations
            if it % 50 == 0 and it != 0:
                sum_sq = 0
                for i in range(1, self.world.ni - 1):
                    for j in range(1, self.world.nj - 1):
                        R = (-phi[i, j] * (2 * idx2 + 2 * idy2) +
                             rho[i, j] / Const.EPS_0 +
                             idx2 * (phi[i - 1, j] + phi[i + 1, j]) +
                             idy2 * (phi[i, j - 1] + phi[i, j + 1]))

                        sum_sq += R * R

                L2 = sqrt(sum_sq / (self.world.ni * self.world.nj))
                if L2 < self.tolerance:
                    converged = True
                    break

        if not converged:
            print(f"GS failed to converge, L2={L2}")
        return converged

    def compute_ef(self):
        phi = self.world.phi  # reference to world.phi
        dh = self.world.get_dh()  # get cell spacing
        dx, dy = dh[0], dh[1]

        for i in range(self.world.ni):
            for j in range(self.world.nj):
                ef = self.world.ef.data[i, j]  # reference to (i,j) ef vector

                # x component, efx
                if i == 0:  # forward difference
                    ef[0] = (3 * phi(i, j) - 4 * phi(i + 1, j) + phi(i + 2, j)) / (2 * dx)
                elif i == self.world.ni - 1:  # backward difference
                    ef[0] = (4 * phi(i - 1, j) - phi(i - 2, j) - 3 * phi(i, j)) / (2 * dx)
                else:  # central difference
                    ef[0] = (phi(i - 1, j) - phi(i + 1, j)) / (2 * dx)

                # y component, efy
                if j == 0:  # forward difference
                    ef[1] = (3 * phi(i, j) - 4 * phi(i, j + 1) + phi(i, j + 2)) / (2 * dy)
                elif j == self.world.nj - 1:  # backward difference
                    ef[1] = (4 * phi(i, j - 1) - phi(i, j - 2) - 3 * phi(i, j)) / (2 * dy)
                else:  # central difference
                    ef[1] = (phi(i, j - 1) - phi(i, j + 1)) / (2 * dy)

        # print("Finished computing electric field.")
