# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:00:43 2024

@author: marce
"""
from fields import *
from world import *
import numpy as np
from math import sqrt,exp
import time
# from numba import njit

class PotentialSolver:
    def __init__(self, world, max_solver_it, tolerance, w=1.9):
        self.world = world
        self.max_solver_it = max_solver_it
        self.tolerance = tolerance
        self.w = w  # SOR relaxation factor
        
        # Initialize reference values with default placeholders
        self.phi0 = 0.0
        self.Te0 = 0.0
        self.n0 = 0.0

    def set_reference_values(self, phi0, Te0, n0):
        self.phi0 = phi0
        self.Te0 = Te0
        self.n0 = n0

    def solve(self):
        phi = self.world.phi.data
        rho = self.world.rho.data  # rho contains only ion contribution

        # precompute 1/(dx^2)
        dh = self.world.get_dh()
        idx2 = 1.0 / (dh[0] * dh[0])
        idy2 = 1.0 / (dh[1] * dh[1])
        idz2 = 1.0 / (dh[2] * dh[2])

        L2 = 0  # norm
        converged = False

        # solve potential
        for it in range(self.max_solver_it):
            for i in range(self.world.ni):
                for j in range(self.world.nj):
                    for k in range(self.world.nk):
                        # skip over solid (fixed) nodes = Dirichlet boundaries
                        if self.world.objectid.data[i, j, k] > 0:
                            continue

                        if i == 0:
                            phi[i, j, k] = phi[i + 1, j, k]
                        elif i == self.world.ni - 1:
                            phi[i, j, k] = phi[i - 1, j, k]
                        elif j == 0:
                            phi[i, j, k] = phi[i, j + 1, k]
                        elif j == self.world.nj - 1:
                            phi[i, j, k] = phi[i, j - 1, k]
                        elif k == 0:
                            phi[i, j, k] = phi[i, j, k + 1]
                        elif k == self.world.nk - 1:
                            phi[i, j, k] = phi[i, j, k - 1]
                        else:  # standard internal open node
                            # evaluate electron density from the Boltzmann relationship
                            ne = self.n0 * np.exp((phi[i, j, k] - self.phi0) / self.Te0)

                            phi_new = ((rho[i, j, k] - Const.QE * ne) / Const.EPS_0 +
                                       idx2 * (phi[i - 1, j, k] + phi[i + 1, j, k]) +
                                       idy2 * (phi[i, j - 1, k] + phi[i, j + 1, k]) +
                                       idz2 * (phi[i, j, k - 1] + phi[i, j, k + 1])) / (2 * idx2 + 2 * idy2 + 2 * idz2)

                            # SOR
                            phi[i, j, k] = phi[i, j, k] + self.w * (phi_new - phi[i, j, k])

            # check for convergence every 25 iterations
            if it % 25 == 0:
                sum_sq = 0
                for i in range(self.world.ni):
                    for j in range(self.world.nj):
                        for k in range(self.world.nk):
                            # skip over solid (fixed) nodes
                            if self.world.objectid.data[i, j, k] > 0:
                                continue

                            R = 0
                            if i == 0:
                                R = phi[i, j, k] - phi[i + 1, j, k]
                            elif i == self.world.ni - 1:
                                R = phi[i, j, k] - phi[i - 1, j, k]
                            elif j == 0:
                                R = phi[i, j, k] - phi[i, j + 1, k]
                            elif j == self.world.nj - 1:
                                R = phi[i, j, k] - phi[i, j - 1, k]
                            elif k == 0:
                                R = phi[i, j, k] - phi[i, j, k + 1]
                            elif k == self.world.nk - 1:
                                R = phi[i, j, k] - phi[i, j, k - 1]
                            else:
                                # evaluate electron density from the Boltzmann relationship
                                ne = self.n0 * np.exp((phi[i, j, k] - self.phi0) / self.Te0)
                                R = (-phi[i, j, k] * (2 * idx2 + 2 * idy2 + 2 * idz2) +
                                     (rho[i, j, k] - Const.QE * ne) / Const.EPS_0 +
                                     idx2 * (phi[i - 1, j, k] + phi[i + 1, j, k]) +
                                     idy2 * (phi[i, j - 1, k] + phi[i, j + 1, k]) +
                                     idz2 * (phi[i, j, k - 1] + phi[i, j, k + 1]))

                            sum_sq += R * R

                L2 = np.sqrt(sum_sq / (self.world.ni * self.world.nj * self.world.nk))

                # Output the iteration, L2 norm, and tolerance
                print(f"Iteration: {it}\t L2: {L2}\t Tolerance: {self.tolerance}")

                if L2 < self.tolerance:
                    converged = True
                    break

        if not converged:
            print(f"GS failed to converge, L2={L2}")

        return converged

    # simple gauss seidel form grounded box example
    # @njit
    def solve_GS(self):
        phi = self.world.phi  # reference to world.phi
        rho = self.world.rho  # reference to world.rho
    
        # precompute 1/(dx^2), 1/(dy^2), 1/(dz^2)
        dh = self.world.get_dh()
        idx2 = 1.0 / (dh[0] * dh[0])
        idy2 = 1.0 / (dh[1] * dh[1])
        idz2 = 1.0 / (dh[2] * dh[2])
        
                # Example variables (these would normally come from your code)
        # phi = self.world.phi
        # rho = self.world.rho
        # dh = self.world.get_dh()
        ni, nj, nk = self.world.ni, self.world.nj, self.world.nk
        
    

        # Print out the precomputed values
        # print(f"idx2: {idx2}, idy2: {idy2}, idz2: {idz2}")
        
        # Print out the value of rho at position (10, 10, 10)
        # print(f"rho(10, 10, 10): {rho[10, 10, 10]}")
    
        L2 = 0  # norm
        converged = False
        # solve potential
        for it in range(self.max_solver_it):
            #print(f"Iteration {it + 1}/{self.max_solver_it}")
            #time.sleep(0.01)
            for i in range(1, self.world.ni - 1):
                #(f"test {i} out of {self.world.ni - 1}")
                for j in range(1, self.world.nj - 1):
                    for k in range(1, self.world.nk - 1):
                        # standard internal open node
                        phi_new = (rho[i, j, k] / Const.EPS_0 +
                                   idx2 * (phi[i - 1, j, k] + phi[i + 1, j, k]) +
                                   idy2 * (phi[i, j - 1, k] + phi[i, j + 1, k]) +
                                   idz2 * (phi[i, j, k - 1] + phi[i, j, k + 1])) / \
                                  (2 * idx2 + 2 * idy2 + 2 * idz2)
    
                        # SOR
                        phi.w_at(i, j, k, phi[i, j, k] + self.w  * (phi_new - phi[i, j, k]))
          
            # check for convergence every 25 iterations
            if it % 50 == 0 and it != 0:
                #print('check for convergence')
                sum_sq = 0
                for i in range(1, self.world.ni - 1):
                    for j in range(1, self.world.nj - 1):
                        for k in range(1, self.world.nk - 1):
                            R = (-phi[i, j, k] * (2 * idx2 + 2 * idy2 + 2 * idz2) +
                                 rho[i, j, k] / Const.EPS_0 +
                                 idx2 * (phi[i - 1, j, k] + phi[i + 1, j, k]) +
                                 idy2 * (phi[i, j - 1, k] + phi[i, j + 1, k]) +
                                 idz2 * (phi[i, j, k - 1] + phi[i, j, k + 1]))
    
                            sum_sq += R * R
    
                L2 = sqrt(sum_sq / (self.world.ni * self.world.nj * self.world.nk))
                if L2 < self.tolerance:
                    converged = True
                    break
    
        if not converged:
            print(f"GS failed to converge, L2={L2}")
        return converged

        
    # @njit
    def compute_ef(self):
        phi = self.world.phi  # reference to world.phi
        dh = self.world.get_dh()  # get cell spacing
        dx, dy, dz = dh[0], dh[1], dh[2]
        # print(f"dx: {dx}, dy: {dy}, dz: {dz}")

        for i in range(self.world.ni):
            for j in range(self.world.nj):
                for k in range(self.world.nk):
                    ef = self.world.ef.data[i, j, k]  # reference to (i,j,k) ef vector
    
                    # Debug prints for grid spacing
    
                    # Print the value of phi at (10, 10, 10)
                    # if i == 10 and j == 10 and k == 10:
                        # print(f"phi(10, 10, 10): {phi(i, j, k)}")
    
                    # x component, efx
                    if i == 0:  # forward difference
                        ef[0] = (3 * phi(i, j, k) - 4 * phi(i + 1, j, k) + phi(i + 2, j, k)) / (2 * dx)
                    elif i == self.world.ni - 1:  # backward difference
                        ef[0] = (4 * phi(i - 1, j, k) - phi(i - 2, j, k) - 3 * phi(i, j, k)) / (2 * dx)
                    else:  # central difference
                        ef[0] = (phi(i - 1, j, k) - phi(i + 1, j, k)) / (2 * dx)
    
                    # y component, efy
                    if j == 0:  # forward difference
                        ef[1] = (3 * phi(i, j, k) - 4 * phi(i, j + 1, k) + phi(i, j + 2, k)) / (2 * dy)
                    elif j == self.world.nj - 1:  # backward difference
                        ef[1] = (4 * phi(i, j - 1, k) - phi(i, j - 2, k) - 3 * phi(i, j, k)) / (2 * dy)
                    else:  # central difference
                        ef[1] = (phi(i, j - 1, k) - phi(i, j + 1, k)) / (2 * dy)
    
                    # z component, efz
                    if k == 0:  # forward difference
                        ef[2] = (3 * phi(i, j, k) - 4 * phi(i, j, k + 1) + phi(i, j, k + 2)) / (2 * dz)
                    elif k == self.world.nk - 1:  # backward difference
                        ef[2] = (4 * phi(i, j, k - 1) - phi(i, j, k - 2) - 3 * phi(i, j, k)) / (2 * dz)
                    else:  # central difference
                        ef[2] = (phi(i, j, k - 1) - phi(i, j, k + 1)) / (2 * dz)
    
        # print("Finished computing electric field.")
