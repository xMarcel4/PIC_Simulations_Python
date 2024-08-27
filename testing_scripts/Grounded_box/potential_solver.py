# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:00:43 2024

@author: marce
"""
from fields import *
from world import *
import numpy as np
from math import sqrt
import time

class PotentialSolver:
    def __init__(self, world, max_it, tol):
        self.world = world
        self.max_solver_it = max_it
        self.tolerance = tol

    def solve(self):
        phi = self.world.phi  # reference to world.phi
        rho = self.world.rho  # reference to world.rho

        # precompute 1/(dx^2), 1/(dy^2), 1/(dz^2)
        dh = self.world.get_dh()
        idx2 = 1.0 / (dh[0] * dh[0])
        idy2 = 1.0 / (dh[1] * dh[1])
        idz2 = 1.0 / (dh[2] * dh[2])

        L2 = 0  # norm
        converged = False

        # solve potential
        for it in range(self.max_solver_it):
            #print(f"Iteration {it + 1}/{self.max_solver_it}")
            time.sleep(0.01)
            for i in range(1, self.world.ni - 1):
                #(f"test {i} out of {self.world.ni - 1}")
                for j in range(1, self.world.nj - 1):
                    for k in range(1, self.world.nk - 1):
                        # standard internal open node
                        phi_new = (rho(i, j, k) / Const.EPS_0 +
                                   idx2 * (phi(i - 1, j, k) + phi(i + 1, j, k)) +
                                   idy2 * (phi(i, j - 1, k) + phi(i, j + 1, k)) +
                                   idz2 * (phi(i, j, k - 1) + phi(i, j, k + 1))) / \
                                  (2 * idx2 + 2 * idy2 + 2 * idz2)

                        # SOR
                        phi.w_at(i, j, k, phi(i, j, k) + 1.4 * (phi_new - phi(i, j, k)))

            # check for convergence every 25 iterations
            if it % 25 == 0 and it != 0:
                #print('check for convergence')
                sum_sq = 0
                for i in range(1, self.world.ni - 1):
                    for j in range(1, self.world.nj - 1):
                        for k in range(1, self.world.nk - 1):
                            R = (-phi(i, j, k) * (2 * idx2 + 2 * idy2 + 2 * idz2) +
                                 rho(i, j, k) / Const.EPS_0 +
                                 idx2 * (phi(i - 1, j, k) + phi(i + 1, j, k)) +
                                 idy2 * (phi(i, j - 1, k) + phi(i, j + 1, k)) +
                                 idz2 * (phi(i, j, k - 1) + phi(i, j, k + 1)))

                            sum_sq += R * R

                L2 = sqrt(sum_sq / (self.world.ni * self.world.nj * self.world.nk))
                if L2 < self.tolerance:
                    converged = True
                    break

        if not converged:
            #print(f"GS failed to converge, L2={L2}")
            return converged

    def compute_ef(self):
        phi = self.world.phi  # reference to world.phi
        dh = self.world.get_dh()  # get cell spacing
        dx, dy, dz = dh[0], dh[1], dh[2]

        for i in range(self.world.ni):
            for j in range(self.world.nj):
                for k in range(self.world.nk):
                    ef = self.world.ef.data[i, j, k]  # reference to (i,j,k) ef vector

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
