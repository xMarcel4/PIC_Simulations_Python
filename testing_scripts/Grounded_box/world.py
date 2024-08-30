# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:30:51 2024

@author: marce
"""

import numpy as np
from fields import *
import numpy as np
import time
from vec3 import *
from numba import jit

class Const:
    EPS_0 = 8.8541878e-12   # C/V/m, vacuum permittivity
    QE = 1.602176565e-19     # C, electron charge
    AMU = 1.660838921e-27    # kg, atomic mass unit
    ME = 9.10938215e-31      # kg, electron mass
    K = 1.380648e-23         # J/K, Boltzmann constant
    PI = 3.141592653         # pi
    EvToK = QE / K           # 1 eV in K ~ 11604


class World:
    def __init__(self, ni, nj, nk):
        self.ni, self.nj, self.nk = ni, nj, nk
        self.phi = Field(ni, nj, nk)       # potential (scalar field)
        self.rho = Field(ni, nj, nk)       # charge density (scalar field)
        self.ef = Field(ni, nj, nk, 3)     # electric field components (vector field)
        self.node_vol = Field(ni, nj, nk)  # node volumes (scalar field)

        # Initialize mesh-related attributes
        self.x0 = np.zeros(3)  # mesh origin
        self.dh = np.zeros(3)  # cell spacing
        self.xm = np.zeros(3)  # mesh max bound
        self.xc = np.zeros(3)  # domain centroid

        # # # Initialize time-related attributes
        # self.dt = 2e-10        # size of time step
        # self.num_ts = 1        # number of time steps
        # self.ts = -1           # current time step
        # self.time = 0          # current simulation time
        # self.time_start = time.time()  # store the start time

    def set_extents(self, x0, xm):
        self.x0 = np.array(x0)
        self.xm = np.array(xm)

        # Compute cell spacing
        self.dh = (self.xm - self.x0) / (np.array([self.ni, self.nj, self.nk]) - 1.0)

        # Compute centroid
        self.xc = 0.5 * (self.x0 + self.xm)

        self.compute_node_volumes()
    
        # @nijt
    def x_to_l(self, x):
        """Convert a physical position x to logical coordinates."""
        lc = Vec3(
            (x[0] - self.x0[0]) / self.dh[0],
            (x[1] - self.x0[1]) / self.dh[1],
            (x[2] - self.x0[2]) / self.dh[2]
        )
        return lc


    
    
    def compute_node_volumes(self):
        for i in range(self.ni):    # loop over nodes
            for j in range(self.nj):
                for k in range(self.nk):
                    V = self.dh[0] * self.dh[1] * self.dh[2]  # standard volume
                    if i == 0 or i == self.ni - 1: V *= 0.5    # adjust on boundaries
                    if j == 0 or j == self.nj - 1: V *= 0.5
                    if k == 0 or k == self.nk - 1: V *= 0.5
                    self.node_vol.w_at(i, j, k, V)

    # @njit
    def compute_charge_density(self, species_list):
        self.rho = Field(self.ni, self.nj, self.nk)  # Reset charge density field to zero
        # print("Compute Charge Density:")
        
        for sp in species_list:  # Loop over species
            # Print the species name and charge
            # print(f"Processing species: {sp.name} with charge: {sp.charge}")
            
            if sp.charge == 0:
                continue  # Skip neutrals
            
            
            # Add the scaled density to rho
            self.rho += sp.charge * sp.den  # This should now work correctly
            
            # Print the updated `rho` value at the specific node (10, 10, 10)
            # print(f"rho at (10, 10, 10) after processing {sp.name}: {self.rho[10, 10, 10]}")
        
        # Optionally, print the final `rho` value at the specific node (10, 10, 10)
        # print(f"Final rho at (10, 10, 10): {self.rho[10, 10, 10]}")



    # @nijt
    def get_pe(self):
        pe = 0.0
        for i in range(self.ni):
            for j in range(self.nj):
                for k in range(self.nk):
                    # ef is a vector, so square each component and sum them up to get the magnitude squared
                    ef_node = self.ef[i, j, k]
                    ef2 = np.dot(ef_node, ef_node)  # This gives a scalar (magnitude squared)
                    pe += ef2 * self.node_vol[i, j, k]
                    
                    # Print the values at the specific node (10, 10, 10)
                    if i == 4  and j == 4 and k == 4:
                        print(f"ef_node at (10, 10, 10): {ef_node}")
                        print(f"ef2 at (10, 10, 10): {ef2}")
                        print(f"node_vol at (10, 10, 10): {self.node_vol[i, j, k]}")
                        print(f"Partial PE at (10, 10, 10): {ef2 * self.node_vol[i, j, k]}")
    
        return 0.5 * Const.EPS_0 * pe



    def get_x0(self):
        return self.x0

    def get_xm(self):
        return self.xm

    def get_xc(self):
        return self.xc

    def get_dh(self):
        return self.dh

    def get_nn(self):
        return self.ni, self.nj, self.nk
    
    def in_bounds(self, pos):
        """Check if the position is within the domain bounds."""
        for i in range(3):
            if pos[i] < self.x0[i] or pos[i] > self.xm[i]:
                return False
        return True

    def set_time(self, dt, num_ts):
        # Set the time step and number of time steps
        self.dt = dt
        self.num_ts = num_ts
        
        # Initialize the current time step and simulation time
        self.ts = 0  # Start at the first time step
        self.time = 0  # Reset simulation time to 0
        
        # Store the start time
        self.time_start = time.time()  # Record the start time of the simulation
    

    def get_ts(self):
        return self.ts

    def get_time(self):
        return self.time

    def get_dt(self):
        return self.dt

    def is_last_time_step(self):
        return self.ts == self.num_ts - 1

    def advance_time(self):
        """Advances to the next time step."""
        self.time += self.dt
        self.ts += 1
        return self.ts <= self.num_ts - 1

    def get_wall_time(self):
        """Returns the elapsed time since the simulation started."""
        return time.time() - self.time_start

