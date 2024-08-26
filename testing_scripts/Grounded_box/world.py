# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:30:51 2024

@author: marce
"""

import numpy as np
from fields_2 import Field
from fields_2 import *

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

    def set_extents(self, x0, xm):
        self.x0 = np.array(x0)
        self.xm = np.array(xm)

        # Compute cell spacing
        self.dh = (self.xm - self.x0) / (np.array([self.ni, self.nj, self.nk]) - 1.0)

        # Compute centroid
        self.xc = 0.5 * (self.x0 + self.xm)

        self.compute_node_volumes()

    def compute_node_volumes(self):
        for i in range(self.ni):    # loop over nodes
            for j in range(self.nj):
                for k in range(self.nk):
                    V = self.dh[0] * self.dh[1] * self.dh[2]  # standard volume
                    if i == 0 or i == self.ni - 1: V *= 0.5    # adjust on boundaries
                    if j == 0 or j == self.nj - 1: V *= 0.5
                    if k == 0 or k == self.nk - 1: V *= 0.5
                    self.node_vol.w_at(i, j, k, V)

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
