# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:30:51 2024

@author: marce
"""

import numpy as np
from fields_2 import Field
from fields_2 import *

class World:
    def __init__(self, ni, nj, nk):
        self.ni, self.nj,self.nk = nj,ni,nk
        self.nn = (ni, nj, nk)

        # Initialize fields as numpy arrays
        self.phi = Field(ni, nj, nk)       # potential
        self.rho = Field(ni, nj, nk)       # charge density
        self.ef = Field(ni, nj, nk)        # electric field components
        self.node_vol = Field(ni, nj, nk)  # node volumes

        print("World Constructor!")

    def set_extents(self, _x0, _xm):
        self.x0 = np.array(_x0)
        self.xm = np.array(_xm)

        # Compute cell spacing
        self.dh = (self.xm - self.x0) / (np.array(self.nn) - 1.0)

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
                    self.node_vol[i, j, k] = V

    def XtoL(self, x):
        x = np.array(x)
        lc = (x - self.x0) / self.dh
        return lc

    # Getter methods
    def get_x0(self):
        return self.x0

    def get_xm(self):
        return self.xm

    def get_xc(self):
        return self.xc

    def get_dh(self):
        return self.dh
