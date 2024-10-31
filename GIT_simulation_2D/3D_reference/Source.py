# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:33:52 2024

@author: marce
"""

import numpy as np

class ColdBeamSource:
    def __init__(self, species, world, v_drift, den):
        self.sp = species       # Reference to the injected species
        self.world = world      # Reference to world
        self.v_drift = v_drift  # Mean drift velocity
        self.den = den          # Injection density

    def sample(self):
        dh = self.world.get_dh()
        x0 = self.world.get_x0()

        # Area of the XY plane, A=Lx*Ly
        Lx = dh[0] * (self.world.ni - 1)
        Ly = dh[1] * (self.world.nj - 1)
        A = Lx * Ly

        # Compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
        num_real = self.den * self.v_drift * A * self.world.get_dt()
        # print(f"Num Real:   {num_real}")

        # Number of simulation particles
        num_sim = int(num_real / self.sp.mpw0 + np.random.rand())
        # print(f"Num Sim:   {num_sim}")

        # Inject particles
        for _ in range(num_sim):
            pos = np.array([x0[0] + np.random.rand() * Lx, 
                            x0[1] + np.random.rand() * Ly, 
                            x0[2]])
            vel = np.array([0, 0, self.v_drift])
            self.sp.add_particle(pos, vel, self.sp.mpw0)
