# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:57:53 2024

@author: marce
"""

import numpy as np

class ColdBeamSource2D:
    def __init__(self, species, world, v_drift, den):
        self.sp = species       # Reference to the injected species
        self.world = world      # Reference to world (2D)
        self.v_drift = v_drift  # Mean drift velocity in the Y direction
        self.den = den          # Injection density

    def sample(self):
        dh = self.world.get_dh()  # Get cell spacing
        x0 = self.world.get_x0()  # Get origin position

        # Length in X dimension
        Lx = dh[0] * (self.world.ni - 1)
        # In 2D, the area is simply the length in the X dimension (since we're assuming drift along Y)
        A = Lx

        # Compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
        num_real = self.den * self.v_drift * A * self.world.get_dt()

        # Number of simulation particles
        num_sim = int(num_real / self.sp.mpw0 + np.random.rand())

        # Inject particles
        for _ in range(num_sim):
            # Random X position, Y is set to the starting position at the bottom boundary (injection)
            pos = np.array([x0[0] + np.random.rand() * Lx, x0[1]])
            # Velocity: no motion in X, drift velocity in Y
            vel = np.array([0, self.v_drift])
            # Add particle to the species
            self.sp.add_particle(pos, vel, self.sp.mpw0)
