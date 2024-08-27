# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 09:31:47 2024

@author: marce
"""

import numpy as np
from world import *
from fields import *
from vec3 import *
import random

class Particle:
    def __init__(self, pos, vel, mpw):
        self.pos = Vec3(*pos)  # position as Vec3
        self.vel = Vec3(*vel)  # velocity as Vec3
        self.mpw = mpw         # macroparticle weight


class Species:
    def __init__(self, name, mass, charge, world):
        self.name = name          # species name
        self.mass = mass          # particle mass in kg
        self.charge = charge      # particle charge in C
        self.world = world        # reference to the World object
        self.den = Field(world.ni, world.nj, world.nk)  # number density
        self.particles = []       # list for storing particles

    def get_np(self):
        return len(self.particles)  # simulation particles

    def get_real_count(self):
        mpw_sum = sum(part.mpw for part in self.particles)
        return mpw_sum

    def get_momentum(self):
        mom = Vec3(0, 0, 0)
        for part in self.particles:
            mom += part.mpw * part.vel
        return self.mass * mom

    def get_ke(self):
        ke = 0
        for part in self.particles:
            v2 = part.vel[0] ** 2 + part.vel[1] ** 2 + part.vel[2] ** 2
            ke += part.mpw * v2
        return 0.5 * self.mass * ke

    def advance(self):
        dt = self.world.get_dt()
        x0 = self.world.get_x0()
        xm = self.world.get_xm()

        for part in self.particles:
            lc = self.world.x_to_l(part.pos)
            ef_part = self.world.ef.gather(lc)

            part.vel += ef_part * (dt * self.charge / self.mass)
            part.pos += part.vel * dt

            for i in range(3):
                if part.pos[i] < x0[i]:
                    part.pos[i] = 2 * x0[i] - part.pos[i]
                    part.vel[i] *= -1.0
                elif part.pos[i] >= xm[i]:
                    part.pos[i] = 2 * xm[i] - part.pos[i]
                    part.vel[i] *= -1.0

    def compute_number_density(self):
        print(f"Computing number density for {self.name}...")
    
        # Reset density to 0
        self.den = Field(self.world.ni, self.world.nj, self.world.nk)
    
        # Loop over particles and compute the number density
        for idx, part in enumerate(self.particles):
            lc = self.world.x_to_l(part.pos)
            self.den.scatter(lc, part.mpw)
    
            # Print progress every 1000 particles (adjust as needed)
            if (idx + 1) % 100 == 0 or (idx + 1) == len(self.particles):
                print(f"Processed {idx + 1}/{len(self.particles)} particles")
    
        # Divide by node volume to get the final density
        self.den /= self.world.node_vol
    
        print(f"Finished computing number density for {self.name}")



    def add_particle(self, pos, vel, mpw):
        if not self.world.in_bounds(pos):
            return
        lc = self.world.x_to_l(pos)
        ef_part = self.world.ef.gather(lc)
        vel -= self.charge / self.mass * ef_part * (0.5 * self.world.get_dt())
        self.particles.append(Particle(pos, vel, mpw))

    def load_particles_box(self, x1, x2, num_den, num_sim):
        box_vol = np.prod([x2[i] - x1[i] for i in range(3)])
        num_real = num_den * box_vol
        mpw = num_real / num_sim

        for idx in range(num_sim):
            pos = [x1[i] + random.random() * (x2[i] - x1[i]) for i in range(3)]
            # Initial velocity is 0
            vel = Vec3(0, 0, 0)
            self.add_particle(pos, vel, mpw)
    
            # Print progress every 1000 particles
            if (idx + 1) % 100000 == 0 or (idx + 1) == num_sim:
                print(f"Loaded {idx + 1}/{num_sim} particles")

    def load_particles_box_qs(self, x1, x2, num_den, num_sim):
        box_vol = np.prod([x2[i] - x1[i] for i in range(3)])
        num_sim_tot = np.prod([num_sim[i] - 1 for i in range(3)])
        num_real = num_den * box_vol
        mpw = num_real / num_sim_tot

        di, dj, dk = [(x2[i] - x1[i]) / (num_sim[i] - 1) for i in range(3)]

        for i in range(num_sim[0]):
            for j in range(num_sim[1]):
                for k in range(num_sim[2]):
                    pos = [
                        x1[0] + i * di,
                        x1[1] + j * dj,
                        x1[2] + k * dk
                    ]

                    if pos[0] == x2[0]: pos[0] -= 1e-4 * di
                    if pos[1] == x2[1]: pos[1] -= 1e-4 * dj
                    if pos[2] == x2[2]: pos[2] -= 1e-4 * dk

                    w = 1
                    if i == 0 or i == num_sim[0] - 1: w *= 0.5
                    if j == 0 or j == num_sim[1] - 1: w *= 0.5
                    if k == 0 or k == num_sim[2] - 1: w *= 0.5

                    vel = Vec3(0, 0, 0)
                    self.add_particle(pos, vel, mpw * w)
