# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:56:38 2024

@author: marce
"""

import numpy as np
from world_2d import *  # Assuming world_2d contains the 2D World class
from fields_2d import *  # Assuming fields_2d contains the 2D Field class
from vec2 import *  # Assuming vec2 is the 2D version of vec3
import random


class Particle2D:
    def __init__(self, pos, vel, mpw):
        self.pos = Vec2(*pos)  # position as Vec2
        self.vel = Vec2(*vel)  # velocity as Vec2
        self.mpw = mpw         # macroparticle weight
        


class Species2D:
    def __init__(self, name, mass, charge, world, mpw0):
        self.name = name          # species name
        self.mass = mass          # particle mass in kg
        self.charge = charge      # particle charge in C
        self.world = world        # reference to the 2D World object
        self.mpw0 = mpw0          # default macroparticle weight
        self.den = Field2D(world.ni, world.nj)  # number density in 2D
        self.particles = []       # list for storing particles
        


    def get_np(self):
        return len(self.particles)  # simulation particles

    def get_real_count(self):
        mpw_sum = sum(part.mpw for part in self.particles)
        return mpw_sum

    def get_momentum(self):
        mom = Vec2(0, 0)
        for part in self.particles:
            mom += part.mpw * part.vel
        return self.mass * mom

    def get_ke(self):
        ke = 0
        for part in self.particles:
            v2 = part.vel[0] ** 2 + part.vel[1] ** 2  # Only x and y components in 2D
            ke += part.mpw * v2
        return 0.5 * self.mass * ke
    
    # @nijt
    def advance(self):

        dt = self.world.get_dt()
        x0 = self.world.get_x0()
        xm = self.world.get_xm()

        for part in self.particles:
            # Get logical coordinate of particle's position
            lc = self.world.x_to_l(part.pos)

            # Get electric field at particle position
            ef_part = self.world.ef.gather(lc)

            # Update velocity from F = qE = ma
            part.vel += ef_part * (dt * self.charge / self.mass)

            # Update position from v = dx/dt
            part.pos[0] += part.vel[0] * dt
            part.pos[1] += part.vel[1] * dt

            # Check if particle is out of bounds
            if not self.world.in_bounds(part.pos):
                self.world.particles_out_of_bounds += 1
                part.mpw = 0  # Mark the particle as dead

            # Check if the particle intersects the plane
            if self.world.in_plane_with_hole(part.pos):
                self.world.particles_in_plane += 1
                part.mpw = 0  # Mark the particle as dead

        # Remove dead particles (with zero macroparticle weight)
        self.particles = [part for part in self.particles if part.mpw > 0]

        # Output the number of out-of-bounds particles
        #print(f"Particles out of bounds: {particles_out_of_bounds}")

    def compute_number_density(self):
        # Initialize the density field to zero
        self.den = Field2D(self.world.ni, self.world.nj)  # Reset density to zero

        # Scatter particle weights into the density field
        for part in self.particles:
            lc = self.world.x_to_l(part.pos)
            self.den.scatter(lc, part.mpw)

        # Perform element-wise division by the node volumes directly within the function
        for i in range(self.world.ni):
            for j in range(self.world.nj):
                if self.world.node_vol.data[i, j] != 0:
                    self.den.data[i, j] /= self.world.node_vol.data[i, j]
                else:
                    self.den.data[i, j] = 0  # Handle division by zero

    def add_particle(self, pos, vel, mpw):
        if not self.world.in_bounds(pos):
            print("Not in Bounds")
            return
        
        lc = self.world.x_to_l(pos)
        ef_part = self.world.ef.gather(lc)
        
        vel = np.array(vel, dtype=np.float64)  # Ensure vel is a numpy array with float type
        
        # Update the velocity
        vel -= self.charge / self.mass * ef_part * (0.5 * self.world.get_dt())
        
        # Append the particle to the list
        self.particles.append(Particle2D(pos, vel, mpw))
    
    def load_particles_box(self, x1, x2, num_den, num_sim):
        box_vol = np.prod([x2[i] - x1[i] for i in range(2)])  # 2D volume (area)
        num_real = num_den * box_vol
        mpw = num_real / num_sim

        for idx in range(num_sim):
            pos = [x1[i] + random.random() * (x2[i] - x1[i]) for i in range(2)]
            vel = Vec2(0, 0)  # Initial velocity is 0 in 2D
            self.add_particle(pos, vel, mpw)

            if (idx + 1) % 1000 == 0 or (idx + 1) == num_sim:
                print(f"Loaded {idx + 1}/{num_sim} particles")

    def load_particles_box_qs(self, x1, x2, num_den, num_sim):
        box_area = np.prod([x2[i] - x1[i] for i in range(2)])  # 2D area
        num_sim_tot = np.prod([num_sim[i] - 1 for i in range(2)])  # Total number of simulation particles
        num_real = num_den * box_area
        mpw = num_real / num_sim_tot  # Macroparticle weight
    
        di, dj = [(x2[i] - x1[i]) / (num_sim[i] - 1) for i in range(2)]  # Particle grid spacing

        particle_count = 0
        for i in range(num_sim[0]):
            for j in range(num_sim[1]):
                pos = [
                    x1[0] + i * di,
                    x1[1] + j * dj,
                ]

                # Adjust particles on max faces back into the domain
                if pos[0] == x2[0]: pos[0] -= 1e-4 * di
                if pos[1] == x2[1]: pos[1] -= 1e-4 * dj

                w = 1  # Relative weight
                if i == 0 or i == num_sim[0] - 1: w *= 0.5
                if j == 0 or j == num_sim[1] - 1: w *= 0.5

                vel = Vec2(0, 0)  # Initialize velocity as stationary
                self.add_particle(pos, vel, mpw * w)

                particle_count += 1
                if particle_count % 1000 == 0:
                    print(f"Loaded {particle_count}/{num_sim_tot} particles")
