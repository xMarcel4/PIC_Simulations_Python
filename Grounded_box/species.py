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
# from numba import jit

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
    
    # @nijt
    def advance(self):
        dt = self.world.get_dt()
        x0 = self.world.get_x0()
        xm = self.world.get_xm()
        # print(f"dt: {dt}")
        # print(f"x0: {x0}")
        # print(f"xm: {xm}")  

        # Set the particle index to debug
        debug_particle_index = 1000
    
        # Iterate over all particles
        for p_idx, part in enumerate(self.particles):
            
            # Check if this is the particle you want to debug
            # if p_idx == debug_particle_index:
                # print(f"Debugging Particle {p_idx}:")
                # print(f"Initial position: ({part.pos[0]}, {part.pos[1]}, {part.pos[2]})")
                # print(f"Initial velocity: ({part.vel[0]}, {part.vel[1]}, {part.vel[2]})")
            
            # Get logical coordinate of particle's position
            lc = self.world.x_to_l(part.pos)
    
            # Get electric field at particle position
            ef_part = self.world.ef.gather(lc)
    
            # if p_idx == debug_particle_index:
            #     print(f"Logical coordinates: ({lc[0]}, {lc[1]}, {lc[2]})")
            #     print(f"Electric field at position: ({ef_part[0]}, {ef_part[1]}, {ef_part[2]})")
            #     print(f"Charge: {self.charge} \t mass: {self.mass}")

            
            # Update velocity from F = qE = ma
            part.vel += ef_part * (dt * self.charge / self.mass)
    
            # if p_idx == debug_particle_index:
            #     print(f"Updated velocity: ({part.vel[0]}, {part.vel[1]}, {part.vel[2]})")
            
            # Update position from v = dx/dt
            # part.pos += part.vel * dt
            part.pos[0] = part.pos[0] + part.vel[0] *dt
            part.pos[1] = part.pos[1] + part.vel[1] *dt
            part.pos[2] = part.pos[2] + part.vel[2] *dt

    
            # if p_idx == debug_particle_index:
                # print(f"Updated position: ({part.pos[0]}, {part.pos[1]}, {part.pos[2]})")
            
            # Reflect particles leaving the domain
            for i in range(3):
                if part.pos[i] < x0[i]:
                    # if p_idx == debug_particle_index:
                    #     print(f"Particle out of lower bound in dimension {i}:")
                    #     print(f"Original Position: {part.pos[i]}, Velocity: {part.vel[i]}")
                    
                    # Reflecting the particle
                    part.pos[i] = 2 * x0[i] - part.pos[i]
                    part.vel[i] *= -1.0
                    # if p_idx == debug_particle_index:
                    #     print(f"Reflected Position: {part.pos[i]}, Velocity: {part.vel[i]}")
                    
                elif part.pos[i] >= xm[i]:
                    # if p_idx == debug_particle_index:
    
                    #     print(f"Particle out of upper bound in dimension {i}:")
                    #     print(f"Original Position: {part.pos[i]}, Velocity: {part.vel[i]}")
                        
                    # Reflecting the particle
                    part.pos[i] = 2 * xm[i] - part.pos[i]
                    part.vel[i] *= -1.0

                    # if p_idx == debug_particle_index:
                    #     print(f"Reflected Position: {part.pos[i]}, Velocity: {part.vel[i]}")

    
            # if p_idx == debug_particle_index:
            #     print(f"Final position after reflection: ({part.pos[0]}, {part.pos[1]}, {part.pos[2]})")
            #     print(f"Final velocity after reflection: ({part.vel[0]}, {part.vel[1]}, {part.vel[2]})")
    

    
    # maybe add old compute number density with /= 
    
    # def compute_number_density(self):
    #     print("Compute Number Density:")
    #     # Initialize the density field to zero
    #     self.den = Field(self.world.ni, self.world.nj, self.world.nk)  # Reset density to zero
    #     print("Initialized density field to zero.")
    
    #     # Scatter particle weights into the density field
    #     for part in self.particles:
    #         lc = self.world.x_to_l(part.pos)
    #         self.den.scatter(lc, part.mpw)  # Assuming scatter is defined
    #     print("Scattered particle weights into the density field.")
    
    #     # Print the density at the specific node (10, 10, 10) after scattering
    #     print(f"Density at node (10, 10, 10) after scattering: {self.den.data[10, 10, 10]}")
    #     print(f"Node volume at node (10, 10, 10): {self.world.node_vol.data[10, 10, 10]}")

    
    #     # Perform element-wise division by the node volumes using __itruediv__
    #     try:
    #         print("Performing element-wise division by the node volumes.")
    #         self.den /= self.world.node_vol
    #     except Exception as e:
    #         print(f"Error during division: {e}")
    #         raise
    
    #     # Print the density at the specific node (10, 10, 10) after division
    #     print(f"Density at node (10, 10, 10) after division by node volumes: {self.den.data[10, 10, 10]}")
    
    #     # print("Finished computing number density.")

         

       # Perform element-wise division by the node volumes directly within the function

    # @nijt
    def compute_number_density(self):
        # print("Compute Number Density:")

        # Initialize the density field to zero
        self.den = Field(self.world.ni, self.world.nj, self.world.nk)  # Reset density to zero
        # print("Initialized density field to zero.")

            # Set the index of the particle you want to debug
        debug_particle_index = 1000
        current_particle_index = 0
    
        # Scatter particle weights into the density field
        for part in self.particles:
            lc = self.world.x_to_l(part.pos)
    
            # # Check if this is the particle you want to debug
            # if current_particle_index == debug_particle_index:
            #     print(f"Debugging Particle {current_particle_index}:")
            #     print(f"Physical position: ({part.pos[0]}, {part.pos[1]}, {part.pos[2]})")
            #     print(f"Logical coordinates: ({lc[0]}, {lc[1]}, {lc[2]})")
            #     print(f"Macroparticle weight (mpw): {part.mpw}")
    
            # Perform the scatter operation
            self.den.scatter(lc, part.mpw)  # Assuming scatter is defined
    
            current_particle_index += 1
    
        # # Print the density at the specific node (10, 10, 10) after scattering
        # print(f"Density at node (10, 10, 10) after scattering: {self.den.data[10, 10, 10]}")
        # print(f"Node volume at node (10, 10, 10): {self.world.node_vol.data[10, 10, 10]}")

        try:
            # Direct element-wise division without using __itruediv__
            for i in range(self.world.ni):
                for j in range(self.world.nj):
                    for k in range(self.world.nk):
                        if self.world.node_vol.data[i, j, k] != 0:
                            self.den.data[i, j, k] /= self.world.node_vol.data[i, j, k]
                        else:
                            self.den.data[i, j, k] = 0  # Handle division by zero
        except Exception as e:
            print(f"Error during manual division: {e}")
            raise
        # print(f"Density at node (10, 10, 10) after division by node volumes: {self.den.data[10, 10, 10]}")



    
    def add_particle(self, pos, vel, mpw):
        if not self.world.in_bounds(pos):
            print("Not in Bounds")
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
        num_sim_tot = np.prod([num_sim[i] - 1 for i in range(3)])  # Total number of simulation particles
        num_real = num_den * box_vol
        mpw = num_real / num_sim_tot  # Macroparticle weight
    
        # Compute particle grid spacing
        di, dj, dk = [(x2[i] - x1[i]) / (num_sim[i] - 1) for i in range(3)]
    
        particle_count = 0  # To keep track of particle loading progress
    
        for i in range(num_sim[0]):
            for j in range(num_sim[1]):
                for k in range(num_sim[2]):
                    pos = [
                        x1[0] + i * di,
                        x1[1] + j * dj,
                        x1[2] + k * dk
                    ]
    
                    # Adjust particles on max faces back into the domain
                    if pos[0] == x2[0]: pos[0] -= 1e-4 * di
                    if pos[1] == x2[1]: pos[1] -= 1e-4 * dj
                    if pos[2] == x2[2]: pos[2] -= 1e-4 * dk
    
                    w = 1  # Relative weight
                    if i == 0 or i == num_sim[0] - 1: w *= 0.5
                    if j == 0 or j == num_sim[1] - 1: w *= 0.5
                    if k == 0 or k == num_sim[2] - 1: w *= 0.5
    
                    vel = Vec3(0, 0, 0)  # Initialize velocity as stationary
                    self.add_particle(pos, vel, mpw * w)  # Add to the particles list
    
                    particle_count += 1
                    # if particle_count % 1000 == 0:
                    #     print(f"Loaded {particle_count}/{num_sim_tot} particles")
