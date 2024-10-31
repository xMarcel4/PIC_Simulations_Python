# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:05:41 2024

@author: marce
"""
import numpy as np
from fields_2d import *
import numpy as np
import time
from vec2 import *
import os

class Const:
    EPS_0 = 8.8541878e-12   # C/V/m, vacuum permittivity
    QE = 1.602176565e-19     # C, electron charge
    AMU = 1.660838921e-27    # kg, atomic mass unit
    ME = 9.10938215e-31      # kg, electron mass
    K = 1.380648e-23         # J/K, Boltzmann constant
    PI = 3.141592653         # pi
    EvToK = QE / K           # 1 eV in K ~ 11604


class World2D:
    def __init__(self, ni, nj):
        self.ni, self.nj = ni, nj
        self.phi = Field2D(ni, nj)       # potential (scalar field)
        self.rho = Field2D(ni, nj)       # charge density (scalar field)
        self.ef = Field2D(ni, nj, 2)     # electric field components (vector field, now 2D)
        self.node_vol = Field2D(ni, nj)  # node volumes (scalar field)

        # Initialize mesh-related attributes
        self.x0 = Vec2(0.0, 0.0)  # mesh origin in 2D (x, y)
        self.dh = Vec2(0.0, 0.0)  # cell spacing in 2D (dx, dy)
        self.xm = Vec2(0.0, 0.0)  # mesh max bound in 2D (x_max, y_max)
        self.xc = Vec2(0.0, 0.0)  # domain centroid in 2D (x_c, y_c)

        self.objectid = Field2D(ni, nj)  # Field to flag fixed nodes or objects
        
        # Initialize sphere-related attributes for 2D (now a circle)
        self.circle_x0 = Vec2(0.0, 0.0)  # Circle centroid
        self.circle_rad2 = 0             # Circle radius squared
        
        self.objects = []
        
        self.particles_out_of_bounds = 0
        self.particles_in_plane = 0
            
        
    def add_plane_with_hole(self, y_position, thickness, phi_plane, hole_half_width):
        """
        Adds a plane with a hole to the world. Each plane's properties are stored in a list.
        """
        # Store plane properties in a dictionary or list
        plane = {
            'y_position': y_position,
            'thickness': thickness,
            'phi_plane': phi_plane,
            'hole_half_width': hole_half_width
        }
        
        # Store all planes in a list
        if not hasattr(self, 'planes'):
            self.planes = []
        
        self.planes.append(plane)
        
        # Optional: Also update the grid object data here for the new plane
        half_thickness = thickness / 2
        x_center = (self.x0[0] + self.xm[0]) / 2
        
        for i in range(self.ni):
            for j in range(self.nj):
                x, y = self.pos(i, j)
                if y_position - half_thickness <= y <= y_position + half_thickness:
                    if not (x_center - hole_half_width <= x <= x_center + hole_half_width):
                        self.objectid.data[i, j] = 1
                        self.phi.data[i, j] = phi_plane
    

    
    def add_plane(self, y_position, thickness, phi_plane):
        """
        Adds a plane (horizontal line) to the 2D world at a specific Y position with a given thickness.
        
        Parameters:
        - y_position: The Y-coordinate where the center of the plane should be placed.
        - thickness: The thickness of the plane (distance in the Y-axis).
        - phi_plane: The potential to set within the plane.
        """
        half_thickness = thickness / 2  # Half of the thickness for easy range checking
    
        # Loop over all nodes in the 2D grid
        for i in range(self.ni):
            for j in range(self.nj):
                x, y = self.pos(i, j)  # Get the (x, y) coordinates of the current node
    
                # Check if the current node is within the plane's thickness range along the Y-axis
                if y_position - half_thickness <= y <= y_position + half_thickness:
                    self.objectid.data[i, j] = 1  # Set object flag for the plane
                    self.phi.data[i, j] = phi_plane  # Set potential inside the plane
    
                        
    def add_inlet(self):
        """
        Marks the bottom boundary (j=0) as a 0V Dirichlet boundary.
        """
        for i in range(self.ni):
            #self.objectid.data[i, 0] = 2  # Mark the bottom boundary
            self.phi.data[i, 0] = 0.0     # Set the potential to 0V
    

    def add_circle(self, center, radius, phi_circle):
        """
        Adds a circle-shaped object to the 2D world.
    
        Parameters:
        - center: Tuple with the (x, y) coordinates of the center of the circle.
        - radius: The radius of the circle.
        - phi_circle: The potential to set inside the circle.
        """
        radius2 = radius ** 2  # Square of the radius for distance comparison
    
        # Loop over all nodes in the 2D grid
        for i in range(self.ni):
            for j in range(self.nj):
                x, y = self.pos(i, j)  # Get the (x, y) coordinates of the current node
                
                # Calculate the distance from the center of the circle
                dist2 = (x - center[0]) ** 2 + (y - center[1]) ** 2
                
                # Check if the current node is inside the circle
                if dist2 <= radius2:
                    self.objectid.data[i, j] = 1  # Set object flag
                    self.phi.data[i, j] = phi_circle  # Set potential
                    

    def pos(self, i, j):
        """Convert logical coordinates (i, j) to physical position."""
        logical_coordinates = Vec2(i, j)
        x0_vec = Vec2(*self.x0)
        dh_vec = Vec2(*self.dh)
        physical_position = x0_vec + logical_coordinates * dh_vec
        return physical_position

    def in_circle(self, x):
        """Check if a point x is inside the circle (now 2D)."""
        return (x - self.circle_x0).dot(x - self.circle_x0) <= self.circle_rad2
    
    def in_plane_with_hole(self, x):
        """
        Check if a particle's position intersects any plane with a hole.
        
        Parameters:
        - x: The particle position (Vec3).
        
        Returns:
        - True if the particle intersects any plane and is outside the hole, False otherwise.
        """
        # Ensure we have planes to check
        if not hasattr(self, 'planes'):
            return False
        
        # Iterate over each plane
        for plane in self.planes:
            y_position = plane['y_position']
            thickness = plane['thickness']
            hole_half_width = plane['hole_half_width']
            x_center = (self.x0[0] + self.xm[0]) / 2
            
            half_thickness = thickness / 2
            
            # Check if particle is within the Y-range of this plane
            if y_position - half_thickness <= x[1] <= y_position + half_thickness:
                # Check if particle is outside the hole in the X-dimension
                if not (x_center - hole_half_width <= x[0] <= x_center + hole_half_width):
                    return True
                    print("In plane")
        
        # If the particle does not intersect any plane
        return False



    def set_extents(self, x0, xm):
        self.x0 = np.array(x0[:2])  # Keep only the 2D extents (x, y)
        self.xm = np.array(xm[:2])

        # Compute cell spacing in 2D
        self.dh = (self.xm - self.x0) / (np.array([self.ni, self.nj]) - 1.0)

        # Compute centroid in 2D
        self.xc = 0.5 * (self.x0 + self.xm)

        self.compute_node_volumes()

    def compute_node_volumes(self):
        for i in range(self.ni):    # loop over nodes
            for j in range(self.nj):
                V = self.dh[0] * self.dh[1]  # standard volume in 2D
                if i == 0 or i == self.ni - 1: V *= 0.5    # adjust on boundaries
                if j == 0 or j == self.nj - 1: V *= 0.5
                self.node_vol.w_at(i, j, V)
    
    def compute_charge_density(self, species_list):
            self.rho = Field2D(self.ni, self.nj)  # Reset charge density field to zero
            
            for sp in species_list:  # Loop over species
                if sp.charge == 0:
                    continue  # Skip neutrals
                self.rho += sp.charge * sp.den  # Add the scaled density to rho
    
    def get_pe(self):
        pe = 0.0
        for i in range(self.ni):
            for j in range(self.nj):
                ef_node = self.ef[i, j]
                ef2 = np.dot(ef_node, ef_node)  # This gives a scalar (magnitude squared)
                pe += ef2 * self.node_vol[i, j]
        return 0.5 * Const.EPS_0 * pe
    
    def x_to_l(self, x):
        """Convert a 2D physical position (x, y) to logical coordinates (i, j)."""
        lc = Vec2(
            (x[0] - self.x0[0]) / self.dh[0],  # Logical coordinate i
            (x[1] - self.x0[1]) / self.dh[1]   # Logical coordinate j
        )
        return lc


    def get_x0(self):
        return self.x0

    def get_xm(self):
        return self.xm

    def get_xc(self):
        return self.xc

    def get_dh(self):
        return self.dh

    def get_nn(self):
        return self.ni, self.nj

    def in_bounds(self, pos):
        """Check if the position is within the domain bounds."""
        for i in range(2):
            if pos[i] < self.x0[i] or pos[i] > self.xm[i]:
                return False
                print("Out of bounds")
        return True

    def set_time(self, dt, num_ts):
        self.dt = dt
        self.num_ts = num_ts
        self.ts = 0  # Start at the first time step
        self.time = 0  # Reset simulation time to 0
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
        self.time += self.dt
        self.ts += 1
        return self.ts <= self.num_ts - 1

    def get_wall_time(self):
        return time.time() - self.time_start
    
    def save_field(self, field_name, project_dir):
        field_data = getattr(self, field_name).data
        filename = os.path.join(project_dir, f'{field_name}.npy')
        np.save(filename, field_data)


    def load_field(self, field_name, project_dir):
        filename = os.path.join(project_dir, f'{field_name}.npy')
        field_data = np.load(filename)
        setattr(self, field_name, Field(self.ni, self.nj))  # Reinitialize the field object
        getattr(self, field_name).data = field_data

    def save_phi(self, project_dir):
        self.save_field('phi', project_dir)

    def load_phi(self, project_dir):
        self.load_field('phi', project_dir)

    def save_ef(self, project_dir):
        self.save_field('ef', project_dir)

    def load_ef(self, project_dir):
        self.load_field('ef', project_dir)
