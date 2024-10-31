# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 10:42:45 2024

@author: marce
"""
import pandas as pd
import matplotlib.pyplot as plt

print("\033[H\033[J")

# Constants for the simulation
ME = 9.10938356e-31  # Electron mass (kg)
QE = 1.602176634e-19  # Elementary charge (C)

# Load particle
m = ME  # Particle mass
q = -QE  # Particle charge
x_euler = 0.0  # Initial position for Euler method
v_euler = 0.0  # Initial velocity for Euler method
x_backward = 0.0  # Initial position for Backward Euler method
v_backward = 0.0  # Initial velocity for Backward Euler method
x_leapfrog = 0.0  # Initial position for Leap-Frog method
v_leapfrog = 0.0  # Initial velocity for Leap-Frog method

# Simulation parameters
dt = 1e-9  # Time step
E = -100.0  # Electric field (V/m)
iterations = 10  # Number of iterations

# Create and open particle trace file for Explicit Euler method
with open("trace_euler.csv", "w") as out_euler:
    if not out_euler:
        print("Failed to create file")
        exit(-1)

    # Write the header
    out_euler.write("t,x,v\n")

    # Particle loop using Explicit Euler method
    for it in range(iterations):
        # Write trace data for Euler method
        out_euler.write(f"{it * dt},{x_euler},{v_euler}\n")

        # Advance velocity and position using Explicit Euler method
        x_euler += v_euler * dt  # x^(k+1) = x^k + v^k * dt
        v_euler += (q / m) * E * dt  # v^(k+1) = v^k + (q/m) * E(x^k) * dt

# Create and open particle trace file for Implicit Euler method
with open("trace_backward.csv", "w") as out_backward:
    if not out_backward:
        print("Failed to create file")
        exit(-1)

    # Write the header
    out_backward.write("t,x,v\n")

    # Particle loop using Implicit Euler method
    for it in range(iterations):
        # Advance velocity using Implicit Euler method
        v_backward += (q / m) * E * dt  # v^(k+1) = v^k + (q/m) * E(x^k) * dt

        # Advance position using Implicit Euler method
        x_backward += v_backward * dt  # x^(k+1) = x^k + v^(k+1) * dt

        # Write trace data for Implicit Euler method
        out_backward.write(f"{it * dt},{x_backward},{v_backward}\n")

# Create and open particle trace file for Leap-Frog method
with open("trace_leapfrog.csv", "w") as out_leapfrog:
    if not out_leapfrog:
        print("Failed to create file")
        exit(-1)

    # Write the header
    out_leapfrog.write("t,x,v\n")

    # Velocity rewind for Leap-Frog method (half-step back)
    v_leapfrog -= 0.5 * (q / m) * E * dt

    # Particle loop using Leap-Frog method
    for it in range(iterations):
        # Write trace data for Leap-Frog method
        out_leapfrog.write(f"{it * dt},{x_leapfrog},{v_leapfrog}\n")

        # Compute future velocity and position using Leap-Frog method
        v_leapfrog += (q / m) * E * dt
        x_leapfrog += v_leapfrog * dt

# Load the data from the CSV files
data_euler = pd.read_csv('trace_euler.csv')
data_backward = pd.read_csv('trace_backward.csv')
data_leapfrog = pd.read_csv('trace_leapfrog.csv')

# Plot position (x) vs velocity (v) for all three methods
plt.figure(figsize=(8, 6))
plt.plot(data_euler['x'], data_euler['v'], marker='o', linestyle='-', label='Explicit Euler Method')
plt.plot(data_backward['x'], data_backward['v'], marker='^', linestyle='-', label='Implicit Euler Method')
plt.plot(data_leapfrog['x'], data_leapfrog['v'], marker='s', linestyle='-', label='Leap-Frog Method')

# Labeling the plot
plt.xlabel('Position (m)')
plt.ylabel('Velocity (m/s)')
plt.title('Position vs Velocity: Euler vs Implicit Euler vs Leap-Frog Method')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
