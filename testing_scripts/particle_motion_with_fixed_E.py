# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 10:42:45 2024

@author: marce
"""
import pandas as pd
import matplotlib.pyplot as plt


# Constants for the simulation
ME = 9.10938356e-31  # Electron mass (kg)
QE = 1.602176634e-19  # Elementary charge (C)

# Load particle
m = ME  # Particle mass
q = -QE  # Particle charge
x = 0.0  # Initial position
v = 0.0  # Stationary initial velocity

# Simulation parameters
dt = 1e-9  # Time step
E = -100.0  # Electric field (V/m)

# Create and open particle trace file
out = open("trace3.csv", "w")
if not out:
    print("Failed to create file")
    exit(-1)

# Write the header
out.write("t,x,v\n")

# Particle loop
for it in range(20):
    # Write trace data
    out.write(f"{it * dt},{x},{v}\n")

    # Advance velocity and position
    x += v * dt
    v += (q / m) * E * dt

# Close the file
out.close()

# Normal exit
print("Simulation complete, data saved to trace.csv")


# Load the data from the CSV files
data1 = pd.read_csv('trace1.csv')

# Plot position (x) vs velocity (v) for all three traces
plt.figure(figsize=(8, 6))
plt.plot(data1['x'], data1['v'], marker='o', linestyle='-', label='Trace 1')


# Labeling the plot
plt.xlabel('Position (m)')
plt.ylabel('Velocity (m/s)')
plt.title('Position vs Velocity for Multiple Traces')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()

