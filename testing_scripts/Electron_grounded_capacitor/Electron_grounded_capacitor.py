# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 12:17:32 2024

@author: Marcel
"""

#Problemstellung:
"""
    - Kondensatorplatten unendlich groß -> symmetrie nur 1D
    - Plattenabstand 10cm
    - Direlect Randbedinung -> Potential an beiden Platten = 0
    - zwischen Platten homogene Ionendichte n = 1e12 1/m  -> Einfluss von Elektron vernachlässigbar
    - Unterteilung in 20 räumliche Domänen -> 21 Knotenpunkte
    - keine Reibung, kein B-Feld nur elektrostatischen Teil der Lorentzkraft E-Feld
    - keine Stöße von Elektronen mit Ionen
    - Elektron wird an beliebigen Punkt zwischen den Platten platziert

"""

# Felder neu berechnen hier nicht da IOnendichte konstant bleibt
print("\033[H\033[J")
import numpy as np
import csv
import matplotlib.pyplot as plt


EPS_0 = 8.85418782e-12 # vacuum permittivity
QE = 1.602176565e-19  # Elementary charge (Coulombs)
ME = 9.10938215e-31  # Electron mass (kg)
max_it = 5000

x_positions = []
ts_values = []
pe_array = []
ke_array = []

nodes = 21
dt = 1e-10         # time increment


# Fields: rho, phi, ef (electric field)
rho = np.full(nodes, QE * 1e12)  # Initialize with ion density * QE
phi = np.zeros(nodes)  # Initialize with zeros
ef = np.zeros(nodes)   # Initialize with zeros

# Initialization
def main():

    #Distance 0,1m
    x0 = 0.0  # Origin
    xd = 0.1  # Opposite end
    dx = (xd - x0) / (nodes - 1) # Calculate node spacing
    x_values = np.linspace(x0, xd, nodes)  # Example x-axis values from 0 to 0.1 m 

    
    # Define Electron
    m = ME  # Mass of the electron
    q = -QE  # Charge of the electron
    x = 4 * dx  # 4 cells from the left edge (at start)
    v = 0.0  # No initial velocity
    
    solvePotentialGS(dx, phi, rho, max_it)
    computeEF(dx, ef, phi,2)
    
    
    # Velocity rewind for Leapfrog method
    li = XtoL(x, dx)  # Convert physical position x to logical coordinate li
    ef_p = gather(li, ef)  # Get electric field at particle position
    v -=0.5 * (q / m) * ef_p * dt

    #NOT so important for simulation, just for demonstration purposes
    phi_max = max(phi)

    
    # Open the file for writing
    with open("trace.csv", "w", newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["time", "x", "v", "KE", "PE"])
        
    # Main loop
    x_old = x
    ts_ary = np.arange(1, 5001)

    for ts in ts_ary:  # Particle loop
    
        # Sample mesh data at particle position
        li = XtoL(x, dx)
        ef_p = gather(li, ef)
    
        # Integrate velocity and position
        x_old = x
        
        v += (q / m) * ef_p * dt
        
        x += v * dt
        x_positions.append(x)

        phi_p = gather(XtoL(0.5 * (x + x_old), dx), phi)  # Interpolated phi at x^(k-0.5)
        pe = q * (phi_p - phi_max) / QE  # Potential energy in eV
        pe_array.append(pe)
        ke = 0.5 * m * v**2 / QE  # Kinetic energy in eV
        ke_array.append(ke)

        # Write to file
        #writer.writerow([ts * dt, x, v, ke, pe])
    
        # Screen output
        if ts == 1 or ts % 1000 == 0:
            print(f"ts: {ts}, x: {x}, v: {v}, ke: {ke}, pe: {pe}")
    
    # Call the function to output the CSV
    outputCSV(x0, dx, phi, rho, ef)
    # Wait for user input before closing (simulating std::cin.get())
    
    plot_simulation_results(x_values, rho, phi, ef, ts_ary *1e-4, x_positions, pe_array, ke_array)
    print("finished :)")
    
    # Normal exit
    return rho,ef,x_positions,ts_values,phi  # This line is not necessary in Python if main is not defined as a function, but can be used if it is.
                      
def XtoL(x, dx,x_0=0):
    return (x - x_0) / dx

def solvePotentialGS(dx, phi, rho, max_it):
    L2 = 0
    dx2 = dx * dx  # Precompute dx^2
    w = 1.4  # Relaxation factor
    ni = len(phi)  # Number of mesh nodes

    # Solve potential
    phi[0] = 0  # Dirichlet boundary on the left
    phi[ni - 1] = 0  # Dirichlet boundary on the right

    for solver_it in range(max_it):
        # Gauss-Seidel method with Successive Over-Relaxation (SOR)
        for i in range(1, ni - 1):
            g = 0.5 * (phi[i - 1] + phi[i + 1] + dx2 * rho[i] / EPS_0)
            phi[i] = phi[i] + w * (g - phi[i])

        # Check for convergence every 50 iterations
        if solver_it % 50 == 0:
            sum_val = 0.0
            # Internal nodes, automatically satisfied on Dirichlet boundaries
            for i in range(1, ni - 1):
                R = -rho[i] / EPS_0 - (phi[i - 1] - 2 * phi[i] + phi[i + 1]) / dx2
                sum_val += R * R

            L2 = np.sqrt(sum_val / ni)
            if L2 < 1e-6:
                print(f"Gauss-Seidel converged after {solver_it} iterations")
                return True

    print(f"Gauss-Seidel failed to converge, L2 = {L2}")
    return False

def gather(li, field):
    i = int(li)  # Cast to integer: number of primary node
    di = li - i  # Obtain fractional part: distance to next node
    # Perform linear interpolation
    return field[i] * (1 - di) + field[i + 1] * di

def computeEF(dx, ef, phi, second_order):
    ni = len(phi)  # Number of mesh nodes

    # Central difference on internal nodes
    for i in range(1, ni - 1):
        ef[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx)

    # Boundaries
    if second_order:
        print("second order")
        # Second order difference at boundaries
        ef[0] = (3 * phi[0] - 4 * phi[1] + phi[2]) / (2 * dx)
        ef[ni - 1] = -(phi[ni - 3] - 4 * phi[ni - 2] + 3 * phi[ni - 1]) / (2 * dx)
    else:
        # First order difference at boundaries
        ef[0] = (phi[0] - phi[1]) / dx
        ef[ni - 1] = (phi[ni - 2] - phi[ni - 1]) / dx
 
def outputCSV(x0, dx, phi, rho, ef):
    try:
        with open("results.csv", "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write header
            writer.writerow(["x", "phi", "rho", "ef"])
            
            # Write data rows
            for i in range(len(phi)):
                x_i = x0 + i * dx  # Calculate the i-th position
                writer.writerow([x_i, phi[i], rho[i], ef[i]])
        
        return True
    
    except IOError:
        print("Could not open output file!")
        return False

def plot_simulation_results(x_values, rho, phi, ef, ts, x_positions, pe_array, ke_array):
    """
    Function to plot the results of the simulation.

    Parameters:
    - x_values: array-like, x-axis values (e.g., positions in meters)
    - rho: array-like, charge density (e.g., in C/m)
    - phi: array-like, potential (e.g., in volts)
    - ef: array-like, electric field (e.g., in V/m)
    - ts: array-like, time values (e.g., in microseconds)
    - x_positions: array-like, positions over time (e.g., in meters)
    - pe_array: array-like, potential energy over time (e.g., in eV)
    - ke_array: array-like, kinetic energy over time (e.g., in eV)
    """
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Plot rho
    axs[0, 0].plot(x_values, rho)
    axs[0, 0].set_title('rho [C/m]')
    axs[0, 0].set_xlabel('x [m]')

    # Plot phi
    axs[1, 0].plot(x_values, phi)
    axs[1, 0].set_title('phi [V]')
    axs[1, 0].set_xlabel('x [m]')

    # Plot ef
    axs[1, 1].plot(x_values, ef)
    axs[1, 1].set_title('ef [V/m]')
    axs[1, 1].set_xlabel('x [m]')

    # Plot x positions over time
    axs[0, 2].plot(ts, x_positions)
    axs[0, 2].set_title('x [m]')
    axs[0, 2].set_xlabel('t [us]')

    pe_array = np.array(pe_array)
    ke_array = np.array(ke_array)
    total_energy = pe_array + ke_array
    
    # Plot energies over time
    axs[1, 2].plot(ts, ke_array, label='KE')
    axs[1, 2].plot(ts, pe_array, label='PE')
    axs[1, 2].plot(ts, total_energy, label='E_total')
    axs[1, 2].set_title('Energy [eV]')
    axs[1, 2].set_xlabel('t [us]')
    axs[1, 2].legend()

    # Hide the empty subplot
    axs[0, 1].axis('off')

    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    main()
