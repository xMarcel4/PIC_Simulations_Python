# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:07:28 2024

@author: marce
"""

import os
import shutil

# Define the source and destination folders
source_folder = r"C:\Users\Marcel\OneDrive\Dokumente\GitHub_repositories\GIT_particle_trajectory_simulation\testing_scripts\GIT_simulation_2D\Test_loop__xenon_v5000_scale_1000_den13\\\\results"

# Construct the full path to the source and destination files
destination_folder = os.path.join(source_folder, "renamed_files")

# Create the destination folder if it doesn't exist
os.makedirs(destination_folder, exist_ok=True)

# Iterate through all files in the source folder
for filename in os.listdir(source_folder):
    if filename.endswith(".vti"):
        # Extract the ts number from the filename
        ts_index = filename.find("ts")
        if ts_index != -1:
            ts_number = filename[ts_index + 2:ts_index + 7]  # Extract the ts number

            # Construct the new filename in the format "fields_00000.vti"
            new_filename = f"fields_{ts_number}.vti"

            # Construct the full path to the source and destination files
            src_file = os.path.join(source_folder, filename)
            dst_file = os.path.join(destination_folder, new_filename)

            # Copy the file to the new destination with the new name
            shutil.copy(src_file, dst_file)

            print(f"Renamed and copied {filename} to {new_filename}")

print("All files have been renamed and copied.")
