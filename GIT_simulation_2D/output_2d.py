import os
import numpy as np
from datetime import datetime
from world_2d import *

class Output2D:
    
    @staticmethod
    def fields_without_species_2d(world, filename_prefix, project_dir=""):
        output_dir = os.path.join(project_dir, "results")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
        timestamp = datetime.now().strftime("%d_%m_%Y_%H%M%S")
        filename = os.path.join(output_dir, f"{filename_prefix}_{timestamp}.vti")
    
        try:
            with open(filename, 'w') as out:
                # Ensure valid VTK header
                out.write('<?xml version="1.0"?>\n')
                out.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
                x0 = world.get_x0()
                dh = world.get_dh()
                ni, nj = world.ni, world.nj
                
                # Write ImageData attributes
                out.write(f'<ImageData Origin="{x0[0]} {x0[1]}" ')
                out.write(f'Spacing="{dh[0]} {dh[1]}" ')
                # Correctly define the WholeExtent for 2D (third dimension is 0)
                out.write(f'WholeExtent="0 {ni - 1} 0 {nj - 1} 0 0">\n')
                out.write('<PointData>\n')
    
                # Debugging lengths
                print(f"object_id length: {len(world.objectid.data.flatten())}")
                print(f"phi length: {len(world.phi.data.flatten())}")
                print(f"rho length: {len(world.rho.data.flatten())}")
                print(f"ef length: {len(world.ef.data.flatten())}")
    
                world.objectid.data = world.objectid.data.astype(int)
    
                # Write the object_id data
                out.write('<DataArray Name="object_id" NumberOfComponents="1" format="ascii" type="Int32">\n')
                out.write(Output2D.field_data_to_string(world.objectid))
                out.write('</DataArray>\n')
    
                out.write('<DataArray Name="NodeVol" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.node_vol))
                out.write('</DataArray>\n')
    
                out.write('<DataArray Name="phi" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.phi))
                out.write('</DataArray>\n')
    
                out.write('<DataArray Name="rho" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.rho))
                out.write('</DataArray>\n')
    
                # Write the ef data with 2 components for 2D fields
                out.write('<DataArray Name="ef" NumberOfComponents="2" format="ascii" type="Float64">\n')
                out.write(Output2D.vector_field_data_to_string(world.ef))
                out.write('</DataArray>\n')
    
                out.write('</PointData>\n')
                out.write('</ImageData>\n')
                out.write('</VTKFile>\n')  # Ensure the closing tag for VTKFile
    
            print(f"Data successfully written to {filename}")
    
        except Exception as e:
            print(f"Could not write data to {filename}: {e}")

    @staticmethod
    def field_data_to_string(field):
        """Converts scalar field data to a string in VTK format."""
        data_str = ""
        for j in range(field.nj):
            for i in range(field.ni):
                data_str += f"{field(i, j)} "
        return data_str + "\n"
    
    @staticmethod
    def vector_field_data_to_string(field):
        """Converts vector field data to a string in VTK format."""
        data_str = ""
        for j in range(field.nj):
            for i in range(field.ni):
                data_str += f"{field(i, j, 0)} {field(i, j, 1)} "
        return data_str + "\n"

    @staticmethod
    def create_results_dir(output_dir="results"):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as error:
            return 0
    
    @staticmethod
    def screen_output(world, species_list):
        print(f"ts: {world.get_ts()}", end="")
        
        for sp in species_list:
            # Display the number of particles for each species
            print(f"\t {sp.name}: {sp.get_np():.3f}", end="")
    
            # Print the number of particles out of bounds and in the plane
            print(f"\t Particles out of bounds: {sp.world.particles_out_of_bounds}", end="")
            print(f"\t Particles in plane: {sp.world.particles_in_plane}", end="")
        
        print("")  # Print a newline at the end


    @staticmethod
    def diag_output(world, species_list, project_dir=""):
        diag_file = os.path.join(project_dir, "runtime_diags.csv")
        write_header = not os.path.exists(diag_file)
        
        with open(diag_file, 'a') as f_diag:
            if write_header:
                header = "ts,time,wall_time"
                for sp in species_list:
                    header += f",mp_count.{sp.name},real_count.{sp.name},px.{sp.name},py.{sp.name},KE.{sp.name}"
                header += ",PE,total_E\n"
                f_diag.write(header)
            
            line = f"{world.get_ts()},{world.get_time()},{world.get_wall_time()}"
            total_KE = 0
            for sp in species_list:
                KE = sp.get_ke()
                total_KE += KE
                mom = sp.get_momentum()
                line += f",{sp.get_np()},{sp.get_real_count()},{mom[0]},{mom[1]},{KE}"
            
            PE = world.get_pe()
            total_E = PE + total_KE
            line += f",{PE},{total_E}\n"
            f_diag.write(line)

        
    @staticmethod
    def fields(world, species_list, filename_prefix="fields", project_dir=""):
        # Create the results directory
        output_dir = os.path.join(project_dir, "results")
        Output2D.create_results_dir(output_dir)
        
        # Create a timestamp for the filename
        timestamp = datetime.now().strftime("%d_%m_%Y_%H%M%S")
        filename = os.path.join(output_dir, f"{filename_prefix}_{timestamp}.vti")
        
        try:
            # Open the file for writing
            with open(filename, 'w') as out:
                # Write VTK header
                out.write('<?xml version="1.0"?>\n')
                out.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
                
                # Get world grid spacing and origin
                x0 = world.get_x0()
                dh = world.get_dh()
                ni, nj = world.ni, world.nj
                
                # Define the ImageData for 2D
                out.write(f'<ImageData Origin="{x0[0]} {x0[1]}" ')
                out.write(f'Spacing="{dh[0]} {dh[1]}" ')
                out.write(f'WholeExtent="0 {ni - 1} 0 {nj - 1} 0 0">\n')
                out.write('<PointData>\n')
    
                world.objectid.data = world.objectid.data.astype(int)
    
                # Write the object_id data
                out.write('<DataArray Name="object_id" NumberOfComponents="1" format="ascii" type="Int32">\n')
                out.write(Output2D.field_data_to_string(world.objectid))
                out.write('</DataArray>\n')
    
                # Write the node volume data (2D)
                out.write('<DataArray Name="NodeVol" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.node_vol))
                out.write('</DataArray>\n')
    
                # Write the electric potential (phi) data
                out.write('<DataArray Name="phi" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.phi))
                out.write('</DataArray>\n')
    
                # Write the charge density (rho) data
                out.write('<DataArray Name="rho" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output2D.field_data_to_string(world.rho))
                out.write('</DataArray>\n')
    
                # Write the electric field data (2D vector field)
                out.write('<DataArray Name="ef" NumberOfComponents="2" format="ascii" type="Float64">\n')
                out.write(Output2D.vector_field_data_to_string(world.ef))
                out.write('</DataArray>\n')
    
                # Loop through each species and write the number density for each species
                for sp in species_list:
                    out.write(f'<DataArray Name="nd.{sp.name}" NumberOfComponents="1" format="ascii" type="Float64">\n')
                    out.write(Output2D.field_data_to_string(sp.den))
                    out.write('</DataArray>\n')
    
                # Close PointData and ImageData sections
                out.write('</PointData>\n')
                out.write('</ImageData>\n')
                out.write('</VTKFile>\n')  # Close the VTKFile tag
    
            # Print success message
            print(f"Data successfully written to {filename}")
    
        except IOError as e:
            print(f"Error writing to file {filename}: {e}")
