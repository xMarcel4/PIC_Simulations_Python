import os
import numpy as np
from datetime import datetime
from world import *
class Output:

    @staticmethod
    def create_results_dir(output_dir="results"):
        try:
            os.makedirs(output_dir, exist_ok=True)
            # print(f"Directory '{output_dir}' created successfully.")
        except OSError as error:
            return 0
            # print(f"Failed to create directory '{output_dir}': {error}")

    @staticmethod
    def fields(world, species_list, filename_prefix="fields",project_dir=""):
        # The output directory is predefined as 'results'
        output_dir = os.path.join(project_dir, "results")
        
        # Ensure the output directory exists
        Output.create_results_dir(output_dir)
        
        timestamp = datetime.now().strftime("%d_%m_%Y_%H%M%S")

        # Construct the filename using the timestep
        filename = os.path.join(output_dir, f"{filename_prefix}_{timestamp}.vti")
        
        try:
            with open(filename, 'w') as out:
                # Write the VTK header for ImageData
                out.write('<VTKFile type="ImageData">\n')
                
                x0 = world.get_x0()
                dh = world.get_dh()

                # Write the ImageData attributes
                out.write(f'<ImageData Origin="{x0[0]} {x0[1]} {x0[2]}" ')
                out.write(f'Spacing="{dh[0]} {dh[1]} {dh[2]}" ')
                out.write(f'WholeExtent="0 {world.ni - 1} 0 {world.nj - 1} 0 {world.nk - 1}">\n')

                # Start PointData section
                out.write('<PointData>\n')

                # Write node_vol data
                out.write('<DataArray Name="NodeVol" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.node_vol))
                out.write('</DataArray>\n')

                # Write phi data
                out.write('<DataArray Name="phi" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.phi))
                out.write('</DataArray>\n')

                # Write rho data
                out.write('<DataArray Name="rho" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.rho))
                out.write('</DataArray>\n')

                # Write ef data (vector field)
                out.write('<DataArray Name="ef" NumberOfComponents="3" format="ascii" type="Float64">\n')
                out.write(Output.vector_field_data_to_string(world.ef))
                out.write('</DataArray>\n')

                # Write species number densities
                for sp in species_list:
                    out.write(f'<DataArray Name="nd.{sp.name}" NumberOfComponents="1" format="ascii" type="Float64">\n')
                    out.write(Output.field_data_to_string(sp.den))
                    out.write('</DataArray>\n')

                # Close tags
                out.write('</PointData>\n')
                out.write('</ImageData>\n')
                out.write('</VTKFile>\n')

            #print(f"Data successfully written to {filename}")

        except IOError as e:
            print(f"Error writing to file {filename}: {e}")
            
    def fields_without_species(world, filename_prefix, project_dir=""):
    # The output directory is predefined as 'results'
        output_dir = os.path.join(project_dir, "results")
        
        # Ensure the output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
        # Create a timestamp for the filename
        timestamp = datetime.now().strftime("%d_%m_%Y_%H%M%S")
    
        # Construct the filename using the prefix and timestamp
        filename = os.path.join(output_dir, f"{filename_prefix}_{timestamp}.vti")
    
        try:
            with open(filename, 'w') as out:
                # Write the VTK header for ImageData
                out.write('<VTKFile type="ImageData">\n')
                
                x0 = world.get_x0()
                dh = world.get_dh()
    
                # Write the ImageData attributes
                out.write(f'<ImageData Origin="{x0[0]} {x0[1]} {x0[2]}" ')
                out.write(f'Spacing="{dh[0]} {dh[1]} {dh[2]}" ')
                out.write(f'WholeExtent="0 {world.ni - 1} 0 {world.nj - 1} 0 {world.nk - 1}">\n')
    
                # Start PointData section
                out.write('<PointData>\n')
    
                # Write node_vol data
                out.write('<DataArray Name="NodeVol" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.node_vol))
                out.write('</DataArray>\n')
    
                # Write phi data
                out.write('<DataArray Name="phi" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.phi))
                out.write('</DataArray>\n')
    
                # Write rho data
                out.write('<DataArray Name="rho" NumberOfComponents="1" format="ascii" type="Float64">\n')
                out.write(Output.field_data_to_string(world.rho))
                out.write('</DataArray>\n')
    
                # Write ef data (vector field)
                out.write('<DataArray Name="ef" NumberOfComponents="3" format="ascii" type="Float64">\n')
                out.write(Output.vector_field_data_to_string(world.ef))
                out.write('</DataArray>\n')
    
                # Close tags
                out.write('</PointData>\n')
                out.write('</ImageData>\n')
                out.write('</VTKFile>\n')
    
            print(f"Data successfully written to {filename}")
    
        except Exception as e:
            print(f"Could not write data to {filename}: {e}")
    
        
    @staticmethod
    def screen_output(world, species_list):
        output = f"ts: {world.get_ts()}"
        for sp in species_list:
            real_count = int(round(sp.get_real_count()))  # Ensure it's an integer
            output += f"  {sp.name}: {real_count}"
        print(output)

    @staticmethod
    def diag_output(world, species_list, project_dir=""):
        # print(f"Iteration: {world.get_ts()}")
        diag_file = os.path.join(project_dir, "runtime_diags.csv")
        
        # Check if the file exists to write the header
        write_header = not os.path.exists(diag_file)
        
        with open(diag_file, 'a') as f_diag:
            if write_header:
                # Write header
                header = "ts,time,wall_time"
                for sp in species_list:
                    header += f",mp_count.{sp.name},real_count.{sp.name},px.{sp.name},py.{sp.name},pz.{sp.name},KE.{sp.name}"
                header += ",PE,total_E\n"
                f_diag.write(header)
            
            # Write data
            line = f"{world.get_ts()},{world.get_time()},{world.get_wall_time()}"
            
            total_KE = 0
            for sp in species_list:
                KE = sp.get_ke()
                total_KE += KE
                mom = sp.get_momentum()
                line += f",{sp.get_np()},{sp.get_real_count()},{mom[0]},{mom[1]},{mom[2]},{KE}"
            
            PE = world.get_pe()
            total_E = PE + total_KE
            line += f",{PE},{total_E}\n"
            f_diag.write(line)
        
        # # Print the energy values to the console
        # print(f"Total KE: {total_KE}")
        # print(f"PE: {PE}")
        # print(f"Total Energy: {total_E}")

       
    # @staticmethod
    # def diag_output(world, species_list, project_dir=""):
    #     diag_file = os.path.join(project_dir, "runtime_diags.csv")
        
    #     # Check if the file exists to write the header
    #     write_header = not os.path.exists(diag_file)
        
    #     with open(diag_file, 'a') as f_diag:
    #         if write_header:
    #             # Write header
    #             header = "ts,time,wall_time"
    #             for sp in species_list:
    #                 header += f",mp_count.{sp.name},real_count.{sp.name},px.{sp.name},py.{sp.name},pz.{sp.name},KE.{sp.name}"
    #             header += ",PE,total_E\n"
    #             f_diag.write(header)
            
    #         # Write data
    #         line = f"{world.get_ts()},{world.get_time()},{world.get_wall_time()}"
            
    #         total_KE = 0
    #         for sp in species_list:
    #             KE = sp.get_ke()
    #             total_KE += KE
    #             mom = sp.get_momentum()
    #             line += f",{sp.get_np()},{sp.get_real_count()},{mom[0]},{mom[1]},{mom[2]},{KE}"
            
    #         PE = world.get_pe()
    #         line += f",{PE},{PE + total_KE}\n"
    #         f_diag.write(line)

    @staticmethod
    def field_data_to_string(field):
        """Converts scalar field data to a string in VTK format."""
        data_str = ""
        for k in range(field.nk):
            for j in range(field.nj):
                for i in range(field.ni):
                    data_str += f"{field(i, j, k)} "
        return data_str + "\n"

    @staticmethod
    def vector_field_data_to_string(field):
        """Converts vector field data to a string in VTK format."""
        data_str = ""
        for k in range(field.nk):
            for j in range(field.nj):
                for i in range(field.ni):
                    data_str += f"{field(i, j, k, 0)} {field(i, j, k, 1)} {field(i, j, k, 2)} "
        return data_str + "\n"
