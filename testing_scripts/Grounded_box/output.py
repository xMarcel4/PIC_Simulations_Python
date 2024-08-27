import os
import numpy as np
from datetime import datetime

class Output:

    @staticmethod
    def create_results_dir(output_dir="results"):
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Directory '{output_dir}' created successfully.")
        except OSError as error:
            print(f"Failed to create directory '{output_dir}': {error}")

    @staticmethod
    def fields(world, species_list, filename_prefix="fields"):
        # The output directory is predefined as 'results'
        output_dir = "results"
        
        # Ensure the output directory exists
        Output.create_results_dir(output_dir)

        # Construct the filename using the timestep
        filename = os.path.join(output_dir, f"{filename_prefix}_{world.getTs():05d}.vti")
        
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

            print(f"Data successfully written to {filename}")

        except IOError as e:
            print(f"Error writing to file {filename}: {e}")

    @staticmethod
    def screen_output(world, species_list):
        print(f"ts: {world.getTs()}")
        for sp in species_list:
            print(f"\t{sp.name}: {sp.getNp()}")
        print("")

    @staticmethod
    def diag_output(world, species_list):
        diag_file = "runtime_diags.csv"
        
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
            line = f"{world.getTs()},{world.getTime()},{world.getWallTime()}"
            
            total_KE = 0
            for sp in species_list:
                KE = sp.getKE()
                total_KE += KE
                mom = sp.getMomentum()
                line += f",{sp.getNp()},{sp.getRealCount()},{mom[0]},{mom[1]},{mom[2]},{KE}"
            
            PE = world.getPE()
            line += f",{PE},{PE + total_KE}\n"
            f_diag.write(line)

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
