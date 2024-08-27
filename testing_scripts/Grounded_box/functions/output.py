import os
import numpy as np

class Output:

    @staticmethod
    def create_results_dir(output_dir="results"):
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Directory '{output_dir}' created successfully.")
        except OSError as error:
            print(f"Failed to create directory '{output_dir}': {error}")

    @staticmethod
    def fields(world, output_dir="results"):
        # Ensure the output directory exists
        Output.create_results_dir(output_dir)

        # Construct the full path for the output file
        filename = os.path.join(output_dir, "fields.vti")

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

        except IOError as e:
            print(f"Error writing to file {filename}: {e}")

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
