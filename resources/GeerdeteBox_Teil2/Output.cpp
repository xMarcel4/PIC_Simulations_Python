#include "Output.h"

// saves output in VTK format
void Output::fields(World& world, std::vector<Species>& species) {
	std::stringstream name;
	name << "results\\fields_" << std::setfill('0') << std::setw(5) << world.getTs() << ".vti";	// Windows
	//name << "results/fields_" << std::setfill('0') << std::setw(5) << world.getTs() << ".vti";// Linux (others)
	
	// open output file
	std::ofstream out(name.str());
	if (!out.is_open()) { std::cerr << "Could not open " << name.str() << std::endl; return; }

	// ImageData is a VTK format for structured Cartesian meshes
	out << "<VTKFile type=\"ImageData\">\n";
	double3 x0 = world.getX0();
	double3 dh = world.getDh();
	out << "<ImageData Origin=\"" << x0[0] << " " << x0[1] << " " << x0[2] << "\" ";
	out << "Spacing=\"" << dh[0] << " " << dh[1] << " " << dh[2] << "\" ";
	out << "WholeExtent=\"0 " << world.ni - 1 << " 0 " << world.nj - 1 << " 0 " << world.nk - 1 << "\">\n";

	// output data stored on nodes (point data)
	out << "<PointData>\n";

	// node volumes, scalar
	out << "<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out << world.node_vol;	// use the overloaded << operator
	out << "</DataArray>\n";

	// potential, scalar
	out << "<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out << world.phi;
	out << "</DataArray>\n";

	// charge density, scalar
	out << "<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out << world.rho;
	out << "</DataArray>\n";

	// electric field, 3 component vector
	out << "<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out << world.ef;	// uses overloaded << from Field_ and vec3
	out << "</DataArray>\n";

	// species number densities
	for (Species& sp : species) {
		out << "<DataArray Name=\"nd." << sp.name << "\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out << sp.den;
		out << "</DataArray>\n";
	}

	// close the tags
	out << "</PointData>\n";
	out << "</ImageData>\n";
	out << "</VTKFile>\n";
	out.close();
}

// writes information to the screen
void Output::screenOutput(World& world, std::vector<Species>& species) {
	std::cout << "ts: " << world.getTs();
	for (Species& sp : species)
		std::cout << "\t " << sp.name << ": " << sp.getNp();
	std::cout << std::endl;
}

namespace Output { std::ofstream f_diag; }	// file handle
void Output::diagOutput(World& world, std::vector<Species>& species) {
	using namespace Output;	// to get acces to f_diag
	if (!f_diag.is_open()) {			// if file not open
		f_diag.open("runtime_diags.csv");
		f_diag << "ts,time,wall_time";	// write header
		for (Species& sp : species)		// species specific header
			f_diag << ",mp_count." << sp.name << ",real_count." << sp.name << ",px." << sp.name << ",py." << sp.name << ",pz." << sp.name << ",KE." << sp.name;
		f_diag << ",PE,total_E" << std::endl;
	}

	f_diag << world.getTs() << "," << world.getTime();
	f_diag << "," << world.getWallTime();

	// write out species kinetic energy and momentum
	double tot_KE = 0;
	for (Species& sp : species) {
		double KE = sp.getKE();	// species kinetic energy
		tot_KE += KE;			// accumulate total kinetic energy
		double3 mom = sp.getMomentum();	// momentum

		f_diag << "," << sp.getNp() << "," << sp.getRealCount() << "," << mom[0] << "," << mom[1] << "," << mom[2] << "," << KE;
	}

	// write out total potential and kinetic energy
	double PE = world.getPE();
	f_diag << "," << PE << "," << (PE + tot_KE);

	f_diag << "\n";	// use \n to avoid flush to disc
	if (world.getTs() % 10 == 0) f_diag.flush();	// periodically write
}