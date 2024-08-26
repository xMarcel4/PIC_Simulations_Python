#include <iostream>

#include "vec3.h"
#include "Fields_.h"
#include "World.h"
#include "Species.h"
#include "PotentialSolver.h"
#include "Output.h"

int main() {
	// initialize domain
	World world(21, 21, 21);	// mesh size
	world.setExtents({ -0.1,-0.1,0.0 }, { 0.1,0.1,0.2 });	// start and end of world
	world.setTime(2e-9, 1000);	// time step size and number of time steps

	// set up particle species
	std::vector<Species> species;
	species.reserve(2);	// pre-allocate space for two species
	species.push_back(Species("O+", 16 * Const::AMU, Const::QE, world));
	species.push_back(Species("e-", Const::ME, -1 * Const::QE, world));

	// initialize potential solver and solve initial potential
	PotentialSolver solver(world, 10000, 1e-4);
	solver.solve();

	// obtain initial electric field
	solver.computeEF();

	// create particles
	//// normal start
	//int np_ions = 80000;	// number of simulation ions
	//int np_eles = 10000;	// number of simulation electrons
	//species[0].loadParticlesBox(world.getX0(), world.getXm(), 1e11, np_ions);	// ions
	//species[1].loadParticlesBox(world.getX0(), world.getXc(), 1e11, np_eles);	// electrons
	// quiet start
	int3 np_ions_grid = { 41,41,41 };
	int3 np_eles_grid = { 21,21,21 };
	species[0].loadParticlesBoxQS(world.getX0(), world.getXm(), 1e11, np_ions_grid);	// ions
	species[1].loadParticlesBoxQS(world.getX0(), world.getXc(), 1e11, np_eles_grid);	// electrons

	// main loop
	while (world.advanceTime()) {
		// move particles
		for (Species& sp : species) {
			sp.advance();
			sp.computeNumberDensity();
		}

		// compute charge densitiy
		world.computeChargeDensity(species);

		// update potential
		solver.solve();

		// obtain electric field
		solver.computeEF();

		// screen and file output
		if (world.getTs() % 10 == 0 || world.isLastTimeStep()) Output::screenOutput(world, species);
		Output::diagOutput(world, species);

		// periodically write out results
		if (world.getTs() % 10 == 0 || world.isLastTimeStep())
			Output::fields(world, species);
	}

	// output run time
	std::cout << "Simulation took " << world.getWallTime() << " seconds" << std::endl;

	return 0;	// indicate normal exit
}