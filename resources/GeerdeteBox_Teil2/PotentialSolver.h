#pragma once
#include <iostream>

#include "World.h"
#include "Fields_.h"

class PotentialSolver {
public:
	// constructor, sets members to given inputs
	PotentialSolver(World &world, int max_it, double tol):
		world(world), max_solver_it(max_it), tolerance(tol) { }

	// solves potential using Gauss-Seidel and SOR
	bool solve();

	// compute electric field = -gradient(phi)
	void computeEF();

protected:
	World& world;
	unsigned max_solver_it;	// maximum number of solver iterations
	double tolerance;		// solver tolerance
};