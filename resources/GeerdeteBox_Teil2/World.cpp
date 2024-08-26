#include "World.h"

// constructor
World::World(int ni, int nj, int nk) :
	ni{ ni }, nj{ nj }, nk{ nk }, nn{ ni,nj,nk },
	phi(ni, nj, nk), rho(ni, nj, nk), node_vol(ni, nj, nk), ef(ni, nj, nk) {
		time_start = std::chrono::high_resolution_clock::now();	// save starting time point
}

Rnd rnd;	// create an instance of a Rnd object

double3 World::XtoL(double3 x) const {
	double3 lc;
	lc[0] = (x[0] - x0(0)) / dh(0);
	lc[1] = (x[1] - x0(1)) / dh(1);
	lc[2] = (x[2] - x0(2)) / dh(2);
	return lc;
}

double World::getWallTime() {
	auto time_now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_delta = time_now - time_start;
	return time_delta.count();
}

void World::computeNodeVolumes() {
	for (int i = 0; i < ni; i++)	// loop over nodes
		for (int j = 0; j < nj; j++)
			for (int k = 0; k < nk; k++) {
				double V = dh[0] * dh[1] * dh[2];		// standard volume
				if (i == 0 || i == ni - 1) V *= 0.5;	// adjust on boundaries
				if (j == 0 || j == nj - 1) V *= 0.5;
				if (k == 0 || k == nk - 1) V *= 0.5;
				node_vol[i+j*ni+k*ni*nj] = V;
			}
}

void World::setExtents(const double3 &_x0, const double3 &_xm) {
	x0 = _x0;	// set our copy of the origin
	xm = _xm;	// do the same for xmax

	for (int i = 0; i < 3; i++)
		dh[i] = (xm[i] - x0[i])/(static_cast<double>(nn(i)) - 1.0);

	// compute centroid
	xc = 0.5 * (x0 + xm);

	computeNodeVolumes();
}

void World::computeChargeDensity(std::vector<Species>& species) {
	rho = 0;
	for (Species& sp : species) {		// loop over species
		if (sp.charge == 0) continue;	// don't bother with neutrals
		rho += sp.charge * sp.den;		// accumulate density scaled by charge
	}
}

double World::getPE() {
	double pe = 0;
	for (int index = 0; index < ni * nj * nk; index++) {
		double3 efn = ef(index);	// ef at this node
		double ef2 = efn[0] * efn[0] + efn[1] * efn[1] + efn[2] * efn[2];
		pe += ef2 * node_vol(index);
	}
	return 0.5 * Const::EPS_0 * pe;
}
