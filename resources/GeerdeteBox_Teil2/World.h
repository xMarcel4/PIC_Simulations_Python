#pragma once
#include <vector>
#include <random>
#include <chrono>

#include "vec3.h"
#include "Fields_.h"
#include "Species.h"

struct Particle;
class Species;
/* to avoid circular dependencies use forward declaration */

namespace Const {
	const double EPS_0 = 8.8541878e-12;	// C/V/m, vac. permittivity
	const double QE = 1.602176565e-19;	// C, electron charge
	const double AMU = 1.660838921e-27;	// kg, atomic mass unit
	const double ME = 9.10938215e-31;	// kg, electron mass
	const double K = 1.380648e-23;		// J/K, Boltzmann constant
	const double PI = 3.141592653;		// pi
	const double EvToK = QE / K;		// 1eV in K ~ 11604
}

class World {
public:
	World(int ni, int nj, int nk);	// constructor

	double3 getX0() const { return double3(x0); }
	double3 getXm() const { return double3(xm); }
	double3 getXc() const { return double3(xc); }
	double3 getDh() const { return double3(dh); }

	double3 XtoL(double3 x) const;
	
	bool inBounds(double3 pos) {
		for (int i = 0; i < 3; i++)
			if (pos[i]<x0[i] || pos[i]>xm[i]) return false;
		return true;
	}
	
	// sets time step and number of time steps
	void setTime(double dt, int num_ts) {this->dt = dt; this->num_ts = num_ts;}

	// functions for accesing time information
	int getTs() const { return ts; }
	double getTime() const { return time; }
	double getDt() const { return dt; }
	bool isLastTimeStep() const { return ts == num_ts - 1; }
	double getWallTime();	// returns elapsed time in seconds

	// advances to the next time step, return true while time remains
	bool advanceTime() { time += dt; ts++; return ts <= num_ts - 1; }

	// sets the mesh span, also recomputes cell spacing
	void setExtents(const double3& _x0, const double3& _xm);

	void computeChargeDensity(std::vector<Species>& species);

	// return potential energy
	double getPE();
	
	const int3 nn;		// number of nodes
	const int ni, nj, nk;	// number of nodes in individual variables

	Field phi;	// potential
	Field rho;	// charge density
	Field3 ef;	// electric field components ??? Wie soll man das initialisieren? Muss eigentlich Field3 sein

	Field node_vol;	// node volumes

protected:
	void computeNodeVolumes();

	double3 x0;	// mesh origin
	double3 dh;	// cell spacing
	double3 xm;	// mesh max bound (opposite corner of x0)
	double3 xc;	// domain centroid

	double dt = 2e-10;	// size time step
	int num_ts = 1;		// number of time steps
	double time = 0;	// current simulation time
	int ts = -1;			// current time step
	std::chrono::time_point<std::chrono::high_resolution_clock>time_start;	// time at simulation start
};


class Rnd {	// object for sampling random numbers
public:
	// constructor: set initial random seed and distribution limits
	Rnd() : mt_gen{ std::random_device() () }, rnd_dist{ 0,1.0 } {}
	double operator() () { return rnd_dist(mt_gen); }

protected:
	std::mt19937 mt_gen;	// random number generator
	std::uniform_real_distribution<double> rnd_dist;	// distribution
};

extern Rnd rnd;	// type Rnd object called rnd defined somewhere