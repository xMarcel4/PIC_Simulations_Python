#pragma once
#include <vector>

#include "World.h"

class World;	// forward declaration
class Rnd;

struct Particle {
	double3 pos;	// position
	double3 vel;	// velocity
	double mpw;		// macroparticle weight
	
	Particle(double3 x, double3 v, double mpw) :	// constructor
		pos{ x[0], x[1], x[2] }, vel{ v[0], v[1], v[2] }, mpw{ mpw } { }
};


class Species {
public:
	// constructor
	Species(std::string name, double mass, double charge, World& world);
	/* had to be moved out of class definition (avoidance of circular dependencies) */

	// returns the number of particles
	size_t getNp() const { return particles.size(); }	// simulation particles
	double getRealCount();	// real particles

	// returns momentum
	double3 getMomentum();

	// returns kinetic energy
	double getKE();

	// moves all particles using electric field ef[]
	void advance();

	// compute number density
	void computeNumberDensity();

	// adds a new particle
	void addParticle(double3 pos, double3 vel, double mpw);

	// loads num_sim simulated particles in a x1-x2 box with num_den number density
	void loadParticlesBox(double3 x1, double3 x2, double num_den, int num_sim);	// random distribution
	void loadParticlesBoxQS(double3 x1, double3 x2, double num_den, int3 num_sim);	// quiet start

	const std:: string name;	// species name
	const double mass;			// particle mass in kg
	const double charge;		// particle charge in C

	Field den;							// number density
	std::vector<Particle> particles;	// array for storing particles

protected:
	World& world;	// reference to the World object
};