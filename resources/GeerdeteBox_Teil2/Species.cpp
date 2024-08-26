#include "Species.h"

Species::Species(std::string name, double mass, double charge, World& world):
		name{ name }, mass{ mass }, charge{ charge }, den(world.ni, world.nj, world.nk), world{ world } { }

double Species::getRealCount() {
	double mpw_sum = 0;
	for (Particle& part : particles) mpw_sum += part.mpw;
	return mpw_sum;
}

double3 Species::getMomentum() {
	double3 mom;
	for (Particle& part : particles) mom += part.mpw * part.vel;
	return mass*mom;
}

double Species::getKE() {
	double ke = 0;
	for (Particle& part : particles) {
		double v2 = part.vel[0] * part.vel[0] + part.vel[1] * part.vel[1] + part.vel[2] * part.vel[2];
		ke += part.mpw * v2;
	}
	return 0.5*mass*ke;
}

void Species::advance() {
	// get the time step
	double dt = world.getDt();

	// get mesh bounds
	double3 x0 = world.getX0();
	double3 xm = world.getXm();

	// loop over all particles
	for (Particle& part : particles) {
		// get logical coordinate of particle's postition
		double3 lc = world.XtoL(part.pos);

		// electric field at particle position
		double3 ef_part = world.ef.gather(lc);

		// update velocity from F = qE
		part.vel += ef_part * (dt * charge / mass);

		// update position from v = dx/dt
		part.pos += part.vel * dt;

		// reflect particles leaving the domain
		for (int i = 0; i < 3; i++) {
			if (part.pos[i] < x0[i]) { part.pos[i] = 2 * x0[i] - part.pos[i]; part.vel[i] *= -1.0; }
			else if (part.pos[i] >= xm[i]) { part.pos[i] = 2 * xm[i] - part.pos[i]; part.vel[i] *= -1.0; }
		}
	}
}

void Species::computeNumberDensity() {
	den = 0;	// set all values to zero
	for (Particle &part:particles) {		// loop over particles
		double3 lc = world.XtoL(part.pos);	// get logical coordinates
		den.scatter(lc, part.mpw);			// deposit weight
	}
	den /= world.node_vol;					// divide by node volume
}

void Species::addParticle(double3 pos, double3 vel, double mpw) {
	// don't do anything if pos outside domain bounds [x0,xd)
	if (!world.inBounds(pos)) return;

	// get particle logical coordinate
	double3 lc = world.XtoL(pos);

	// evaluate electric field at particle position
	double3 ef_part = world.ef.gather(lc);

	// rewind velocity by 0.5*dt*ef
	vel -= charge / mass * ef_part * (0.5 * world.getDt());

	// add to list
	particles.emplace_back(pos, vel, mpw);
}

void Species::loadParticlesBox(double3 x1, double3 x2, double num_den, int num_sim) {
	double box_vol = (x2[0] - x1[0]) * (x2[1] - x1[1]) * (x2[2] - x1[2]);	// box volume
	double num_real = num_den * box_vol;	// number of real particles
	double mpw = num_real / num_sim;	// macroparticle weight

	// preallocate memory for particles
	particles.reserve(num_sim);

	// load particles on an equally spaced grid
	for (int p = 0; p < num_sim; p++) {
		// sample random position
		double3 pos;
		pos[0] = x1[0] + rnd() * (x2[0] - x1[0]);
		pos[1] = x1[1] + rnd() * (x2[1] - x1[1]);
		pos[2] = x1[2] + rnd() * (x2[2] - x1[2]);

		// set initial velocity
		double3 vel{ 0,0,0 };	// stationary particle

		addParticle(pos, vel, mpw);	// add a new particle to the array
	}
}

void Species::loadParticlesBoxQS(double3 x1, double3 x2, double num_den, int3 num_sim) {
	double box_vol = (x2[0] - x1[0]) * (x2[1] - x1[1]) * (x2[2] - x1[2]);	// box volume
	int num_sim_tot = (num_sim[0] - 1) * (num_sim[1] - 1) * (num_sim[2] - 1);	// total number of simulation particles
	double num_real = num_den * box_vol;	// number of real particles
	double mpw = num_real / num_sim_tot;	// macroparticle weight

	// compute particle grid spacing
	double di = (x2[0] - x1[0]) / (static_cast<double>(num_sim[0]) - 1);
	double dj = (x2[1] - x1[1]) / (static_cast<double>(num_sim[1]) - 1);
	double dk = (x2[2] - x1[2]) / (static_cast<double>(num_sim[2]) - 1);
	
	// preallocate memory for particles
	particles.reserve(num_sim_tot);

	// load particles on an equally spaced grid
	for (int i = 0; i < num_sim[0]; i++)
		for (int j = 0; j < num_sim[1]; j++)
			for (int k = 0; k < num_sim[2]; k++) {
				double3 pos;
				pos[0] = x1[0] + i * di;
				pos[1] = x1[1] + j * dj;
				pos[2] = x1[2] + k * dk;

				// shift particles on max faces back to the domain
				if (pos[0] == x2[0]) pos[0] -= 1e-4 * di;
				if (pos[1] == x2[1]) pos[1] -= 1e-4 * dj;
				if (pos[2] == x2[2]) pos[2] -= 1e-4 * dk;

				double w = 1;	// relative weight
				if (i == 0 || i == num_sim[0] - 1) w *= 0.5;
				if (j == 0 || j == num_sim[1] - 1) w *= 0.5;
				if (k == 0 || k == num_sim[2] - 1) w *= 0.5;

				// add rewind
				double vel[3] = { 0,0,0 };	// particle is stationary

				addParticle(pos, vel, mpw * w);	// add to array
			}
}