#include "PotentialSolver.h"

bool PotentialSolver::solve() {
	Field& phi = world.phi;	// refernces to avoid writing world.phi
	Field& rho = world.rho;

	// precompute 1/(dx^2)
	double3 dh = world.getDh();
	double idx2 = 1.0 / (dh[0] * dh[0]);	// 1/dx^2
	double idy2 = 1.0 / (dh[1] * dh[1]);	// 1/dy^2
	double idz2 = 1.0 / (dh[2] * dh[2]);	// 1/dz^2

	double L2 = 0;	// norm
	bool converged = false;

	// solve potential
	for (unsigned it = 0; it < max_solver_it; it++) {
		for (int i = 1; i < world.ni - 1; i++)
			for (int j = 1; j < world.nj - 1; j++)
				for (int k = 1; k < world.nk - 1; k++) {
					// standard internal open node
					double phi_new = (rho(i,j,k) / Const::EPS_0 +
						idx2 * (phi(i - 1,j,k) + phi(i + 1,j,k)) +
						idy2 * (phi(i,j - 1,k) + phi(i,j + 1,k)) +
						idz2 * (phi(i,j,k - 1) + phi(i,j,k + 1))) /
						(2 * idx2 + 2 * idy2 + 2 * idz2);

					// SOR
					phi.w_at(i, j, k) = phi(i, j, k) + 1.4 * (phi_new - phi(i, j, k));
				}

		// check for convergence
		if (it % 25 == 0) {
			double sum = 0;
			for (int i = 1; i < world.ni - 1; i++)
				for (int j = 1; j < world.nj - 1; j++)
					for (int k = 1; k < world.nk - 1; k++) {
						double R = -phi(i,j,k) * (2 * idx2 + 2 * idy2 + 2 * idz2) +
							rho(i,j,k) / Const::EPS_0 +
							idx2 * (phi(i - 1,j,k) + phi(i + 1,j,k)) +
							idy2 * (phi(i,j - 1,k) + phi(i,j + 1,k)) +
							idz2 * (phi(i,j,k - 1) + phi(i,j,k + 1));

						sum += R * R;
					}
			L2 = sqrt(sum / (static_cast<double>(world.ni) * static_cast<double>(world.nj) * static_cast<double>(world.nk)));
			if (L2 < tolerance) { converged = true; break; }
		}
	}

	if (!converged) std::cerr << "GS failed to converge, L2=" << L2 << std::endl;
	return converged;
}

// computes electric field = -gradient(phi) using 2nd order differencing
void PotentialSolver::computeEF() {
	// reference to phi to avoid writing world.phi
	Field& phi = world.phi;

	double3 dh = world.getDh();		// get cell spacing
	double dx = dh[0];
	double dy = dh[1];
	double dz = dh[2];

	for (int i = 0; i < world.ni; i++)				// loop over nodes
		for (int j = 0; j < world.nj; j++)
			for (int k = 0; k < world.nk; k++) {
				double3& ef = world.ef[i + j * world.ni + k * world.ni * world.nj];	// ref to (i,j,k) ef vec3

				// x component, efx
				if (i == 0)					// forward difference
					ef[0] = (3 * phi(i, j, k) - 4 * phi(i + 1, j, k) + phi(i + 2, j, k)) / (2 * dx);
				else if (i == world.ni - 1)	//backward difference
					ef[0] = (4 * phi(i - 1, j, k) - phi(i - 2, j, k) - 3 * phi(i, j, k)) / (2 * dx);
				else						// central difference
					ef[0] = (phi(i - 1, j, k) - phi(i + 1, j, k)) / (2 * dx);

				// y component, efy
				if (j == 0)					// forward difference
					ef[1] = (3 * phi(i, j, k) - 4 * phi(i, j + 1, k) + phi(i, j + 2, k)) / (2 * dy);
				else if (j == world.nj - 1)	//backward difference
					ef[1] = (4 * phi(i, j - 1, k) - phi(i, j - 2, k) - 3 * phi(i, j, k)) / (2 * dy);
				else						// central difference
					ef[1] = (phi(i, j - 1, k) - phi(i, j + 1, k)) / (2 * dy);

				// z component, efz
				if (k == 0)					// forward difference
					ef[2] = (3 * phi(i, j, k) - 4 * phi(i, j, k + 1) + phi(i, j, k + 2)) / (2 * dz);
				else if (k == world.nk - 1)	//backward difference
					ef[2] = (4 * phi(i, j, k - 1) - phi(i, j, k - 2) - 3 * phi(i, j, k)) / (2 * dz);
				else						// central difference
					ef[2] = (phi(i, j, k - 1) - phi(i, j, k + 1)) / (2 * dz);
			}
}