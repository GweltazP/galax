#include <cmath>

#include "Model_CPU_naive.hpp"

Model_CPU_naive
::Model_CPU_naive(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
}

void Model_CPU_naive
::step()
{
	std::fill(accelerationsx.begin(), accelerationsx.end(), 0);
	std::fill(accelerationsy.begin(), accelerationsy.end(), 0);
	std::fill(accelerationsz.begin(), accelerationsz.end(), 0);

	Octotree Octree = new Octotree({0, 0, 0, max(particles.x, particles.y, particles.z)}); // REGARDER COMMENT CALCULER LA TAILLE DE LA BOITE DE BARNES-HUT

	for (int i = 0; i < n_particles; i++)
	{
		Point p = new Point({particles.x[i], particles.y[i], particles.z[i], initstate.masses[i]});
		Octree.include(p);
	}
	//POSSIBILITER DE PARALLÉLLISME SUR LES BOUCLE FOR, 8 THREADS POUR LE PARCOUR DE L'ARBRE/SOUS-ARBRE COMME ACCÈDE PARTIE DE MÉMOIRE NON CONCURENTE

}