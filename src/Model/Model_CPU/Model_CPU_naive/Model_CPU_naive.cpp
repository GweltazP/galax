#include <cmath>

#include "Model_CPU_naive.hpp"
#include "Octotree.hpp"

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

	Boundary b = new Boundary({0.0f, 0.0f, 0.0f, max(particles.x, particles.y, particles.z)}); // REGARDER COMMENT CALCULER LA TAILLE DE LA BOITE DE BARNES-HUT
	//TO FIX

	Octotree Octree = new Octotree(b); 

	for (int i = 0; i < n_particles; i++)
	{
		Point p = new Point({particles.x[i], particles.y[i], particles.z[i], initstate.masses[i]});
		Octree.include(p);
	}

	for (int i = 0; i < n_particles; i++)
	{
		for (int j; j < n_particles; j++)
		{
			/*PSEUDO CODE
			Si distance (particule[i] -> CoM du noeud) < theta
				On fait le calcule avec le centre de masse
			Sinon
				On descende plus bas dans l'arbre pour faire le calcule.
			*/
		}
	}
	//POSSIBILITER DE PARALLÉLLISME SUR LES BOUCLE FOR, 8 THREADS POUR LE PARCOUR DE L'ARBRE/SOUS-ARBRE COMME ACCÈDE PARTIE DE MÉMOIRE NON CONCURENTE

	Octree.~Octotree();

}