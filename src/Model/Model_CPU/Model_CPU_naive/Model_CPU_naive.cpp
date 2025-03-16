#include <cmath>

#include "Model_CPU_naive.hpp"
#include <iostream>
#include <xsimd/xsimd.hpp>
#include <omp.h>

namespace xs = xsimd;
using b_type = xs::batch<float, xs::avx2>;

void Octree::delete_tree(OctreeNode* node) {
    if (!node) return;
    for (int i = 0; i < 8; i++) {
        delete_tree(node->children[i]);
    }
    delete node;
}

void Octree::insert_particle(int index, const Particles& particles, const Initstate& initstate) {
    insert_recursive(root, index, particles, initstate);
}

void Octree::insert_recursive(OctreeNode* node, int index, const Particles& particles, const Initstate& initstate) {
    float px = particles.x[index];
    float py = particles.y[index];
    float pz = particles.z[index];
    if (node->particle_index == -1 && node->children[0] == nullptr) {
        node->particle_index = index;
        node->x = px;
        node->y = py;
        node->z = pz;
        node->mass = initstate.masses[index];
    } else {
        if (node->children[0] == nullptr) {
            float halfSize = node->size / 2.0f;
            node->subdivide();
            int oldIndex = node->particle_index;
            node->particle_index = -1;
            // insert_recursive(root, oldIndex, particles, initstate.masses[index]);
        }
        int childIndex = (px > node->x) | ((py > node->y) << 1) | ((pz > node->z) << 2);
        insert_recursive(node->children[childIndex], index, particles, initstate);
    
        node->x = (node->x * node->mass + px * initstate.masses[index]) / (node->mass + initstate.masses[index]);
        node->y = (node->y * node->mass + py * initstate.masses[index]) / (node->mass + initstate.masses[index]);
        node->z = (node->z * node->mass + pz * initstate.masses[index]) / (node->mass + initstate.masses[index]);
        node->mass += initstate.masses[index];
    }
}

void Octree::compute_force(int index, std::vector<float>& ax, std::vector<float>& ay, std::vector<float>& az, const Particles& particles, const Initstate& initstate) {
    compute_force_recursive(root, index, ax, ay, az, particles, initstate);
}

void Octree::compute_force_recursive(OctreeNode* node, int index, std::vector<float>& ax, std::vector<float>& ay, std::vector<float>& az, const Particles& particles, const Initstate& initstate) {
    // std::cout << "Calcul récursif sur la particule " << index << std::endl;
	// std::cout << "Index du nœud courrant " << (node->particle_index) << std::endl;
	// std::cout << "Acceleration x/y/z " << ax[index] << "/" << ay[index] << "/" << az[index] << std::endl;
	if (!node || (node->particle_index == index) || (node->mass == 0)) {
        // std::cout << "La particule n'applique pas de force sur elle même " << std::endl;
		// std::cout << std::endl;
		return;
    }

	// std::cout << "x " << particles.x[index] << std::endl;
    // std::cout << "y " << particles.y[index] << std::endl;
    // std::cout << "z " << particles.z[index] << std::endl;
	// std::cout << "Centre du noeud courrant :" << (node->center_x) << "/" << (node->center_y) << "/" << (node->center_y) << "/" << (node->size) << std::endl;

	if(node->contains(particles.x[index], particles.y[index], particles.z[index]) && (node->children[0] != nullptr)){
		// std::cout << "La particule " << index <<" est dans le nœud" << std::endl;
		// std::cout << std::endl;
		for (int i = 0; i < 8; i++) {
            compute_force_recursive(node->children[i], index, ax, ay, az, particles, initstate);
        }
	}

	// std::cout << " CoM_x " << node->x << std::endl;
	// std::cout << " CoM_y " << node->y << std::endl;
	// std::cout << " CoM_z " << node->z << std::endl;
	// std::cout << " Masse " << node->mass << std::endl;

    float dx = node->x - particles.x[index];
    float dy = node->y - particles.y[index];
    float dz = node->z - particles.z[index];
    float dist2 = dx * dx + dy * dy + dz * dz;

	// std::cout << " dx " << dx << std::endl;
	// std::cout << " dy " << dy << std::endl;
	// std::cout << " dz " << dz << std::endl;
	// std::cout << " distance carré " << dist2 << std::endl;


    if (dist2 < 1.0) {
        dist2 = 10.0;
    }

    float dist = std::sqrt(dist2);

	// std::cout << " distance " << dist << std::endl;
	// std::cout << " taille " << node->size << std::endl;
	// std::cout << " ratio " << (node->size / dist) << std::endl;

    if ((node->size / dist) < theta || node->children[0] == nullptr) {
        // std::cout << " On approxime au centre de masse " << std::endl;
		float force = 10.0 * (1/(dist * dist * dist));
        ax[index] += node->mass * dx * force;
        ay[index] += node->mass * dy * force;
        az[index] += node->mass * dz * force;
        // std::cout << " force " << node->mass * dist << std::endl;
        // std::cout << " acceleration x/y/z " << ax[index] << "/" << ay[index] << "/" << az[index] << std::endl;
		// std::cout << std::endl;
    } else {
		// std::cout << " On approxime pas au centre de masse " << std::endl;
		// std::cout << std::endl;
        for (int i = 0; i < 8; i++) {
            compute_force_recursive(node->children[i], index, ax, ay, az, particles, initstate);
        }
    }
}

Model_CPU_naive::Model_CPU_naive(const Initstate& initstate, Particles& particles)
    : Model_CPU(initstate, particles) {}

void Model_CPU_naive::step() {

    float max_x = -std::numeric_limits<float>::infinity();
    float max_y = -std::numeric_limits<float>::infinity();
    float max_z = -std::numeric_limits<float>::infinity();
    float min_x = std::numeric_limits<float>::infinity();
    float min_y = std::numeric_limits<float>::infinity();
    float min_z = std::numeric_limits<float>::infinity();

    for (int i = 0; i < n_particles; i++) {
        if(particles.x[i] > max_x) max_x = particles.x[i];
        if(particles.y[i] > max_x) max_y = particles.y[i];
        if(particles.z[i] > max_x) max_z = particles.z[i];
        if(particles.x[i] < min_x) min_x = particles.x[i];
        if(particles.y[i] < min_x) min_y = particles.y[i];
        if(particles.z[i] < min_x) min_z = particles.z[i];
    }

    float dist = std::sqrt(std::pow(max_x-min_x, 2) + std::pow(max_y-min_y, 2) + std::pow(max_z-min_z, 2));
    // std::cout << "dist" << dist << std::endl;
    Octree tree(std::sqrt(std::pow(max_x-min_x, 2) + std::pow(max_y-min_y, 2) + std::pow(max_z-min_z, 2))); // MODIFIER TAILLE

    
    for (int i = 0; i < n_particles; i++) {
        tree.insert_particle(i, particles, initstate);
    }

    // #pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        accelerationsx[i] = 0;
        accelerationsy[i] = 0;
        accelerationsz[i] = 0;
        tree.compute_force(i, accelerationsx, accelerationsy, accelerationsz, particles, initstate);
        // std::cout << "############" << std::endl;
        // std::cout << " acceleration totale x/y/z " << accelerationsx[i] << "/" << accelerationsy[i] << "/" << accelerationsz[i] << std::endl;
        // std::cout << "############" << std::endl;      
    }

    for (int i = 0; i < n_particles; i++) {
        velocitiesx[i] += accelerationsx[i] * 2.0f;
        velocitiesy[i] += accelerationsy[i] * 2.0f;
        velocitiesz[i] += accelerationsz[i] * 2.0f;
        particles.x[i] += velocitiesx[i] * 0.1f;
        particles.y[i] += velocitiesy[i] * 0.1f;
        particles.z[i] += velocitiesz[i] * 0.1f;
    }
}
