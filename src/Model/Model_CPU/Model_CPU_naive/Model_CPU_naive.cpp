#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>
#include <vector>
#include <memory>
#include "Model_CPU_fast.hpp"
#include <omp.h>

// Structure d'un nœud de l'Octree
struct OctreeNode {
    float x, y, z;  // Centre de masse
    float mass;      // Masse totale du nœud
    float size;      // Taille de la boîte englobante
    int particle_index; // Indice de la particule si c'est une feuille (-1 sinon)
    std::vector<std::unique_ptr<OctreeNode>> children;
    
    OctreeNode(float _x, float _y, float _z, float _size)
        : x(_x), y(_y), z(_z), size(_size), mass(0), particle_index(-1) {
        children.resize(8, nullptr);
    }
};

// Classe Octree
class Octree {
public:
    std::unique_ptr<OctreeNode> root;
    float theta = 0.5f; // Paramètre de précision
    
    Octree(float size) {
        root = std::make_unique<OctreeNode>(0.0f, 0.0f, 0.0f, size);
    }

    void insert_particle(int index, const Particles& particles, const std::vector<float>& masses);
    void compute_force(int index, float& ax, float& ay, float& az, const Particles& particles, const std::vector<float>& masses);

private:
    void insert_recursive(OctreeNode* node, int index, const Particles& particles, const std::vector<float>& masses);
    void compute_force_recursive(OctreeNode* node, int index, float& ax, float& ay, float& az, const Particles& particles, const std::vector<float>& masses);
};

// Insère une particule dans l'Octree
void Octree::insert_particle(int index, const Particles& particles, const std::vector<float>& masses) {
    insert_recursive(root.get(), index, particles, masses);
}

void Octree::insert_recursive(OctreeNode* node, int index, const Particles& particles, const std::vector<float>& masses) {
    float px = particles.x[index];
    float py = particles.y[index];
    float pz = particles.z[index];

    if (node->particle_index == -1 && node->children[0] == nullptr) {
        node->particle_index = index;
        node->x = px;
        node->y = py;
        node->z = pz;
        node->mass = masses[index];
    } else {
        if (node->children[0] == nullptr) {
            float halfSize = node->size / 2.0f;
            for (int i = 0; i < 8; i++) {
                float offsetX = (i & 1) ? halfSize : -halfSize;
                float offsetY = (i & 2) ? halfSize : -halfSize;
                float offsetZ = (i & 4) ? halfSize : -halfSize;
                node->children[i] = std::make_unique<OctreeNode>(node->x + offsetX, node->y + offsetY, node->z + offsetZ, halfSize);
            }
            int oldIndex = node->particle_index;
            node->particle_index = -1;
            insert_recursive(root.get(), oldIndex, particles, masses);
        }
        int childIndex = (px > node->x) | ((py > node->y) << 1) | ((pz > node->z) << 2);
        insert_recursive(node->children[childIndex].get(), index, particles, masses);
        node->mass += masses[index];
        node->x = (node->x * node->mass + px * masses[index]) / (node->mass + masses[index]);
        node->y = (node->y * node->mass + py * masses[index]) / (node->mass + masses[index]);
        node->z = (node->z * node->mass + pz * masses[index]) / (node->mass + masses[index]);
    }
}

void Octree::compute_force(int index, float& ax, float& ay, float& az, const Particles& particles, const std::vector<float>& masses) {
    compute_force_recursive(root.get(), index, ax, ay, az, particles, masses);
}

void Octree::compute_force_recursive(OctreeNode* node, int index, float& ax, float& ay, float& az, const Particles& particles, const std::vector<float>& masses) {
    if (!node || (node->particle_index == index && node->particle_index != -1)) {
        return;
    }

    float dx = node->x - particles.x[index];
    float dy = node->y - particles.y[index];
    float dz = node->z - particles.z[index];
    float dist2 = dx * dx + dy * dy + dz * dz;
    
    if (dist2 < 1e-4) {
        return;
    }

    float dist = std::sqrt(dist2);
    if (node->size / dist < theta || node->children[0] == nullptr) {
        float force = node->mass / (dist2 * dist);
        ax += dx * force;
        ay += dy * force;
        az += dz * force;
    } else {
        for (const auto& child : node->children) {
            compute_force_recursive(child.get(), index, ax, ay, az, particles, masses);
        }
    }
}

Model_CPU_fast::Model_CPU_fast(const Initstate& initstate, Particles& particles)
    : Model_CPU(initstate, particles) {}

void Model_CPU_fast::step() {
    Octree tree(100.0f); // MODIFIER TAILLE
    
    for (int i = 0; i < n_particles; i++) {
        tree.insert_particle(i, particles, initstate.masses);
    }

    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        float ax = 0, ay = 0, az = 0;
        tree.compute_force(i, ax, ay, az, particles, initstate.masses);
        accelerationsx[i] = ax;
        accelerationsy[i] = ay;
        accelerationsz[i] = az;
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

#endif // GALAX_MODEL_CPU_FAST