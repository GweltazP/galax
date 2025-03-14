#ifndef MODEL_CPU_NAIVE_HPP_
#define MODEL_CPU_NAIVE_HPP_

#include "../Model_CPU.hpp"

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

class Model_CPU_naive : public Model_CPU
{
public:
    Model_CPU_naive(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_naive() = default;

    virtual void step();
};
#endif // MODEL_CPU_NAIVE_HPP_