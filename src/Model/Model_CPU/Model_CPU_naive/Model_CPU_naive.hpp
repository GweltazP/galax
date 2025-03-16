#ifndef MODEL_CPU_NAIVE_HPP_
#define MODEL_CPU_NAIVE_HPP_

#include "../Model_CPU.hpp"

// Structure d'un nœud de l'Octree

struct OctreeNode {
    float x, y, z;  // Centre de masse
    float center_x, center_y, center_z; // Centre de la boite
    float mass;      // Masse totale du nœud
    float size;      // Taille de la boîte englobante
    int particle_index; // Indice de la particule si c'est une feuille (-1 sinon)
    OctreeNode* children[8];
    
    OctreeNode(float _x, float _y, float _z, float _size) : x(0), y(0), z(0), center_x(_x), center_y(_y), center_z(_z), size(_size), mass(0), particle_index(-1)
    {
        for (int i = 0; i < 8; i++)
        {
            children[i] = nullptr;
        }
    }

    bool contains(float _x, float _y, float _z)
    {
        return (
            _x >= x - size / 2 && _x <= x + size / 2 &&
            _y >= y - size / 2 && _y <= y + size / 2 &&
            _z >= z - size / 2 && _z <= z + size / 2
        );
    }

    void subdivide()
    {
        float halfSize = size / 2.0f;
        for (int i = 0; i < 8; i++)
        {
            float offsetX = (i & 1) ? halfSize : -halfSize;
            float offsetY = (i & 2) ? halfSize : -halfSize;
            float offsetZ = (i & 4) ? halfSize : -halfSize;
            children[i] = new OctreeNode(center_x + offsetX, center_y + offsetY, center_z + offsetZ, halfSize);
        };
    };
};

// Classe Octree
class Octree {
public:
    OctreeNode* root;
    float theta = 0.5f; // Paramètre de précision
    
    Octree(float size) {
        root = new OctreeNode(0.0f, 0.0f, 0.0f, size);
    }

    ~Octree() {
        delete_tree(root);
    }

    void insert_particle(int index, const Particles& particles, const Initstate& initstate);
    void compute_force(int index, std::vector<float>& ax, std::vector<float>& ay, std::vector<float>& az, const Particles& particles, const Initstate& initstate);

private:
    void insert_recursive(OctreeNode* node, int index, const Particles& particles, const Initstate& initstate);
    void compute_force_recursive(OctreeNode* node, int index, std::vector<float>& ax, std::vector<float>& ay, std::vector<float>& az, const Particles& particles, const Initstate& initstate);
    void delete_tree(OctreeNode* node);
};

class Model_CPU_naive : public Model_CPU
{
public:
    Model_CPU_naive(const Initstate& initstate, Particles& particles);

    virtual ~Model_CPU_naive() = default;

    virtual void step();
};
#endif // MODEL_CPU_NAIVE_HPP_