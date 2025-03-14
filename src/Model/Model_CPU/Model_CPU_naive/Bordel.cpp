// Structure d'un nœud de l'Octree
struct OctreeNode {
    float x, y, z;  // Centre de masse
    float mass;      // Masse totale du nœud
    float size;      // Taille de la boîte englobante
    int particle_index; // Indice de la particule si c'est une feuille (-1 sinon)
    OctreeNode* children[8];
    
    OctreeNode(float _x, float _y, float _z, float _size)
        : x(_x), y(_y), z(_z), size(_size), mass(0), particle_index(-1) {
        for (int i = 0; i < 8; i++) {
            children[i] = nullptr;
        }
    }
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

    void insert_particle(int index, const Particles& particles, const float* masses);
    void compute_force(int index, float& ax, float& ay, float& az, const Particles& particles, const float* masses);

private:
    void insert_recursive(OctreeNode* node, int index, const Particles& particles, const float* masses);
    void compute_force_recursive(OctreeNode* node, int index, float& ax, float& ay, float& az, const Particles& particles, const float* masses);
    void delete_tree(OctreeNode* node);
};

void Octree::delete_tree(OctreeNode* node) {
    if (!node) return;
    for (int i = 0; i < 8; i++) {
        delete_tree(node->children[i]);
    }
    delete node;
}

void Octree::insert_particle(int index, const Particles& particles, const float* masses) {
    insert_recursive(root, index, particles, masses);
}

void Octree::insert_recursive(OctreeNode* node, int index, const Particles& particles, const float* masses) {
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
                node->children[i] = new OctreeNode(node->x + offsetX, node->y + offsetY, node->z + offsetZ, halfSize);
            }
            int oldIndex = node->particle_index;
            node->particle_index = -1;
            insert_recursive(root, oldIndex, particles, masses);
        }
        int childIndex = (px > node->x) | ((py > node->y) << 1) | ((pz > node->z) << 2);
        insert_recursive(node->children[childIndex], index, particles, masses);
        node->mass += masses[index];
        node->x = (node->x * node->mass + px * masses[index]) / (node->mass + masses[index]);
        node->y = (node->y * node->mass + py * masses[index]) / (node->mass + masses[index]);
        node->z = (node->z * node->mass + pz * masses[index]) / (node->mass + masses[index]);
    }
}

void Octree::compute_force(int index, float& ax, float& ay, float& az, const Particles& particles, const float* masses) {
    compute_force_recursive(root, index, ax, ay, az, particles, masses);
}

void Octree::compute_force_recursive(OctreeNode* node, int index, float& ax, float& ay, float& az, const Particles& particles, const float* masses) {
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
        for (int i = 0; i < 8; i++) {
            compute_force_recursive(node->children[i], index, ax, ay, az, particles, masses);
        }
    }
}