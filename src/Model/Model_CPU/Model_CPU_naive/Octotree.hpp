#ifndef OCTOTREE_H
#define OCTOTREE_H

#include <iostream>
#include <vector>

// Structure pour représenter un point en 3D
struct Point {
    float x, y, z;
    float mass;
};

// Structure pour représenter une région cubique
struct Boundary {
    float x, y, z, size;
    float contains(const Point& p) const
    {
        return (p.x >= x - size / 2 && p.x <= x + size / 2 &&
                p.y >= y - size / 2 && p.y <= y + size / 2 &&
                 p.z >= z - size / 2 && p.z <= z + size / 2);
    };
};

// Classe pour l'Octotree
class Octotree {
private:
    static constexpr int CAPACITY = 1;
    Boundary boundary;
    Point* body = nullptr;
    bool divided = false;
    Octotree* children[8] = {nullptr};
    Point centerOfMass = {0, 0, 0, 0};

    void subdivide();

public:
    Octotree(Boundary boundary);
    ~Octotree();

    bool insert(Point p);
};

#endif // OCTOTREE_H
