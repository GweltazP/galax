#include <cmath>

#include "Octotree.hpp"

void Octotree::subdivide()
{
    float new_size = self.boundary.size / 2;
    for(i = 0; i < 8; i++)
    {
        float newX = self.boundary.x - (-1)^(i) * (new_size / 2);
        float newY = self.boundary.y - (-1)^(i/4) * (new_size / 2);
        float newZ = self.boundary.z - (-1)^(i/2) * (new_size / 2);
        self.children[i] = new Octotree({newX, newY, newZ, new_size});
    }
    self.divided = true;
};

bool Octotree::insert(Point p)
{
    if(!self.boundary.contains(p))
    {
        return false;
    }
    if(self.body == nullptr)
    {
        self.body = new Point(p);
        self.centerOfMass = p;
        return true;
    }
    if(!self.divided)
    {
        self.subdivide();
    }

    float newCoM_X = (self.centerOfMass.x + p.x) / 2;
    float newCoM_Y = (self.centerOfMass.y + p.y) / 2;
    float newCoM_Z = (self.centerOfMass.z + p.z) / 2;
    float newMass = self.centerOfMass.mass + p.mass;

    self.centerOfMass = new Point({newCoM_X, newCoM_Y, newCoM_Z, newMass});

    for(i = 0; i < 8; i++)
    {
        self.children[i].insert(p);
        return true;
    }
};

void Octotree::~Octotree()
{
    for (int i = 0; i < 8; i++)
    {
        delete children[i];
    }
}