#ifndef PBS_NEIGHBOR_SEARCH_H
#define PBS_NEIGHBOR_SEARCH_H

#include <vector>
#include <map>
#include "Particle.h"

/// Base class for all neighborhood searches
class NeighborSearch {

public:

    virtual ~NeighborSearch() {}

    /// Get neighbors of a particle within a given radius
    virtual std::vector<Particle*> getNeighbors(Particle* particle, double radius) const = 0;

    /// Get particles that are within a given radius of a point pos
    virtual std::vector<Particle*> getNeighbors(Eigen::Vector3d& pos, double radius) const = 0;

    /// Reset underlying data structure
    virtual void reset() = 0;

    /// Add the given particle to the underlying data structure
    virtual void addParticle(Particle* particle) = 0;
};


#endif //PBS_NEIGHBOR_SEARCH_H