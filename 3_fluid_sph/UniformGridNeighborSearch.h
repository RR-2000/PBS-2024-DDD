#ifndef PBS_UNIFORM_GRID_NEIGHBOR_SEARCH_H
#define PBS_UNIFORM_GRID_NEIGHBOR_SEARCH_H

#include <map>
#include <tuple>
#include "NeighborSearch.h"

/// Efficient neighborhood search implementation
class UniformGridNeighborSearch : public NeighborSearch {

public:
    UniformGridNeighborSearch(double gridWidth);

    void reset() override;

    std::vector<Particle*> getNeighbors(Particle* particle, double radius) const override;

    std::vector<Particle*> getNeighbors(Eigen::Vector3d& pos, double radius) const override;

    void addParticle(Particle *particle) override;

private:
    double gridWidth;

    std::map<int, std::map<int, std::map<int, std::vector<Particle*>>>> uniformGrid;

    std::tuple<int, int, int> getIndex(Eigen::Vector3d& position) const;
};


#endif //PBS_UNIFORM_GRID_NEIGHBOR_SEARCH_H