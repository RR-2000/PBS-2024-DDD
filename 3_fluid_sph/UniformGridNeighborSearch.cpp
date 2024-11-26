#include "UniformGridNeighborSearch.h"
#include <algorithm>

UniformGridNeighborSearch::UniformGridNeighborSearch(double gridWidth) : gridWidth(gridWidth) {}

std::vector<Particle*> UniformGridNeighborSearch::getNeighbors(Eigen::Vector3d& pos, double radius) const {
    int x, y, z;
    std::tie(x, y, z) = getIndex(pos);

    double epsilon = gridWidth / 10000;
    double radiusSquared = pow(radius, 2);

    std::vector<Particle*> neighbors;
    int delta = int(ceil(radius / gridWidth));
    for (int dx = -delta; dx <= delta; ++dx) {
        auto gridX = uniformGrid.find(x + dx);
        if (gridX == uniformGrid.end()) {
            continue;
        }
        for (int dy = -delta; dy <= delta; ++dy) {
            auto gridXY = gridX->second.find(y + dy);
            if (gridXY == gridX->second.end()) {
                continue;
            }
            for (int dz = -delta; dz <= delta; ++dz) {
                auto gridXYZ = gridXY->second.find(z + dz);
                if (gridXYZ == gridXY->second.end()) {
                    continue;
                }
                for (auto neighbor : gridXYZ->second) {
                    if ((neighbor->m_pos - pos).squaredNorm() <= radiusSquared + epsilon) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }
    }
    return neighbors;
}

std::vector<Particle*> UniformGridNeighborSearch::getNeighbors(Particle* particle, double radius) const {
    particle->m_neighbors.clear();
    particle->m_neighbors = getNeighbors(particle->m_pos, radius);
    return particle->m_neighbors;
}

void UniformGridNeighborSearch::addParticle(Particle *particle) {
    int x, y, z;
    std::tie(x, y, z) = getIndex(particle->m_pos);

    uniformGrid[x][y][z].push_back(particle);
}

std::tuple<int, int, int> UniformGridNeighborSearch::getIndex(Eigen::Vector3d& position) const {
    int x = int(floor(position(0) / gridWidth));
    int y = int(floor(position(1) / gridWidth));
    int z = int(floor(position(2) / gridWidth));
    return {x, y, z};
}

void UniformGridNeighborSearch::reset() {
    uniformGrid.clear();
}


