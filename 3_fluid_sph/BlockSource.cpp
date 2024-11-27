#include "BlockSource.h"

#include <utility>

BlockSource::BlockSource(Fluid *fluid)
        : Source(fluid),
          m_initialDimension(10, 10, 10),
          m_initialSpacing(.1),
          m_initialOffset(1.5, 2.5, 1.5),
          m_initialVelocity(0.0, 0.0, 0.0) {}

BlockSource::BlockSource(
        Fluid *fluid,
        Eigen::Vector3i initialDimension,
        Eigen::Vector3d initialOffset,
        Eigen::Vector3d initialVelocity)
        : Source(fluid),
          m_initialDimension(std::move(initialDimension)),
          m_initialOffset(std::move(initialOffset)),
          m_initialVelocity(std::move(initialVelocity)),
          m_initialSpacing(std::cbrt(fluid->m_particleMass / fluid->m_restDensity) * 1.01) {}


void BlockSource::init() {
    auto &d = m_initialDimension;

    for (int x = 0; x < d.x(); ++x) {
        for (int y = 0; y < d.y(); ++y) {
            for (int z = 0; z < d.z(); ++z) {
                Particle particle;
                particle.m_mass = m_fluid->m_particleMass;
                particle.m_pos << x, y, z;
                particle.m_pos *= m_initialSpacing;
                particle.m_pos += m_initialOffset;
                particle.m_vel = m_initialVelocity;
                particle.m_color << 1. * x / d.x(), 1. * y / d.y(), 1. * z / d.z();
                particle.m_fluid = m_fluid;
                m_fluid->m_particles.push_back(particle);
            }
        }
    }
}

std::string BlockSource::toString() {
    return "Block Source (" + m_fluid->m_name + ")";
}