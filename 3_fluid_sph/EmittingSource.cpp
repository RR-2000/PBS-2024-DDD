#include "EmittingSource.h"
#include <iostream>

EmittingSource::EmittingSource(Fluid *fluid, const Eigen::Vector3d& position, const Eigen::Vector3d& velocity)
     : Source(fluid), m_position(position), m_particleVelocity(velocity) {

    m_particlesPerSecond = 20;
    m_maxParticles = 3000;

    m_x = 4;
    m_y = 4;
    m_z = 1;
    m_spacing = 1.1 * std::cbrt(fluid->m_particleMass / fluid->m_restDensity);

    m_particleCount = 0;
}

void EmittingSource::init() {
    m_particleCount = 0;
}

bool EmittingSource::advance(double t) {

    if (m_maxParticles >= 0 && m_particleCount >= m_maxParticles) {
        return false;
    }

    int targetCount = int(floor(t * m_particlesPerSecond));

    int diff = targetCount - m_particleCount;

    for (int i = 0; i < diff; ++i) {
        Particle particle;
        particle.m_mass = m_fluid->m_particleMass;
        particle.m_vel = m_particleVelocity;
        particle.m_fluid = m_fluid;
        for(int x = 0; x < m_x; x++) {
            for(int y = 0; y < m_y; y++) {
                for(int z = 0; z < m_z; z++) {
                    particle.m_pos = m_position + Eigen::Vector3d(x*m_spacing, y*m_spacing, z*m_spacing);
                    m_fluid->m_particles.push_back(particle);
                    ++m_particleCount;
                }
            }
        }

    }
    return true;
}

std::string EmittingSource::toString() {
    return "Emitting Source (" + m_fluid->m_name +")";
}

