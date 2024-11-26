#include "GeneratingSource.h"

GeneratingSource::GeneratingSource(Fluid *fluid) : Source(fluid) {
    m_particlesPerSecond = 50;
    m_maxParticles = 3000;

    m_position << 0, 5, 0;
    m_positionStdDeviation = 0.1;
    m_particleVelocity << 0.1, -0.5, 0.1;
    m_particleVelocityStdDeviation = 0;

    m_particleCount = 0;

    std::random_device rd;
    m_gen = std::mt19937(rd());

}

void GeneratingSource::init() {
    m_particleCount = 0;
}

bool GeneratingSource::advance(double t) {
    if (m_maxParticles >= 0 && m_particleCount == m_maxParticles) {
        return false;
    }

    int targetCount = int(floor(t * m_particlesPerSecond));

    int diff = targetCount - m_particleCount;

    for (int i = 0; i < diff; ++i) {
        Particle particle;
        particle.m_mass = m_fluid->m_particleMass;
        particle.m_pos = m_position + random(m_positionStdDeviation);
        particle.m_vel = m_particleVelocity + random(m_particleVelocityStdDeviation);
        particle.m_fluid = m_fluid;
        m_fluid->m_particles.push_back(particle);
        ++m_particleCount;
    }
    return true;
}

Eigen::Vector3d GeneratingSource::random(double variance) {
    std::normal_distribution<float> d(0, variance);
    return {d(m_gen), d(m_gen), d(m_gen)};
}


std::string GeneratingSource::toString() {
    return "Generating Source (" + m_fluid->m_name +")";
}

