#include "StateSource.h"

#include <iostream>

StateSource::StateSource(Fluid* fluid, std::vector<Particle> particles) : Source(fluid), m_particles(particles) {}

void StateSource::init() {
    m_fluid->m_particles.insert(m_fluid->m_particles.end(), m_particles.begin(), m_particles.end());
}

std::string StateSource::toString() {
    return "State Source (" + m_fluid->m_name +")";
}