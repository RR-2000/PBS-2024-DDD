#include "CustomSource.h"
#include "Particle.h"
#include <fstream>
#include <iostream>

CustomSource::CustomSource(Fluid* fluid, std::string& file_path) : Source(fluid) {
    m_file_path = file_path;
}

void CustomSource::init() {
    // try to read file
    std::ifstream file(m_file_path);
    if(!file) {
        std::cout << "Failed to open file " << m_file_path << std::endl;
    }
    double x, y, z;
    while(file >> x >> y >>z) {
        Particle particle;
        particle.m_mass = m_fluid->m_particleMass;
        particle.m_pos << x, y, z;
        particle.m_fluid = m_fluid;
        m_fluid->m_particles.push_back(particle);
    }
    std::cout << "Finished reading particles \n";
    file.close();
}