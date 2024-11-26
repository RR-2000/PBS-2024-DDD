#ifndef PBS_FLUID_H
#define PBS_FLUID_H

#include <vector>
#include "Particle.h"

/// Holds all fluid properties and contains all particles for this fluid
class Fluid {

public:

    /// All particles of the fluid
    std::vector<Particle> m_particles;

    /// Human readable name of the fluid
    std::string m_name;

    /* SPH parameters */
    bool m_isBoundary;
    double m_viscosity;
    double m_stiffness;
    double m_restDensity;
    double m_particleMass;
    double m_tension;
    double m_tension_thres;
    Eigen::Vector3d m_color;

    Fluid(std::string name, double viscosity, double stiffness, double restDensity,
          double particleMass, Eigen::Vector3d color, bool isBoundary, double tension = 0.1, double m_tension_thres = 1e-6);

    /// Clear particles
    void reset();
};


#endif // PBS_FLUID_H