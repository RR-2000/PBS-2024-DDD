#ifndef PBS_EMITTING_SOURCE_H
#define PBS_EMITTING_SOURCE_H

#include "Source.h"
#include <Eigen/Dense>
#include <random>

/// Continuously generates blocks of particles
class EmittingSource : public Source {

public:
    Eigen::Vector3d m_position;
    Eigen::Vector3d m_particleVelocity;

    int m_maxParticles;
    double m_particlesPerSecond;
    int m_x, m_y, m_z;
    double m_spacing;


    explicit EmittingSource(Fluid *fluid, const Eigen::Vector3d& position=Eigen::Vector3d(0.,5.,0.),
        const Eigen::Vector3d& velocity = Eigen::Vector3d(0.1, -0.5, 0.1));

    void init() override;

    bool advance(double t) override;

    std::string toString() override;

private:
    int m_particleCount;

};


#endif //PBS_EMITTING_SOURCE_H