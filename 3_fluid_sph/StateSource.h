#ifndef PBS_STATE_SOURCE_H
#define PBS_STATE_SOURCE_H

#include "Source.h"
#include "Particle.h"

#include <vector>
#include <Eigen/Dense>

/// Adds vector of particles to scene
class StateSource : public Source {

public:

    StateSource(Fluid *fluid, std::vector<Particle> particles);

    void init() override;

    std::string toString() override;

private:
    std::vector<Particle> m_particles;

};


#endif // PBS_STATE_SOURCE_H