#ifndef PBS_CONSTANTS_H
#define PBS_CONSTANTS_H

#include <Eigen/Dense>

namespace constants {

    const Eigen::Vector3d g(0, -9.81 * 0.1, 0);

    constexpr double taitEq(double k, double rest_density, double density) {
        return k * rest_density * (pow(density/rest_density, 7) - 1) / 7.0;
    }

};



#endif //PBS_CONSTANTS_H