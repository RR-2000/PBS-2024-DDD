#pragma once

#include <math.h>
#include <Eigen/Dense>

/// SPH Kernels
/// All of the following methods assume 0 <= r <= h and
class SPHKernels{

public:
    SPHKernels(double kernelRadius);

    double wPoly6(double r2) const;

    Eigen::Vector3d gwPoly6(const Eigen::Vector3d& r) const;

    double lwPoly6(double r2) const;

    double wSpiky(double r) const;

    Eigen::Vector3d gwSpiky(const Eigen::Vector3d& r) const;

    double lwVisc(double r) const;

private:
    double m_h;
    double m_h2;

    double cPoly6;
    double cgPoly6;
    double cSpiky;
    double cgSpiky;
    double clVisc;
};