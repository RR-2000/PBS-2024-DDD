#include "SPHKernels.h"

SPHKernels::SPHKernels(double kernelRadius) {
    m_h = kernelRadius;
    m_h2 = kernelRadius * kernelRadius;

    cPoly6 = 315 / (64 * M_PI  * pow(m_h, 9));
    cgPoly6 = -945 / (32 * M_PI * pow(m_h, 9));
    cSpiky = 15 / (M_PI * pow(m_h, 6));
    cgSpiky = -45 / (M_PI * pow(m_h, 6));
    clVisc = 45 / (M_PI * pow(m_h, 6));
}

double SPHKernels::wPoly6(double r2) const{
    return  cPoly6 * pow(m_h2 - r2, 3);
}

Eigen::Vector3d SPHKernels::gwPoly6(const Eigen::Vector3d& r) const{
    return cgPoly6 * pow(m_h2-r.squaredNorm(), 2) * r;
}

double SPHKernels::lwPoly6(double r2) const{
    return cgPoly6 * (m_h2-r2) * (3*m_h2 - 7*r2);
}

double SPHKernels::wSpiky(double r) const{
    return cSpiky * pow(m_h - r, 3);
}

Eigen::Vector3d SPHKernels::gwSpiky(const Eigen::Vector3d& r) const{
    return cgSpiky * pow(m_h - r.norm(), 2) * r;
}

double SPHKernels::lwVisc(double r) const{
    return clVisc * (m_h - r);
}