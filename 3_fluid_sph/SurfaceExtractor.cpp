#include "SurfaceExtractor.h"
#include "UniformGridNeighborSearch.h"
#include "Particle.h"

#include <igl/writeOBJ.h>
#include <igl/copyleft/marching_cubes.h>

SurfaceExtractor::SurfaceExtractor() {
    m_res = 0.01;
    m_kernel_radius = 0.2;
    m_isolevel = 0.7;
    m_m = Eigen::Vector3d(-1, 0, -1);
    m_M = Eigen::Vector3d(2, 4, 2);
    m_kernels = new SPHKernels(m_kernel_radius);
}

void SurfaceExtractor::createMesh(std::vector<Particle>& particles, std::string file_path) {
    // initialize uniform grid
    UniformGridNeighborSearch* neighborSearch = new UniformGridNeighborSearch(m_kernel_radius); // take resolution of the simulation

    if(particles.size() <= 0) {
        return;
    }

    // find bounding box
    Eigen::Vector3d m(particles.begin()->m_pos), M(particles.begin()->m_pos);
    for(auto& particle : particles){
        neighborSearch->addParticle(&particle);
    }

    Eigen::Vector3i gridDim;
    gridDim = ((m_M - m_m).array() / m_res).cast<int>();

    int gs = gridDim.x() * gridDim.y() * gridDim.z();

    Eigen::MatrixXd densities(gs, 1);
    densities.setZero();
    Eigen::MatrixXd points(gs, 3);

    int c = 0;
    for(int x = 0; x < gridDim.x(); x++) {
        Eigen::Vector3d pos;
        pos[0] = m_m.x() + x*m_res;
        for(int y = 0; y < gridDim.y(); y++) {
            pos[1] = m_m.y() + y*m_res;
            for(int z = 0; z  < gridDim.z(); z++) {
                pos[2] = m_m.z() + z*m_res;
                points.row(c) = pos.transpose();

                //calculate density using SPH kernel
                std::vector<Particle*> neighbors;
                neighbors = neighborSearch->getNeighbors(pos, m_kernel_radius);
        
                // densities(c) = computeDensity(pos, neighbors);
                densities(c) = computeSDF(pos, neighbors);

                c++;
            }
        }
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;

    // marching cubes
    igl::copyleft::marching_cubes(densities, points, gridDim.x(), gridDim.y(), gridDim.z(), m_isolevel, vertices, faces);

    // write to obj
    igl::writeOBJ(file_path, vertices, faces);
    std::cout << "Save mesh to " << file_path << '\n';

}

double SurfaceExtractor::computeDensity(const Eigen::Vector3d& pos, const std::vector<Particle*>& neighbors) {
    double density = 0.0;
    for(Particle* neighbor : neighbors) {
        double r2 = (pos - neighbor->m_pos).squaredNorm();
        density += neighbor->m_mass * m_kernels->wPoly6(r2) / neighbor->m_density;
    }
    return density;
}

/* According to paper: Particle-based Viscoelastic Fluid Simulation */
double SurfaceExtractor::computeSDF(const Eigen::Vector3d& pos, const std::vector<Particle*>& neighbors) {
    double density = 0.0;
    for(Particle* neighbor : neighbors) {
        double r = (pos - neighbor->m_pos).norm();
        density += std::pow((1 - r/m_kernel_radius), 2);
    }
    return std::sqrt(density);
}