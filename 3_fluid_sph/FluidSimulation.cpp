#include "FluidSimulation.h"
#include "FluidDefinitons.h"
#include "BoxScene.h"
#include <thread>
#include <fstream>

using namespace std;

FluidSimulation::FluidSimulation() : Simulation() {
    m_fluids = fluids::all;
    m_sources = vector<Source*>();
    m_use_particle_color = false;
    m_save_simulation = false;
    m_save_freq = 100;
    m_surface_extractor = new SurfaceExtractor();
    m_mesh_path = "../../mesh";
    m_particles_path = "../../particles";
    m_scene = new BoxScene(Eigen::Vector3d(0,0,0), Eigen::Vector3d(1,1,1));
}

void FluidSimulation::init() {
    reset();
}

void FluidSimulation::resetMembers() {
    for (auto &fluid : m_fluids) {
        fluid->reset();
    }
    for (auto &source : m_sources) {
        source->init();
    }
}

bool FluidSimulation::advance() {
    m_time += m_dt;
    for (auto &source : m_sources) {
        source->advance(m_time);
    }
    m_step++;

    // save mesh
    if(m_save_simulation && (m_step % m_save_freq == 0)) {
        exportMesh();
        exportParticles();
    }

    return false;
}

void FluidSimulation::updateRenderGeometry() {
    int numParticles = 0;
    for (auto &fluid : m_fluids) {
        numParticles += fluid->m_particles.size();
    }

    V.setZero(numParticles, 3);
    C.setZero(numParticles, 3);

    int i = 0;
    if(m_use_particle_color){
        for (auto &fluid : m_fluids) {
            for (auto &particle : fluid->m_particles) {
                V.row(i) = particle.m_pos.transpose();
                C.row(i) = particle.m_color.transpose();
                ++i;
            }
        }
    }else{
        // use fluid color
        for (auto &fluid : m_fluids) {
            for (auto &particle : fluid->m_particles) {
                V.row(i) = particle.m_pos.transpose();
                C.row(i) = fluid->m_color.transpose();
                ++i;
            }
        }
    }
}

void FluidSimulation::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) {
    viewer.data().point_size = 4;
    viewer.data().set_points(V, C);
}

void FluidSimulation::runParallel(int elementCount, const std::function<void(int, int)> &f) {
    int numThreads = std::max(4, (int) std::thread::hardware_concurrency());
    int elementsPerThread = ceil(1. * elementCount / numThreads);
    std::vector<std::thread> threads(numThreads);

    if (elementCount == 0) {
        return;
    }

    for (int i = 0; i < numThreads; ++i) {
        int start = std::min(elementCount - 1, i * elementsPerThread);
        int end = std::min(elementCount, (i + 1) * elementsPerThread);
        threads[i] = std::thread(f, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }
}


void FluidSimulation::getMinMaxParticlePosition(Eigen::Vector3d& minPosition, Eigen::Vector3d& maxPosition) {
    minPosition = Eigen::Vector3d(INFINITY, INFINITY, INFINITY);
    maxPosition = Eigen::Vector3d(-INFINITY, -INFINITY, -INFINITY);
    for (auto &fluid : m_fluids) {
        for (auto &particle : fluid->m_particles) {
            for(int i = 0; i < 3; i++) {
                minPosition[i] = std::min(minPosition[i], particle.m_pos[i]);
                maxPosition[i] = std::max(maxPosition[i], particle.m_pos[i]);
            }
        }
    }
}

void FluidSimulation::setScene(Scene* scene) {
    if(m_scene != NULL) delete m_scene;
    m_scene = NULL;
    m_scene = scene;
    // adjust bounding box for marching cubes
    scene->getMinMax(m_surface_extractor->m_m, m_surface_extractor->m_M);
}

void FluidSimulation::toggleRecording() {
    m_save_simulation = !m_save_simulation;
    if(m_save_simulation) {
        std::cout << "Start recording\n";
    }else{
        std::cout << "Stop recording\n";
    }
}

void FluidSimulation::exportParticles() {
    for(auto &fluid : m_fluids) {
        // no need to save save non existant or fixed particles
        if(fluid->m_particles.size() > 0 && fluid->m_name != "Boundary") {
            ofstream particleFile;
            string filepath = m_particles_path + "/" + fluid->m_name + "_" + std::to_string(m_step) + ".xyz";
            particleFile.open(filepath);
            for(auto& particle : fluid->m_particles) {
                particleFile << particle.m_pos.transpose() << "\n";
            }
            particleFile.close();
            std::cout<< "Written particles to " << filepath << "\n";
        }
    }
}

void FluidSimulation::exportMesh() {
    for(auto &fluid : m_fluids) {
        if(fluid->m_name == "Boundary") continue;
        std::string filePath = m_mesh_path + "/" + fluid->m_name + "_" + std::to_string(m_step) + ".obj";
        m_surface_extractor->createMesh(fluid->m_particles, filePath);
    }
}