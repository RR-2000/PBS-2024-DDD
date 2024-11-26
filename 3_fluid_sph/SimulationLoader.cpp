#include "SimulationLoader.h"
#include "BoxScene.h"
#include "SphSimulation.h"
#include "Fluid.h"
#include "FluidDefinitons.h"
#include "StateSource.h"

#include <fstream>

using namespace std;
// Static function to save the simulation
void SimulationLoader::saveSimulation(std::string exportPath, SphSimulation* sim) {
    ofstream exportFile;
    exportFile.open(exportPath, ios::out | ios::trunc);

    // Save SPH parameters
    exportFile << "sph\n";
    exportFile << sim->m_dt << " " << sim->m_gridWidth << " " << sim->m_kernelRadius << "\n";

    // Save scene information
    auto box_cast = dynamic_cast<BoxScene*>(sim->m_scene);
    if (box_cast != nullptr) {
        Eigen::Vector3d m, M;
        sim->m_scene->getMinMax(m, M);
        exportFile << "box " << m.transpose() << " " << M.transpose() << "\n";
    } else {
        exportFile << "no_scene\n";
    }

    exportFile << sim->m_fluids.size() << "\n";

    // Save particles
    for (auto& fluid : sim->m_fluids) {
        // Save fluid information
        exportFile << fluid->m_name << " " << fluid->m_particles.size() << "\n";
        // Save particle positions and velocities
        for (auto& particle : fluid->m_particles) {
            exportFile << particle.m_pos.transpose() << " " << particle.m_vel.transpose() << "\n";
        }
    }

    exportFile.close();
}

// Static function to load the simulation
SphSimulation* SimulationLoader::loadSimulation(std::string loadPath) {
    SphSimulation* sim = new SphSimulation();

    ifstream loadFile(loadPath);
    if (!loadFile.is_open()) {
        cerr << "Failed to open file: " << loadPath << endl;
        return nullptr;
    }

    // Read SPH parameters
    std::string reader;
    loadFile >> reader;

    if (reader != "sph") {
        cerr << "Invalid simulation file format." << endl;
        loadFile.close();
        delete sim;
        return nullptr;
    }

    // Clear all members and reset
    sim->resetMembers();

    // Read SPH parameters
    loadFile >> sim->m_dt;
    loadFile >> sim->m_gridWidth;
    loadFile >> sim->m_kernelRadius;

    // Set scene information
    loadFile >> reader;

    if (reader == "box") {
        Eigen::Vector3d m, M;
        double x, y, z;

        loadFile >> x >> y >> z;
        m = Eigen::Vector3d(x, y, z);

        loadFile >> x >> y >> z;
        M = Eigen::Vector3d(x, y, z);

        sim->m_scene = new BoxScene(m, M);
    }

    // Restore particles
    int num_fluids;
    loadFile >> num_fluids;
    string fluid_name;
    int num_particles;
    for (int i = 0; i < num_fluids; i++) {
        loadFile >> fluid_name >> num_particles;
        Fluid* fluid = nullptr;

        for (auto temp_fluid : sim->m_fluids) {
            if (fluid_name == temp_fluid->m_name) {
                fluid = temp_fluid;
                break;
            }
        }

        if (fluid == nullptr) {
            cerr << "Unknown fluid type: " << fluid_name << endl;
            continue;
        }

        std::vector<Particle> particles;
        for (int j = 0; j < num_particles; j++) {
            double x, y, z;
            Particle particle;
            particle.m_mass = fluid->m_particleMass;
            loadFile >> x >> y >> z;
            particle.m_pos = Eigen::Vector3d(x, y, z);
            loadFile >> x >> y >> z;
            particle.m_vel = Eigen::Vector3d(x, y, z);
            particle.m_fluid = fluid;
            particles.push_back(particle);
        }

        sim->m_sources.push_back(new StateSource(fluid, particles));
    }

    loadFile.close();

    return sim;
}