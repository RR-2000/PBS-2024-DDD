#ifndef PBS_SIMULATION_LOADER_H
#define PBS_SIMULATION_LOADER_H

#include <string>

class SphSimulation;

/// Class for saving and loading SPH simulation states
class SimulationLoader{
public:
    static void saveSimulation(std::string exportPath, SphSimulation* sim);
    static SphSimulation* loadSimulation(std::string loadPath);
};

#endif