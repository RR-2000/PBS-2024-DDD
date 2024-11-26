#include <igl/writeOFF.h>
#include <thread>
#include "Gui.h"
#include "SphSimulation.h"

#include "BlockSource.h"
#include "EmittingSource.h"
#include "GeneratingSource.h"
#include "CustomSource.h"
#include "FluidDefinitons.h"
#include "BoxScene.h"
#include "SimulationLoader.h"
#include "Simulation.h"


/*
 * GUI for the spring simulation. This time we need additional paramters,
 * e.g. which integrator to use for the simulation and the force applied to the
 * canonball, and we also add some more visualizations (trajectories).
 */
class FluidGui : public Gui {
public:
	// Pointer to the current simulation
    SphSimulation *simulation = nullptr;

	int m_fluid_chooser;

	vector<string> m_fluid_names;

	// Scene boundaries (glass and overall simulation)
    Eigen::Vector3d m_scene_max; // Maximum scene boundaries
    Eigen::Vector3d m_scene_min; // Minimum scene boundaries

    // Paths for external resources
    string m_particles_init_file; // File path for initializing particles (e.g., duck)
    string m_boundary_particles_path; // File path for glass boundary particles

    // New variables for duck object and pouring setup
    Eigen::Vector3d m_duck_position; // Position of the duck object
    Eigen::Vector3d m_duck_velocity; // Velocity of the duck (for movement or dynamics)

    Eigen::Vector3d m_liquid_source_position; // Position of the pouring liquid source
    Eigen::Vector3d m_liquid_source_velocity; // Initial velocity of the pouring liquid




	FluidGui() {
		// Entry point

		// Setup 
		simulation = new SphSimulation();  // Or use a simple SPH solver
		simulation->init();
		setSimulation(simulation);


		// Fluid selection
        m_fluid_chooser = 0;
        for (auto &fluid : fluids::all) {
            m_fluid_names.push_back(fluid->m_name); // Populate fluid names
        }

        // Duck 
        m_duck_position = Eigen::Vector3d(0.0, 1.5, 0.0); // Default position
        m_duck_velocity = Eigen::Vector3d(0.0, 0.0, 0.0); // Default velocity

        // Liquid source setup
        m_liquid_source_position = Eigen::Vector3d(0.0, 3.0, 0.0); // Position above the glass
        m_liquid_source_velocity = Eigen::Vector3d(0.0, -1.0, 0.0); // Downward pouring velocity

		// Add glass boundary (as static particles)
		m_boundary_particles_path = "/home/scai/pbs24/cocktail-mixing/data/cocktail_glass.xyz";  // Precomputed glass particles
		simulation->m_sources.push_back(new CustomSource(fluids::boundary, m_boundary_particles_path));
		simulation->m_sources.back()->init();

		// Add liquid source (pouring water)
		simulation->m_sources.push_back(new EmittingSource(fluids::water, 
			Eigen::Vector3d(0, 3, 0),  // Position above the glass
			Eigen::Vector3d(0, -1, 0) // Velocity (downward)
		));
		simulation->m_sources.back()->init();

		std::cout << "Liquid source added at position: "
          << m_liquid_source_position.transpose()
          << " with velocity: "
          << m_liquid_source_velocity.transpose()
          << std::endl;




		// Add duck object (as dynamic boundary)
		m_particles_init_file = "../../data/duck_object.xyz";  // Need to put the Precomputed duck particles
		simulation->m_sources.push_back(new CustomSource(fluids::boundary, m_particles_init_file));
		simulation->m_sources.back()->init();

		// Setup bounding box for scene
		m_scene_max << 2., 5., 2.;  // Scene limits
		m_scene_min << -2., 0., -2.;
		simulation->setScene(new BoxScene(m_scene_min, m_scene_max));

		// Start GUI
		start();
	}


	void drawSimulationParameterMenu() override {

	// Simulation Settings



    ImGui::InputDouble("Simulation dt", &simulation->m_dt);
    ImGui::InputDouble("Kernel Radius", &simulation->m_kernelRadius);
    ImGui::InputDouble("Grid Width", &simulation->m_gridWidth);
    ImGui::InputDouble("Boundary Repulsion", &simulation->m_boundary_repulsion);

	// position and velocity of the pouring liquid source
	ImGui::InputDouble("Source Pos X", &m_liquid_source_position[0]);
	ImGui::InputDouble("Source Pos Y", &m_liquid_source_position[1]);
	ImGui::InputDouble("Source Pos Z", &m_liquid_source_position[2]);
	ImGui::InputDouble("Source Vel X", &m_liquid_source_velocity[0]);
	ImGui::InputDouble("Source Vel Y", &m_liquid_source_velocity[1]);
	ImGui::InputDouble("Source Vel Z", &m_liquid_source_velocity[2]);



	// Coloring
    ImGui::Checkbox("Particle/Fluid coloring", &simulation->m_use_particle_color);



	// Custom Recording Options
    if (ImGui::CollapsingHeader("Custom Recording")) {
        if (ImGui::Button("Save Sim", ImVec2(-1, 0))) {
            SimulationLoader::saveSimulation("./duck.sim", simulation);
        }
        if (ImGui::Button("Load Sim", ImVec2(-1, 0))) {
            simulation = SimulationLoader::loadSimulation("./duck.sim");
            setSimulation(simulation);
        }
        if (ImGui::Button("Save Particles", ImVec2(-1, 0))) {
            simulation->exportParticles();
        }
        if (ImGui::Button("Save Mesh", ImVec2(-1, 0))) {
            simulation->exportMesh();
        }
    }


	// Boundary Box Configuration
    if (ImGui::CollapsingHeader("Boundary Box")) {
        ImGui::InputDouble("Box min x", &(m_scene_min[0]));
        ImGui::InputDouble("Box min y", &(m_scene_min[1]));
        ImGui::InputDouble("Box min z", &(m_scene_min[2]));
        ImGui::InputDouble("Box max x", &(m_scene_max[0]));
        ImGui::InputDouble("Box max y", &(m_scene_max[1]));
        ImGui::InputDouble("Box max z", &(m_scene_max[2]));

        if (ImGui::Button("Reset Boundary", ImVec2(-1, 0))) {
            simulation->setScene(new BoxScene(m_scene_min, m_scene_max));
        }
    }


	// Glass and Duck Object
    if (ImGui::CollapsingHeader("Glass and Duck Setup")) {
        // Load Glass Boundary
        if (ImGui::Button("Load Glass Boundary", ImVec2(-1, 0))) {
            m_boundary_particles_path = igl::file_dialog_open();
            simulation->m_sources.push_back(new CustomSource(fluids::boundary, m_boundary_particles_path));
            simulation->m_sources.back()->init();
        }

        // Load Duck Object
        if (ImGui::Button("Load Duck Object", ImVec2(-1, 0))) {
            m_particles_init_file = igl::file_dialog_open();
            simulation->m_sources.push_back(new CustomSource(fluids::boundary, m_particles_init_file));
            simulation->m_sources.back()->init();
        }
    }


	// Liquid Source
	if (ImGui::CollapsingHeader("Liquid Source")) {
		// Convert std::vector<std::string> to std::vector<const char*>
		std::vector<const char*> c_strs;
		for (const auto& name : m_fluid_names) {
			c_strs.push_back(name.c_str());
		}

		ImGui::Combo("Fluid Type", &m_fluid_chooser, c_strs.data(), c_strs.size());

		if (ImGui::Button("Add Pouring Liquid Source", ImVec2(-1, 0))) {
			simulation->m_sources.push_back(new EmittingSource(fluids::all[m_fluid_chooser],
															Eigen::Vector3d(0, 3, 0),  // Position
															Eigen::Vector3d(0, -1, 0) // Velocity
			));
			simulation->m_sources.back()->init();
		}
	}


	 // Fluid Properties
    if (ImGui::CollapsingHeader("Fluids")) {
        for (auto &fluid : fluids::all) {
            if (ImGui::CollapsingHeader(fluid->m_name.c_str())) {
                ImGui::InputDouble("Stiffness", &fluid->m_stiffness);
                ImGui::InputDouble("Viscosity", &fluid->m_viscosity);
                ImGui::InputDouble("Particle Mass", &fluid->m_particleMass);
                ImGui::InputDouble("Rest Density", &fluid->m_restDensity);
                ImGui::InputDouble("Surface Tension", &fluid->m_tension);
                ImGui::InputDouble("Surface Tension Threshold", &fluid->m_tension_thres);
            }
        }
    }
	}

	void updateSimulationParameters() override {

		};
			


};

int main(int argc, char* argv[]) {
	// create a new instance of the GUI for the spring simulation
	new FluidGui();

	return 0;
}