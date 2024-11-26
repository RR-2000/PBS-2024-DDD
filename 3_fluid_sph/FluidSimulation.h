#ifndef PBS_PARTICLE_SIMULATION_H
#define PBS_PARTICLE_SIMULATION_H

#include "Simulation.h"
#include "Particle.h"
#include "Scene.h"
#include "Fluid.h"
#include "Source.h"
#include "SurfaceExtractor.h"


/// Base class for fluid Simulations
/// Provides basic functionalities, like managing particles and displaying them
class FluidSimulation : public Simulation {
public:
    explicit FluidSimulation();

    /// All fluids of the simulation
    std::vector<Fluid*> m_fluids;

    /// All sources of the simulation
    std::vector<Source*> m_sources;

    /// If true, particles are drawn in the color specified by the corresponding fluid
	bool m_use_particle_color;

	/// Selected surface extractor
    SurfaceExtractor* m_surface_extractor;

    /// Output path for the mesh export
    std::string m_mesh_path;

    /// Output path for the particle export
    std::string m_particles_path;

    /// Get bounding box of all particles
	void getMinMaxParticlePosition(Eigen::Vector3d& minPosition, Eigen::Vector3d& maxPosition);

	/// Toggle recording
	void toggleRecording();

	/// Export particles
	void exportParticles();

	/// Export mesh
	void exportMesh();

	/// Set a new scene
	void setScene(Scene* scene);


    /* Overrides */
    virtual void init() override;
    virtual void resetMembers() override;
    virtual bool advance() override;

    void updateRenderGeometry() override;
    void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) override;

protected:

    /// Run the given function in parallel over elementCount elements
	static void runParallel(int elementCount, const std::function<void(int, int)>& f);

    /// Scene object
	Scene* m_scene;

private:
    /// Vertex positions
	Eigen::MatrixXd V;

    /// Faces (indices of vertices)
    Eigen::MatrixXi F;

    /// Colors per face or per vertex
    Eigen::MatrixXd C;

    /// If true, the simulation state will be saved periodically
	bool m_save_simulation;

	/// Frequency of simulation state saves
	int m_save_freq;

};

#endif // PBS_PARTICLE_SIMULATION_H