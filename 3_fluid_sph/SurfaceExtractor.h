#ifndef SURFACE_EXTRACTOR
#define SURFACE_EXTRACTOR

#include <vector>
#include <string>
#include <Eigen/Dense>

#include "SPHKernels.h"

class NeighborSearch;
class Particle;

/// Given particle locations from an SPH simulation, calculate the densities and extract an iso surface
class SurfaceExtractor {
public:

	SurfaceExtractor();

	/* Given particle locations from an SPH simulation, calculate the densities and extract an isosurface
	   positions: particle positions of the same fluid
	   file_path: where to save as obj file
	*/
	void createMesh(std::vector<Particle>& particles, std::string file_path);

	double m_res;			// grid width for marching cubes voxel grid
	double m_kernel_radius; // kernel radius for density computation
	double m_isolevel;		// density level for mesh extraction

	Eigen::Vector3d m_m; // min coordinate of the bounding box
	Eigen::Vector3d m_M; // max coordinate of the bounding box

	SPHKernels* m_kernels;
private:
	double computeDensity(const Eigen::Vector3d& pos, const std::vector<Particle*>& neighbors);
	double computeSDF(const Eigen::Vector3d& pos, const std::vector<Particle*>& neighbors);

};


#endif // SURFACE_EXTRACTOR