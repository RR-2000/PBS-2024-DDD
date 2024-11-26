#ifndef PBS_BOX_SCENE_H
#define PBS_BOX_SCENE_H

#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include "Scene.h"

/// A simple analytical bounding box scene
class BoxScene : public Scene {
public:
    BoxScene(const Eigen::Vector3d& m, const Eigen::Vector3d& M);

	void draw(igl::opengl::glfw::Viewer &viewer);

	bool outOfBoundary(const Eigen::Vector3d& pos);

	void getMinMax(Eigen::Vector3d& minCoord, Eigen::Vector3d& maxCoord);

	void getIntersectionPointAndSurfaceNormal(
		const Eigen::Vector3d& pos,
		const Eigen::Vector3d& dir,
		Eigen::Vector3d& intersection,
		Eigen::Vector3d& normal
	);

private:
	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_E;
};

#endif