#include "BoxScene.h"


BoxScene::BoxScene(const Eigen::Vector3d &m, const Eigen::Vector3d &M) {
    // compute vertices
    m_V.resize(8, 3);
    m_V <<
        m(0), m(1), m(2),
            M(0), m(1), m(2),
            M(0), M(1), m(2),
            m(0), M(1), m(2),
            m(0), m(1), M(2),
            M(0), m(1), M(2),
            M(0), M(1), M(2),
            m(0), M(1), M(2);


    m_E.resize(12, 2);
    m_E <<
        0, 1,
            1, 2,
            2, 3,
            3, 0,
            4, 5,
            5, 6,
            6, 7,
            7, 4,
            0, 4,
            1, 5,
            2, 6,
            7, 3;
}

void BoxScene::getMinMax(Eigen::Vector3d& minCoord, Eigen::Vector3d& maxCoord) {
    minCoord = m_V.colwise().minCoeff();
    maxCoord = m_V.colwise().maxCoeff();
}

void BoxScene::draw(igl::opengl::glfw::Viewer &viewer) {
    for (int i = 0; i < m_E.rows(); i++) {
        viewer.data().add_edges(
                m_V.row(m_E(i, 0)),
                m_V.row(m_E(i, 1)),
                Eigen::RowVector3d(1, 0, 0));
    }
}

bool BoxScene::outOfBoundary(const Eigen::Vector3d &pos) {
    return pos.x() > m_V(1, 0) || pos.x() < m_V(0, 0) ||
           pos.y() < m_V(0, 1) || // only check for floor, leave top "open"
           pos.z() > m_V(4, 2) || pos.z() < m_V(0, 2);
}

void BoxScene::getIntersectionPointAndSurfaceNormal(
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d& dir,
    Eigen::Vector3d& intersection,
    Eigen::Vector3d& normal
) {
    double intersection_dt_x = dir.x() == 0 ? INFINITY
        : dir.x() > 0 ? (m_V(1, 0) - pos.x()) / dir.x()
        : (m_V(0, 0) - pos.x()) / dir.x();
    double intersection_dt_y = dir.y() >= 0 ? INFINITY
        : (m_V(0, 1) - pos.y()) / dir.y();
    double intersection_dt_z = dir.z() == 0 ? INFINITY
        : dir.z() > 0 ? (m_V(4, 2) - pos.z()) / dir.z()
        : (m_V(0, 2) - pos.z()) / dir.z();

    double min_dt = std::min({intersection_dt_x, intersection_dt_y, intersection_dt_z});
    intersection = pos + min_dt * dir;
    normal = Eigen::Vector3d(
        (intersection_dt_x == min_dt) * (1 - 2 * (dir.x() > 0)),
        (intersection_dt_y == min_dt) * (1 - 2 * (dir.y() > 0)),
        (intersection_dt_z == min_dt) * (1 - 2 * (dir.z() > 0))
    ).normalized();
}
