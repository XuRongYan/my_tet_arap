#include <iostream>
#include <cinolib/meshes/meshes.h>
#include <vector>
#include "arap_tet.h"

int main() {
    cinolib::Tetmesh<> m("shell_ball_test.vtk");
    Eigen::Matrix3Xd V, U;
    Eigen::Matrix4Xi TET;
    Eigen::VectorXi fix_b(2);
    Eigen::MatrixX3d bc(2, 3);
    std::vector<double > vecV;
    std::vector<int > vecTet;

    for (size_t i = 0; i < m.num_verts(); i++) {
        cinolib::vec3d p = m.vert(i);
        for (size_t j = 0; j < 3; j++) {
            vecV.push_back(p[j]);
        }
    }

    for (size_t i = 0; i < m.num_polys(); i++) {
        auto ids = m.poly_verts_id(i);
        for (size_t j = 0; j < 4; j++) {
            vecTet.push_back(ids[j]);
        }
    }

    V = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), 3, vecV.size() / 3);
    TET = Eigen::Map<Eigen::Matrix4Xi>(vecTet.data(), 4, vecTet.size() / 4);
    auto n = m.vector_vert_normals();
    fix_b << 728, 1782;
    cinolib::vec3d p1, p2;
    p1 = m.vert(fix_b[0]) + n[fix_b[0]];
    p2 = m.vert(fix_b[1]) + n[fix_b[1]];
    bc.row(0) << p1[0], p1[1], p1[2];
    bc.row(1) << p2[0], p2[1], p2[2];
    std::vector<std::pair<int, Eigen::Vector3d>> boundConstrain;
    for (size_t i = 0; i < fix_b.rows(); i++) {
        boundConstrain.emplace_back(fix_b[i], bc.row(i));
    }

    U = xry_mesh::arap_solve(V, TET, boundConstrain);

    for (size_t i = 0; i < U.cols(); i++ ) {
        auto &p = m.vert(i);
        p[0] = U(0, i);
        p[1] = U(1, i);
        p[2] = U(2, i);
    }
    m.save("arap_test.vtk");
    return 0;
}
