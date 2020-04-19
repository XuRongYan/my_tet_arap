#include <iostream>
#include <cinolib/meshes/meshes.h>
#include <vector>
#include "arap.h"

int main() {
    cinolib::Tetmesh<> m("shell_ball_test.vtk");
    Eigen::Matrix3Xd V;
    Eigen::Matrix4Xi TET;
    Eigen::VectorXi fix_b(2);
    Eigen::MatrixX3d bc(2, 3), U;
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
    fix_b << 0, (V.cols() / 6);
    cinolib::vec3d p1, p2;
    p1 = m.vert(fix_b(0, 0)) + 3 * n[fix_b(0, 0)];
    p2 = m.vert(fix_b(1, 0));
    bc.row(0) << p1[0], p1[1], p1[2];
    bc.row(1) << p2[0], p2[1], p2[2];
    xry_mesh::precomputation(V.transpose(), TET.transpose(), fix_b, bc);
    xry_mesh::solve(U);
    for (size_t i = 0; i < U.rows(); i++) {
        auto &p = m.vert(i);
        p[0] = U(i, 0);
        p[1] = U(i, 1);
        p[2] = U(i, 2);
    }
    m.save("arap_test.vtk");
    return 0;
}
