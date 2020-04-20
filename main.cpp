#include <iostream>
#include <cinolib/meshes/meshes.h>
#include <vector>
#include "arap_tet.h"
#include "cino_io.h"
#include "mesh_convert.h"
#include "vtk.h"

Eigen::Matrix3d rotateMatrix(const Eigen::Vector3d &n, double angle) {
    Eigen::MatrixXd E = Eigen::MatrixXd::Identity(3, 3);
    angle = M_PI * (angle / 180);
    Eigen::MatrixXd A(3, 3);
    A << 0, -n[2], n[1],
            n[2], 0, -n[0],
            -n[1], n[0], 0;
    return E * cos(angle) + (1 - cos(angle)) * n * n.transpose() + sin(angle) * A;
}

int main() {
    Eigen::Matrix3Xd V, U;
    Eigen::Matrix3Xi F;
    Eigen::Matrix4Xi TET;
    Eigen::VectorXi fix_b(2);
    Eigen::MatrixX3d bc(2, 3);
    xry_mesh::readTrimesh("bar1.off", V, F);
    xry_mesh::tri2tet(V, F, "tri2tet.vtk");
    std::vector<std::pair<int, Eigen::Vector3d>> boundConstrain;

    Eigen::Vector3d axis = (V.col(49) - V.col(698)).normalized();
    Eigen::Matrix3d R = rotateMatrix(axis, 90);

    for (size_t i = 49; i <= 97; i++) {
        boundConstrain.emplace_back(i, R * V.col(i));
    }
    xry_mesh::readTetmesh("tri2tet.vtk", V, TET);
    U = xry_mesh::arap_solve(V, TET, boundConstrain);
    std::vector<double> vecV;
    std::vector<int> vecTet;
    for (size_t i = 0; i < U.cols(); i++) {
        vecV.push_back(U(0, i));
        vecV.push_back(U(1, i));
        vecV.push_back(U(2, i));
    }
    for (size_t i = 0; i < TET.cols(); i++) {
        vecTet.push_back(TET(0, i));
        vecTet.push_back(TET(1, i));
        vecTet.push_back(TET(2, i));
        vecTet.push_back(TET(3, i));
    }
    ofstream ofs("result.vtk");
    tet2vtk(ofs, vecV.data(), vecV.size() / 3, vecTet.data(), vecTet.size() / 4);
    return 0;
}
