//
// Created by 徐溶延 on 2020/4/14.
//

#ifndef MY_TET_ARAP_ARAP_H
#define MY_TET_ARAP_ARAP_H

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cinolib/meshes/meshes.h>


namespace xry_mesh {
    int precomputation(Eigen::MatrixX3d _V,
                       Eigen::MatrixX4i _TET,
                       Eigen::VectorXi _fix_b,
                       Eigen::MatrixX3d _bc,
                       std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > _idealElem);

    int precomputation(Eigen::MatrixX3d _V,
                       Eigen::MatrixX4i _TET,
                       Eigen::VectorXi _fix_b,
                       Eigen::MatrixX3d _bc);

    int solve(Eigen::MatrixX3d &U);
}


#endif //MY_TET_ARAP_ARAP_H
