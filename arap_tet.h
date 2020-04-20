//
// Created by 徐溶延 on 2020/4/19.
//

#ifndef MY_TET_ARAP_ARAP_TET_H
#define MY_TET_ARAP_ARAP_TET_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace xry_mesh {
    Eigen::Matrix3Xd arap_solve(const Eigen::Matrix3Xd &V,
                                const Eigen::Matrix4Xi &TET,
                                const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                                const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem);

    Eigen::Matrix3Xd arap_solve(const Eigen::Matrix3Xd &V,
                                const Eigen::Matrix4Xi &TET,
                                const std::vector<std::pair<int, Eigen::Vector3d>> &bc);
}


#endif //MY_TET_ARAP_ARAP_TET_H
