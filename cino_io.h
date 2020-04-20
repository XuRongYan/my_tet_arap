//
// Created by 徐溶延 on 2020/4/20.
//

#ifndef MY_TET_ARAP_CINO_IO_H
#define MY_TET_ARAP_CINO_IO_H

#include <iostream>
#include <Eigen/Dense>

namespace xry_mesh {
    int readTrimesh(const std::string &fileName,
                    Eigen::Matrix3Xd &V,
                    Eigen::Matrix3Xi &F);

    int readTetmesh(const std::string &fileName,
                    Eigen::Matrix3Xd &V,
                    Eigen::Matrix4Xi &TET);
}


#endif //MY_TET_ARAP_CINO_IO_H
