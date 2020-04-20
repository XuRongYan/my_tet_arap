//
// Created by 徐溶延 on 2020/4/20.
//

#ifndef GEN_SHELL_SPACE_MESH_CONVERT_H
#define GEN_SHELL_SPACE_MESH_CONVERT_H

#include <Eigen/Dense>
#include <iostream>

namespace xry_mesh {
    int tri2tet(const Eigen::Matrix3Xd &V,
                const Eigen::Matrix3Xi &F,
                const std::string &fileName);
}


#endif //GEN_SHELL_SPACE_MESH_CONVERT_H
