//
// Created by 徐溶延 on 2020/4/20.
//

#include "mesh_convert.h"
#include <vector>
#include <fstream>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include "vtk.h"

namespace xry_mesh {
    int tri2tet(const Eigen::Matrix3Xd &V,
                const Eigen::Matrix3Xi &F,
                const std::string &fileName) {
        Eigen::MatrixX3d TV, Vt;
        Eigen::MatrixX4i TT;
        Eigen::MatrixX3i TF, Ft;
        Vt = V.transpose();
        Ft = F.transpose();
        igl::copyleft::tetgen::tetrahedralize(Vt, Ft, "pq1.414", TV, TT, TF);
        std::vector<double > vecV;
        std::vector<int> vecTet;
        for (size_t i = 0; i < TV.rows(); i++) {
            vecV.push_back(TV(i, 0));
            vecV.push_back(TV(i, 1));
            vecV.push_back(TV(i, 2));
        }

        for (size_t i = 0; i < TT.rows(); i++) {
            vecTet.push_back(TT(i, 0));
            vecTet.push_back(TT(i, 1));
            vecTet.push_back(TT(i, 2));
            vecTet.push_back(TT(i, 3));
        }
        std::ofstream ofs(fileName);
        tet2vtk(ofs, vecV.data(), vecV.size() / 3, vecTet.data(), vecTet.size() / 4);
        return 0;
    }
}