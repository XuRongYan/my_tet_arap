//
// Created by 徐溶延 on 2020/4/20.
//

#include "cino_io.h"
#include <cinolib/meshes/meshes.h>

namespace xry_mesh {
    int readTrimesh(const std::string &fileName,
                    Eigen::Matrix3Xd &V,
                    Eigen::Matrix3Xi &F) {
        std::vector<double > vecV;
        std::vector<int > vecF;
        cinolib::Trimesh<> m(fileName.c_str());
        for (size_t i = 0; i < m.num_verts(); i++) {
            cinolib::vec3d p = m.vert(i);
            for (size_t j = 0; j < 3; j++) {
                vecV.push_back(p[j]);
            }
        }

        for (size_t i = 0; i < m.num_polys(); i++) {
            auto ids = m.poly_verts_id(i);
            for (size_t j = 0; j < 3; j++) {
                vecF.push_back(ids[j]);
            }
        }
        V = Eigen::Map<Eigen::Matrix3Xd>(vecV.data(), 3, vecV.size() / 3);
        F = Eigen::Map<Eigen::Matrix3Xi>(vecF.data(), 3, vecF.size() / 3);
        return 0;
    }

    int readTetmesh(const std::string &fileName,
                    Eigen::Matrix3Xd &V,
                    Eigen::Matrix4Xi &TET) {
        std::vector<double > vecV;
        std::vector<int > vecTet;
        cinolib::Tetmesh<> m(fileName.c_str());
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
        return 0;
    }
}