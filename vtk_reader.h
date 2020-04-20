#ifndef VTK_READER_H_H
#define VTK_READER_H_H

#include <string>
#include <Eigen/Eigen>

namespace io {
int read_mesh_from_vtk(const std::string &file_name,
                       Eigen::Matrix3Xd &V,
                       Eigen::MatrixXi &T,
                       int mesh_type = 1);


int read_mesh_from_vtk(const std::string &file_name,
                       Eigen::Matrix3Xd &V,
                       Eigen::Matrix4Xi &T,
                       Eigen::MatrixXi &P,
                       Eigen::MatrixXi &H);

}


#endif
