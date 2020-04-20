//
// Created by 徐溶延 on 2020/4/19.
//

#include "arap_tet.h"
#include <iostream>

namespace xry_mesh {
    const int w = 1e8;

    int computeIdeal(const Eigen::Matrix3Xd &V,
                     const Eigen::Matrix4Xi &TET,
                     std::vector<Eigen::Matrix<double, 3, 4>> &idealElem) {
        Eigen::Matrix<double, 3, 4> ideal;
        for (size_t i = 0; i < TET.cols(); i++) {
            for (size_t j = 0; j < 4; j++) {
                ideal.col(j) = V.col(TET(j, i));
            }
            idealElem.emplace_back(ideal);
        }
        return 0;
    }

    double computeError(const std::vector<Eigen::Matrix3d> &R,
                        const std::vector<Eigen::Matrix3d> &deformGrad) {
        double error = 0;
        for (size_t i = 0; i < R.size(); i++) {
            error += (deformGrad[i] - R[i]).squaredNorm();
        }
        return error;
    }

    int optimizeRotation(const Eigen::Matrix3d &J, Eigen::Matrix3d &R) {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        R = svd.matrixU() * svd.matrixV().transpose();
        if (fabs(svd.singularValues()(2)) < 1e-8) {
            R.setIdentity();
            return 1;
        }
        if (R.determinant() < 0) {
            Eigen::MatrixXd svdV = svd.matrixV();
            svdV.col(2) = -svdV.col(2);
            R = svd.matrixU() * svdV.transpose();
        }
        return 0;
    }

    int localPhase(const Eigen::Matrix4Xi &TET,
                   const std::vector<Eigen::Matrix3d> &deformGrad,
                   std::vector<Eigen::Matrix3d> &R) {
        Eigen::Matrix3d m;
        for (int i = 0; i < TET.cols(); i++) {
            optimizeRotation(deformGrad[i], R[i]);
        }
        return 0;
    }

    int globalPhase(const Eigen::Matrix4Xi &TET,
                    const std::vector<Eigen::Matrix3d> &R,
                    const Eigen::SparseMatrix<double> &G,
                    const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &llt,
                    Eigen::MatrixX3d &x) {
        Eigen::Matrix3Xd bt(3, 3 * TET.cols());
        for (size_t i = 0; i < TET.cols(); i++) {
            for (size_t j = 0; j < 3; j++) {
                bt.col(3 * i + j) = R[i].col(j);
            }
        }

        Eigen::MatrixX3d Gtb = G.transpose() * bt.transpose();
        for (const auto & pair : bc) {
            Gtb.row(pair.first) += w * pair.second;
        }
        x = llt.solve(Gtb);
        return 0;
    }

    int computeGradients(const Eigen::Matrix3Xd &V,
                         const Eigen::Matrix4Xi &TET,
                         const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem,
                         Eigen::SparseMatrix<double> &G,
                         std::vector<Eigen::Matrix3d> &deformGrad) {
        Eigen::MatrixXd X(3, 4);
        std::vector<Eigen::Matrix<double, 3, 4>> gradPhis;
        std::vector<Eigen::Triplet<double>> triplets;

        X << 1, 0, 0, -1,
             0, 1, 0, -1,
             0, 0, 1, -1;

        for (size_t i = 0; i < TET.cols(); i++) {
            Eigen::Matrix3d delta;
            delta.col(0) = idealElem[i].col(0) - idealElem[i].col(3);
            delta.col(1) = idealElem[i].col(1) - idealElem[i].col(3);
            delta.col(2) = idealElem[i].col(2) - idealElem[i].col(3);
            gradPhis.emplace_back(delta.transpose().inverse() * X);
        }

        for (size_t i = 0; i < TET.cols(); i++) {
            for (size_t j = 0; j < 4; j++) {
                for (size_t k = 0; k < 3; k++) {
                    triplets.emplace_back(3 * i + k, TET(j, i), gradPhis[i](k, j));
                }
            }
        }

        for (size_t i = 0; i < TET.cols(); i++) {
            deformGrad.emplace_back(gradPhis[i] * idealElem[i].transpose());
        }

        G.resize(3 * TET.cols(), V.cols());
        G.setFromTriplets(triplets.begin(), triplets.end());
        return 0;
    }

    int precompute(const Eigen::Matrix3Xd &V,
                   const Eigen::Matrix4Xi &TET,
                   const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                   const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem,
                   Eigen::SparseMatrix<double> &G,
                   std::vector<Eigen::Matrix3d> &deformGrad,
                   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &llt) {
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::SparseMatrix<double> L;
        computeGradients(V, TET, idealElem, G, deformGrad);
        L = G.transpose() * G;
        for (const auto &pair : bc) {
            L.coeffRef(pair.first, pair.first) += w;
        }
        llt.compute(L);
        return 0;
    }

    Eigen::Matrix3Xd arap_solve(const Eigen::Matrix3Xd &V,
                                const Eigen::Matrix4Xi &TET,
                                const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                                const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem) {
        const int max_iter = 100;
        Eigen::SparseMatrix<double> G;  //梯度算子
        std::vector<Eigen::Matrix3d> deformGrad;
        std::vector<Eigen::Matrix3d> R(TET.cols());
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
        Eigen::MatrixX3d x = V.transpose();
        precompute(V, TET, bc, idealElem, G, deformGrad, llt);
        double err0 = -1, err1 = 0;
        size_t iter = 0;
        while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
            err0 = err1;
            localPhase(TET, deformGrad, R);
            globalPhase(TET, R, G, bc, llt, x);
            err1 = computeError(R, deformGrad);
            std::cout << ++iter << ":" << err1 << std::endl;
        }
        return x.transpose();
    }


    Eigen::Matrix3Xd arap_solve(const Eigen::Matrix3Xd &V,
                                const Eigen::Matrix4Xi &TET,
                                const std::vector<std::pair<int, Eigen::Vector3d>> &bc) {
        std::vector<Eigen::Matrix<double, 3, 4>> idealElem;
        computeIdeal(V, TET, idealElem);
        return arap_solve(V, TET, bc, idealElem);
    }
}