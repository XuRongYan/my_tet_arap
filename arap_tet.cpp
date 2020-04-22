//
// Created by 徐溶延 on 2020/4/19.
//

#include "arap_tet.h"
#include <iostream>

namespace xry_mesh {
    const int w = 1e8;

    int updateDeformGrad(const Eigen::Matrix4Xi &TET,
                         const Eigen::SparseMatrix<double> &G,
                         const Eigen::MatrixX3d &x,
                         std::vector<Eigen::Matrix3d> &deformGrad) {
        Eigen::MatrixXd Df = G * x;
        assert(G.rows() == TET.cols() * 3 && G.cols() == x.rows());
        deformGrad.clear();
        for (int i = 0; i < TET.cols(); i++) {
            Eigen::Matrix3d df = Df.block(i * 3, 0, 3, 3);
            deformGrad.emplace_back(df.transpose());
        }
        return 0;
    }

    double computeVol(const Eigen::Matrix<double, 3, 4> &ideal) {
        Eigen::Vector3d a, b, c;
        a = ideal.col(1) - ideal.col(0);
        b = ideal.col(2) - ideal.col(0);
        c = ideal.col(3) - ideal.col(0);
        double res = (a.cross(b)).dot(c) / 6;
        assert(res > 0);
        return res;
    }

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

    double computeError(const Eigen::MatrixX3d &x,
                        const Eigen::Matrix4Xi &TET,
                        const std::vector<Eigen::Matrix3d> &R,
                        const std::vector<Eigen::Matrix3d> &invDeltas,
                        const Eigen::VectorXd &vols) {
        double error = 0;
        for (size_t i = 0; i < R.size(); i++) {
            Eigen::Matrix3d m;
            m.col(0) = x.row(TET(1, i)) - x.row(TET(0, i));
            m.col(1) = x.row(TET(2, i)) - x.row(TET(0, i));
            m.col(2) = x.row(TET(3, i)) - x.row(TET(0, i));
            error += vols[i] * (m * invDeltas[i] - R[i]).squaredNorm();
        }
        return error;
    }

    int optimizeRotation(const Eigen::Matrix3d &J, Eigen::Matrix3d &R) {
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
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
                   const Eigen::MatrixX3d &x,
                   const std::vector<Eigen::Matrix3d> &invDeltas,
                   std::vector<Eigen::Matrix3d> &R) {
        for (int i = 0; i < TET.cols(); i++) {
            Eigen::Matrix3d m;
            m.col(0) = x.row(TET(1, i)) - x.row(TET(0, i));
            m.col(1) = x.row(TET(2, i)) - x.row(TET(0, i));
            m.col(2) = x.row(TET(3, i)) - x.row(TET(0, i));
            optimizeRotation(m * invDeltas[i], R[i]);
        }
        return 0;
    }

    int globalPhase(const Eigen::Matrix4Xi &TET,
                    const std::vector<Eigen::Matrix3d> &R,
                    const Eigen::SparseMatrix<double> &Gt,
                    const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                    const Eigen::VectorXd &vols,
                    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &llt,
                    Eigen::MatrixX3d &x) {

        Eigen::Matrix3Xd bt(3, 3 * TET.cols());
        for (size_t i = 0; i < TET.cols(); i++) {
            for (size_t j = 0; j < 3; j++) {
                bt.col(3 * i + j) = sqrt(vols[i]) * R[i].col(j);
            }
        }

        Eigen::MatrixX3d Gtb = Gt * bt.transpose();

        for (const auto &pair : bc) {
            Gtb.row(pair.first) += w * pair.second;
        }
        x = llt.solve(Gtb);
        return 0;
    }

    int computeGradients(const Eigen::Matrix3Xd &V,
                         const Eigen::Matrix4Xi &TET,
                         const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem,
                         const Eigen::VectorXd &vols,
                         std::vector<Eigen::Matrix3d> &invDeltas,
                         Eigen::SparseMatrix<double> &Gt) {
        std::vector<Eigen::Matrix<double, 3, 4>> gradPhis;
        std::vector<Eigen::Triplet<double>> triplets;


        for (size_t i = 0; i < TET.cols(); i++) {
            Eigen::Matrix3d delta;
            delta.col(0) = idealElem[i].col(1) - idealElem[i].col(0);
            delta.col(1) = idealElem[i].col(2) - idealElem[i].col(0);
            delta.col(2) = idealElem[i].col(3) - idealElem[i].col(0);
            invDeltas.emplace_back(delta.inverse());
        }

        for (size_t i = 0, col = 0; i < TET.cols(); i++) {
            for (size_t j = 0; j < 3; j++, col++) {
                double sum = 0;
                for (size_t k = 0; k < 3; k++) {
                    const double val = invDeltas[i](k, j) * sqrt(vols[i]);
                    sum -= val;
                    triplets.emplace_back(TET(k + 1, i), col, val);
                }
                triplets.emplace_back(TET(0, i), col, sum);
            }
        }

        Gt.resize(V.cols(), 3 * TET.cols());
        Gt.setFromTriplets(triplets.begin(), triplets.end());
        return 0;
    }

    int precompute(const Eigen::Matrix3Xd &V,
                   const Eigen::Matrix4Xi &TET,
                   const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                   const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem,
                   Eigen::VectorXd &vols,
                   Eigen::SparseMatrix<double> &Gt,
                   std::vector<Eigen::Matrix3d> &invDeltas,
                   Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &llt) {
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::SparseMatrix<double> L;

        //compute vols
        vols.resize(TET.cols());
        for (size_t i = 0; i < idealElem.size(); i++) {
            vols[i] = computeVol(idealElem[i]);
        }

        //compute gradient operator(size = 3m * n)
        computeGradients(V, TET, idealElem, vols, invDeltas, Gt);

        //L = GtG
        L = Gt * Gt.transpose();
        //add bound constrains
        for (const auto &pair : bc) {
            L.coeffRef(pair.first, pair.first) += w;
        }
        //cholesky decomposition
        llt.compute(L);
        return 0;
    }

    Eigen::Matrix3Xd arap_solve(const Eigen::Matrix3Xd &V,
                                const Eigen::Matrix4Xi &TET,
                                const std::vector<std::pair<int, Eigen::Vector3d>> &bc,
                                const std::vector<Eigen::Matrix<double, 3, 4>> &idealElem) {
        const int max_iter = 1000;
        assert(idealElem.size() == TET.cols());
        Eigen::SparseMatrix<double> Gt;  //梯度算子
        std::vector<Eigen::Matrix3d> invDeltas;
        std::vector<Eigen::Matrix3d> R(TET.cols());
        Eigen::VectorXd vols;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
        Eigen::MatrixX3d x = V.transpose();
        precompute(V, TET, bc, idealElem, vols, Gt, invDeltas, llt);

        double err0 = -1, err1 = 0;
        size_t iter = 0;

        while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
            err0 = err1;
            localPhase(TET, x, invDeltas, R);
            globalPhase(TET, R, Gt, bc, vols, llt, x);
            err1 = computeError(x, TET, R, invDeltas, vols);
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