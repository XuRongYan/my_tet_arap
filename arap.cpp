//
// Created by 徐溶延 on 2020/4/14.
//

#include "arap.h"
#include <Eigen/Sparse>
#include <utility>

namespace xry_mesh {
    const double w = 1e8;
    int max_iter = 100;//max iteration
    Eigen::VectorXi fix_b;//fixed id
    Eigen::MatrixX3d V;
    Eigen::MatrixX4i TET;
    Eigen::MatrixX3d bc;//fixed position
    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > idealElem;//ideal tets
    Eigen::SparseMatrix<double> Gt, L;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
    std::vector<Eigen::Matrix3d> inv, vecR, vecDf;
    std::vector<Eigen::MatrixXd> vecVolGrad;
    std::vector<Eigen::Triplet<double>> triplets;

    double computeIdeal() {
        Eigen::MatrixX3d tmpV = V;
        Eigen::Matrix<double, 3, 4> ideal;
        for (size_t i = 0; i < TET.rows(); i++) {
            for (size_t j = 0; j < 4; j++) {
                ideal.col(j) = V.row(TET(i, j));
            }
            idealElem.emplace_back(ideal);
        }
        return 0;
    }

    double computeError() {
        double error = 0;
        for (size_t i = 0; i < TET.rows(); i++) {
            error += (vecDf[i] - vecR[i]).squaredNorm();
        }
        return error;
    }

    int computeGradients() {
        Eigen::MatrixXd X(3, 4), gradPhi;
        X << 1, 0, 0, -1,
                0, 1, 0, -1,
                0, 0, 1, -1;

        for (int i = 0; i < TET.rows(); i++) {
            Eigen::MatrixXd delta(3, 3);
            delta.row(0) = (idealElem[i].col(0) - idealElem[i].col(3)).transpose();
            delta.row(1) = (idealElem[i].col(1) - idealElem[i].col(3)).transpose();
            delta.row(2) = (idealElem[i].col(2) - idealElem[i].col(3)).transpose();
            Eigen::MatrixXd df = delta.inverse() * X;
            vecVolGrad.emplace_back(df);
        }

        for (size_t i = 0; i < TET.rows(); i++) {
            for (size_t j = 0; j < 4; j++) {
                for (size_t k = 0; k < 3; k++) {
                    triplets.emplace_back(3 * i + k, TET(i, j), vecVolGrad[i](k, j));
                }
            }
        }

        for (size_t i = 0; i < TET.rows(); i++) {
            vecDf.emplace_back(vecVolGrad[i] * idealElem[i].transpose());
        }

        Eigen::SparseMatrix<double> G(3 * TET.rows(), V.rows());
        G.setFromTriplets(triplets.begin(), triplets.end());
        Gt = G.transpose();
        return 0;
    }

    int optimizeR(const Eigen::MatrixXd &J, Eigen::MatrixXd &R) {
        const int dim = J.rows();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
        R = svd.matrixU() * svd.matrixV().transpose();
        if (fabs(svd.singularValues()(dim - 1)) < 1e-8) {
            R.setIdentity();
            return 1;
        }
        if (R.determinant() < 0) {
            Eigen::MatrixXd svdV = svd.matrixV();
            svdV.col(dim - 1) = -svdV.col(dim - 1);
            R = svd.matrixU() * svdV.transpose();
        }
        return 0;
    }

    //local
    int optimizeRs() {
        Eigen::Matrix3d m;
        Eigen::MatrixXd J(3, 3);
        for (int i = 0; i < TET.rows(); i++) {
            optimizeR(vecDf[i].inverse(), J);
            vecR[i] = J;
        }
        return 0;
    }

    //global
    int optimizeGradients() {
        Eigen::MatrixX3d b(3 * TET.rows(), 3);
        for (size_t i = 0; i < TET.rows(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                b.row(3 * i + j) = vecR[i].col(j);
            }
        }
        Eigen::MatrixX3d Gtb = Gt * b;
        for (size_t i = 0; i < bc.rows(); ++i) {
            Gtb.row(fix_b(i, 0)) += w * bc.row(i);
        }
        Eigen::MatrixX3d x = llt.solve(Gtb);
        cout << (V - x).squaredNorm() << endl;
        V = x;
        return 0;
    }

    int solve(Eigen::MatrixX3d &U) {
        double err0 = -1, err1 = 0;
        size_t iter = 0;
        while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
            err0 = err1;
            optimizeRs();
            optimizeGradients();
            err1 = computeError();
            std::cout << ++iter << ":" << err1 << std::endl;
        }
        U = V;
        return 0;
    }

    int precomputation(Eigen::MatrixX3d _V,
                       Eigen::MatrixX4i _TET,
                       Eigen::VectorXi _fix_b,
                       Eigen::MatrixX3d _bc,
                       std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > _idealElem) {
        xry_mesh::V = std::move(_V);
        xry_mesh::TET = std::move(_TET);
        xry_mesh::fix_b = std::move(_fix_b);
        xry_mesh::bc = std::move(_bc);
        xry_mesh::idealElem = std::move(_idealElem);
        vecR.resize(TET.rows());
        computeGradients();
        L = Gt * Gt.transpose();
        for (size_t i = 0; i < xry_mesh::bc.rows(); ++i) {
                  L.coeffRef(xry_mesh::fix_b(i, 0), xry_mesh::fix_b(i, 0)) += w;
        }
        llt.compute(L);
        return 0;
    }

    int precomputation(Eigen::MatrixX3d _V,
                       Eigen::MatrixX4i _TET,
                       Eigen::VectorXi _fix_b,
                       Eigen::MatrixX3d _bc) {
        xry_mesh::V = std::move(_V);
        xry_mesh::TET = std::move(_TET);
        xry_mesh::fix_b = std::move(_fix_b);
        xry_mesh::bc = std::move(_bc);
        computeIdeal();
        vecR.resize(TET.rows());
        computeGradients();
        L = Gt * Gt.transpose();
        for (size_t i = 0; i < xry_mesh::bc.rows(); ++i) {
            L.coeffRef(xry_mesh::fix_b(i, 0), xry_mesh::fix_b(i, 0)) += w;
        }
        llt.compute(L);
        return 0;
    }
}







