#pragma once
#ifndef ADMM_LINEARSOLVER_H
#define ADMM_LINEARSOLVER_H


#include <Eigen/SparseCholesky>
#include <unordered_map>
#include <memory>
#ifdef EIGEN_USE_MKL_VML
#include <Eigen/PardisoSupport>
#endif

namespace ADMM {

	//
	//	Linear solver base class
	// Base class
	class LinearSolver {
	public:
		typedef Eigen::Matrix<double, 3, 1> Vec3;
		typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
		typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;


		// Updates solver data
		virtual void update_system(const SparseMat& A_) = 0;

		// Solve Ax=b subject to constraints
		virtual int solve(VecX& x, const VecX& b) = 0;

		// Helper class for avoiding divide-by-zero
		static bool is_zero(double x) { return std::abs(x) < std::numeric_limits<double>::min(); }

	}; // end class linear solver


	//
	//	Classic LDLT solver that does not handle constraints (SCA 2016)
	//
	class LDLTSolver : public LinearSolver {
	public:
		typedef Eigen::Matrix<double, 3, 1> Vec3;
		typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
		typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;

#ifdef EIGEN_USE_MKL_VML
		typedef Eigen::PardisoLDLT< Eigen::SparseMatrix<double> > Cholesky;
#else
		typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > Cholesky;
#endif

		SparseMat A;
		std::unique_ptr<Cholesky> m_cholesky;

		LDLTSolver() {
			m_cholesky = std::unique_ptr<Cholesky>(new Cholesky());
		}

		// Does a cholesky factorization on the system matrix
		void update_system(const SparseMat& A_) {
			int dim = A_.rows();
			if (dim != A_.cols() || dim == 0) { throw std::runtime_error("**LDLTSolver Error: Bad dimensions in A"); }
			A = A_;
			m_cholesky->compute(A);
		}

		// Solves for x given A, linear constraints (C), and pinning subspace (P)
		int solve(VecX& x, const VecX& b0) {
			x = m_cholesky->solve(b0);

			return 1;
		}

	}; // end class linear solver

} // ns admm

#endif