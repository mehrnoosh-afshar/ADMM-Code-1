#pragma once
#include <Eigen/Dense>

// Relevent papers:
// Computing the Singular Value Decomposition of 3x3 matrices with minimal branching and elementary floating point operations, McAdams et al.
// Energetically Consistent Invertible Elasticity, Stomakhin et al.
// Invertible Finite Elements For Robust Simulation of Large Deformation, Irving et al.


namespace fsvd {
	typedef  Eigen::Matrix<double, 3, 1> Vec3;
	typedef  Eigen::Matrix<double, 3, 3> Mat3;
	typedef  Eigen::Matrix<double, 2, 1> Vec2;
	typedef  Eigen::Matrix<double, 2, 2> Mat2;

	// Projection, Singular Values, SVD's U, SVD's V
	static inline void signed_svd_3D(const Mat3& F, Vec3& S, Mat3& U, Mat3& V) {
		using namespace Eigen;

		JacobiSVD< Mat3 > svd(F, ComputeFullU | ComputeFullV);
		S = svd.singularValues();
		U = svd.matrixU();
		V = svd.matrixV();
		Mat3 J = Matrix3d::Identity();
		J(2, 2) = -1.0;

		// Check for inversion: U
		if (U.determinant() < 0.0) {
			U = U * J;
			S[2] = -S[2];
		}

		// Check for inversion: V
		if (V.determinant() < 0.0) {
			Mat3 Vt = V.transpose();
			Vt = J * Vt;
			V = Vt.transpose();
			S[2] = -S[2];
		}

	} // end signed svd

	
	static inline void signed_svd_2D(const Mat2 & F, Vec2 & S, Mat2 & U, Mat2 & V) {
		using namespace Eigen;
		
		JacobiSVD< Mat2 > svd(F, ComputeFullU | ComputeFullV);
		S = svd.singularValues();
		U = svd.matrixU();
		V = svd.matrixV();
		Mat2 J = Matrix2d::Identity();
		J(1, 1) = -1.0;


		// Check for inversion: U
		if (U.determinant() < 0.0) {
			U = U * J;
			S[1] = -S[1];
		}

		// Check for inversion: V
		if (V.determinant() < 0.0) {
			Mat2 Vt = V.transpose();
			Vt = J * Vt;
			V = Vt.transpose();
			S[1] = -S[1];
		}

	} // end signed svd

}
