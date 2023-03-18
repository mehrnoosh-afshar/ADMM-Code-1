#pragma once
#include"MaterialModel.h"
#include <LBFGS.h>
#include <iostream>
#include <Eigen/Core>
#include "LBFGSpp/LineSearchBacktracking.h"


class Triangle_Energy_Term
{
protected:
	typedef Eigen::Matrix<int, 3, 1> Vec3i;
	typedef Eigen::Matrix<double, 2, 1> Vec2;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;

	double scalar_weight;
	int dim;
	
	Vec2 S_ref;
	Eigen::MatrixXd Global_Node_vector;
	std::vector<Vec2> verts;
	virtual void prox(VecX& zi) = 0;
	int global_row_index;
	LBFGSpp::LBFGSParam<double> param;
    LBFGSpp::LBFGSSolver<double> solver;
public:
	Vec3i Element_Node_index;
	double area;
	double get_weight() const { return scalar_weight; }
	void get_local_Dmatrix(std::vector< Eigen::Triplet<double> >& Dmatrix_triplets , Eigen::MatrixXd& D_reduction , std::vector<double>& weights);
	Triangle_Energy_Term( const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop);
	void update(const SparseMat& D, const VecX& x, VecX& z, VecX& u);
	
	// Unless derived from uses linear strain (no area conservation)
	//virtual void prox(VecX& zi);
	//virtual double energy(const VecX& F);
	//virtual double gradient(const VecX& F, VecX& grad);


}; //

template<typename Foo>
	void newton(Foo& f, Eigen::VectorXd& x) {
		int maxIter = 1000;
		int iter = 0;
		int n = x.size();
		Eigen::MatrixXd m_hess;
		Eigen::VectorXd m_grad;
		m_hess.resize(n, n);
		m_grad.resize(n);
		auto epsilon = 1e-7;
		auto tau = 0.8;
		auto alpha = 1.0;
		auto beta = 0.1;
		for (; ;) {
			auto fx = f(x, m_grad, m_hess);

			// ****************** MMODIFYING HESSIAN TO ENSURE POSITIVE DEFINITENESS *************/
			auto eigenValues = m_hess.eigenvalues();
			auto minElem = std::min_element(eigenValues.begin(), eigenValues.end(),
                           [](const std::complex<double>& a, const std::complex<double>& b)
                             { return abs(a) < abs(b); });

			auto M = std::max(0.0, (*minElem).real()) * Eigen::MatrixXd::Identity(n, n);
			auto A = m_hess + M;
			// ****************** MMODIFYING HESSIAN TO ENSURE POSITIVE DEFINITENESS END *************/

			// auto begin2 = std::chrono::high_resolution_clock::now();
			// **************** LINE SEARCH CODE ************/
			auto pk = -m_grad;
			Eigen::MatrixXd temp_1_m_hess;
			Eigen::VectorXd temp_1_m_grad;
			Eigen::MatrixXd temp_2_m_hess;
			Eigen::VectorXd temp_2_m_grad;
			while (f(x+ alpha*(pk), temp_1_m_grad, temp_1_m_hess) > (f(x, temp_2_m_grad, temp_2_m_hess) + alpha*beta*(m_grad.dot(pk)))) {
				alpha = tau * alpha;
			}
			// **************** LINE SEARCH CODE END ************/

			if (iter == maxIter) {
				return;
			}
			auto gnorm = m_grad.norm();
			// TODO: (ddesilva) put these values as constants outside of this function
			if(gnorm <= 1e-3 || gnorm <= 1e-5 * x.norm())
            {
                return;
            }
			// x = x - (A.inverse() * m_grad);
			x = A.llt().solve(A*x - m_grad*alpha);
			iter++;
		}
	}


class HyperElastic_Triangle : public Triangle_Energy_Term 
{
public:

	MaterialModel mechanical_object;
	void prox(VecX& zi);
	//double cost_function( Eigen::VectorXd& S, Eigen::VectorXd& grad);

	HyperElastic_Triangle(const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop);

	// currently broken TODO: (ddesilva): fix 
	template<typename Foo>
	double backtrackingLineserach(Foo& f, Eigen::VectorXd& x) {
		auto alpha = 1.0;
		auto rho = 0.8;
		auto c = 1e-4;

		Eigen::MatrixXd m_hess;
		Eigen::VectorXd m_grad;
		Eigen::MatrixXd m_hess2;
		Eigen::VectorXd m_grad2;

		while (f(x + Eigen::VectorXd::Ones(x.size()) * alpha * (-m_grad.dot(x)), m_grad, m_hess) > f(x, m_grad2, m_hess2) + c * alpha * m_grad.dot(x) * (-m_grad.dot(x))) {
			alpha *= rho;
		}
		return alpha;
	}

	// takes 0.23s 
	template<typename Foo>
	void newton(Foo& f, Eigen::VectorXd& x) {
		int maxIter = 100;
		int iter = 0;
		int n = x.size();
		Eigen::MatrixXd m_hess;
		Eigen::VectorXd m_grad;
		m_hess.resize(n, n);
		m_grad.resize(n);
		for (; ;) {
			auto fx = f(x, m_grad, m_hess);

			if (iter == maxIter) {
				return;
			}
			auto gnorm = m_grad.norm();
			// TODO: (ddesilva) put these values as constants outside of this function
			if(gnorm <= 1e-3 || gnorm <= 1e-5 * x.norm())
            {
                return;
            }
			// x = x - (m_hess.inverse() * m_grad);
			x = m_hess.householderQr().solve(m_hess*x - m_grad);
			iter++;
		}
	}

	// takes around 1.9s 
	template<typename Foo>
	void newtonWithLineSearch(Foo& f, Eigen::VectorXd& x) {
		int maxIter = 100; // TODO (ddesilva): intialize these constants outside 
		int iter = 0;
		Eigen::MatrixXd m_hess;
		Eigen::VectorXd m_grad;

		Eigen::VectorXd m_drt;
		Eigen::VectorXd m_xp;
		Eigen::VectorXd m_fx;

		m_xp.noalias() = x;

		LBFGSpp::LBFGSParam<double> m_param; // TODO (ddesilva): intialize these constants outside 
		m_param.epsilon = 1e-3;
		m_param.max_iterations = 200;
		m_param.linesearch = 2;

		double step = 1.0;

		for (; ;) {
			auto fx = f(x, m_grad, m_hess);
			m_drt = - (m_hess.inverse() * m_grad);
			if (iter == 0) {
				step = 1 / m_drt.norm();
			}

			auto gnorm = m_grad.norm();
			if(gnorm <= 1e-3 || gnorm <= 1e-5 * x.norm())
            {
                return;
            }
			if (iter == maxIter) {
				return;
			}

			LBFGSpp::LineSearchBacktracking<double>::LineSearch(f, fx, x, m_grad, m_hess, step, m_drt, m_xp, m_param);
			step = 1.0;
			
			iter++;
		}
	}


class cost_function  {
	typedef Eigen::Matrix<int, 3, 1> Vec3i;
	typedef Eigen::Matrix<double, 2, 1> Vec2; 
public:
	HyperElastic_Triangle HTRI;
	Vec2 S_ref;
	cost_function (const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop,Vec2& S_ref0);
	double operator ()(const Eigen::VectorXd& S, Eigen::VectorXd& grad , Eigen::MatrixXd& Hessian);

};