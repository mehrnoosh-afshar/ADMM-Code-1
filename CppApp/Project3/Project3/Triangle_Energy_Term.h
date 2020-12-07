#pragma once
#include"MaterialModel.h"
#include <LBFGS.h>
#include <iostream>
#include <Eigen/Core>




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
		auto epsilon = 1e-7 * Eigen::MatrixXd::Identity(n, n);
		for (; ;) {
			auto fx = f(x, m_grad, m_hess);
			// auto A = m_hess + epsilon;
			// std::cout << "Iter : " << iter << " start" << std::endl;
			// std::cout << "grad norm" << std::endl;
			// std::cout << m_grad.norm() << std::endl;
			// std::cout << "eigen values of hessian" << std::endl;
			// std::cout << m_hess.eigenvalues() << std::endl;
			// std::cout << "value of cost function" << std::endl;
			// std::cout << fx << std::endl << std::endl;

			if (iter == maxIter) {
				return;
			}
			auto gnorm = m_grad.norm();
			// TODO: (ddesilva) put these values as constants outside of this function
			if(gnorm <= 1e-3 || gnorm <= 1e-5 * x.norm())
            {
                return;
            }
			// x = x - 0.5*(m_hess.inverse() * m_grad);
			x = m_hess.llt().solve(m_hess*x - m_grad);
			iter++;
		}
	}


class HyperElastic_Triangle : public Triangle_Energy_Term 
{
public:

	MaterialModel mechanical_object;
	void prox(VecX& zi);
	//double cost_function( Eigen::VectorXd& S, Eigen::VectorXd& grad);
	HyperElastic_Triangle( const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop);
};

class cost_function  {
	typedef Eigen::Matrix<int, 3, 1> Vec3i;
	typedef Eigen::Matrix<double, 2, 1> Vec2; 
public:
	HyperElastic_Triangle HTRI;
	Vec2 S_ref;
	cost_function (const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop,Vec2& S_ref0);
	double operator ()(const Eigen::VectorXd& S, Eigen::VectorXd& grad , Eigen::MatrixXd& Hessian);

};