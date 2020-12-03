#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <LBFGS.h>


using Eigen::VectorXd;
using namespace LBFGSpp;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const VectorXd& x, VectorXd& grad, Eigen::MatrixXd & Hessian)
    {
        double fx = 0.0;
        //std::vector<T> coefficients; // list of non-zeros coefficients
        Hessian.resize(n, n);
        for (int i = 0; i < n; i += 2)
        {
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
            double hi_i = 2 - 400 * (x[i + 1] - 3 * x[i] * x[i]);
            double hi_ii = -400 * x[i];
            double hii_i = -400 * x[i];
            double hii_ii = 200;
           // coefficients.push_back(T(i, i, hi_i));
           // coefficients.push_back(T(i, i + 1, hi_ii));
           // coefficients.push_back(T(i + 1, i, hii_i));
           // coefficients.push_back(T(i + 1, i + 1, hii_ii));

            Hessian(i, i ) = hi_i;
            Hessian(i, i + 1) = hi_ii;
            Hessian(i+1, i + 1) = hii_i;
            Hessian(i+1, i + 1) = hii_ii;

        }
        //Hessian.resize(n, n);
        //Hessian.setFromTriplets(coefficients.begin(), coefficients.end());
        return fx;
    }
};

int main()
{
    const int n = 10;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Rosenbrock fun(n);

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cin.get();
    return 0;
}
