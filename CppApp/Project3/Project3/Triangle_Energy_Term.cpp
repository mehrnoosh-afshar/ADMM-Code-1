#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include"Triangle_Energy_Term.h"
#include"signedSVD.h"

Triangle_Energy_Term::Triangle_Energy_Term(const Vec3i& Element_Node_index, Eigen::MatrixXd & Global_Node_vector, MaterialModel& mechanical_prop):solver(param)
{
	
	Eigen::Matrix3d D;
	D.row(0) << Global_Node_vector.row(Element_Node_index[0] - 1).col(0), Global_Node_vector.row(Element_Node_index[0] - 1).col(1), 1;
	D.row(1) << Global_Node_vector.row(Element_Node_index[1] - 1).col(0), Global_Node_vector.row(Element_Node_index[1] - 1).col(1), 1;
	D.row(2) << Global_Node_vector.row(Element_Node_index[2] - 1).col(0), Global_Node_vector.row(Element_Node_index[2] - 1).col(1), 1;
	area = 0.5 * D.determinant();

	if (area < 0) {
		throw std::runtime_error("**TriEnergyTerm Error: Inverted initial pose Negative area look at Triangle_Energy_Term.cpp ");
	}
	double k = mechanical_prop.bulk_modulus;
	scalar_weight = std::sqrt(k * area);

	
	for (int i = 0; i < Element_Node_index.size(); i++)
	{ 
		verts.push_back(Global_Node_vector.row(Element_Node_index[i] - 1).transpose());
	}

	dim = 2; // for trianle element 

	param.epsilon = 1e-3;
	param.max_iterations = 200;
	param.linesearch = 2;

}

void Triangle_Energy_Term::get_local_Dmatrix(std::vector< Eigen::Triplet<double> >& Dmatrix_triplets, Eigen::MatrixXd& D_reduction, std::vector<double>& weights) {

	
	global_row_index = weights.size();

	// update weight vector 
	//std::cout << weights.size() << std::endl;
	for (int i = 0; i < (dim*dim); ++i) { weights.push_back(scalar_weight); }
	//std::cout << weights.size() << std::endl;
	//std::cout << dim * dim << std::endl;
	// update D matrix 
	Eigen::Matrix<double, 2, 2> edges; // B matrix F=b*inv(B);
	edges.col(0) = verts[1] - verts[0];
	edges.col(1) = verts[2] - verts[0];

	Eigen::Matrix<double, 2, 2> edges_inv = edges.inverse();

	Eigen::Matrix<double, 3, 2> S; // Selector
	S.setZero();
	S(0, 0) = -1;	S(0, 1) = -1;	
	S(1, 0) = 1;
	S(2, 1) = 1;
	
	Eigen::Matrix<double, 3, 2> D0 = S * edges_inv;
	D_reduction.resize(2, 3);
	D_reduction = D0.transpose(); // Reduction multiplation of this to node corrdinate gives the transpose of F matrix

	const int rows[2] = { 0, dim};
	const int cols[3] = { dim * (Element_Node_index[0] - 1), dim * (Element_Node_index[1] - 1),  dim * (Element_Node_index[2] - 1) };
	for (int r = 0; r < D_reduction.rows(); ++r) {
		for (int c = 0; c < D_reduction.cols(); ++c) {
    		double value = D_reduction(r, c);
			for (int j = 0; j < dim; ++j) {
				Dmatrix_triplets.emplace_back(rows[r] + j, cols[c] + j, value);
			}
		}
	}
}

 void Triangle_Energy_Term::update(const SparseMat& D, const VecX& x, VecX& z, VecX& u) {
	int dof = x.rows();
	int dim_problem = dim*dim;
	VecX Dix = D.block(global_row_index, 0, dim_problem, dof) * x;
	VecX ui = u.segment(global_row_index, dim_problem);
	VecX zi = Dix + ui;
	prox(zi);
	ui += (Dix - zi);
	u.segment(global_row_index, dim_problem) = ui;
	z.segment(global_row_index, dim_problem) = zi;

}

HyperElastic_Triangle::HyperElastic_Triangle(const Vec3i& Element_Node_index0, Eigen::MatrixXd& Global_Node_vector0, MaterialModel& mechanical_prop)
	 :Triangle_Energy_Term(Element_Node_index0, Global_Node_vector0, mechanical_prop)
{
	mechanical_object = mechanical_prop;
	Element_Node_index = Element_Node_index0;
	Global_Node_vector = Global_Node_vector0;
}

void HyperElastic_Triangle::prox(VecX& zi) {

	typedef Eigen::Matrix<double, 4, 1> Vec4;
	typedef Eigen::Matrix<double, 2, 2> Mat2;

	// at the begining zi=dix+ui;
	Mat2 F = Eigen::Map<Mat2>(zi.data());
	Vec2 S; Mat2 U, V;

	fsvd::signed_svd_2D(F, S, U, V);

	//std::cout << S << std::endl;
	// S_per is the refrence S 
	S_ref = S; 

	//std::cout << S_ref << std::endl;
	// If everything is very low, It is collapsed to a point and the minimize
	// will likely fail. So we'll just inflate it a bit.
	const double eps = 1e-6;
	if (std::abs(S[0]) < eps && std::abs(S[1]) < eps ) {
		S[0] = eps; S[1] = eps;
	}

	if (S[1] < 0.0) { S[1] = -S[1]; }


	//Eigen::MatrixXf S0 = S;
	//mechanical_object.density_enery_term(S0);
	//Eigen::VectorXd grad;
	Eigen::VectorXd Sint = S;
	double fm;
	cost_function fcc(Element_Node_index, Global_Node_vector, mechanical_object, S_ref);
	//fm = fcc(Sint, grad);
	//std::cout << grad << std::endl;
	//std::cout << fm << std::endl;

	int niter = solver.minimize(fcc, Sint, fm);
	//std::cout << fm << std::endl;

	Mat2 matp = U * Sint.asDiagonal() * V.transpose();
	zi = Eigen::Map<Vec4>(matp.data());
}

//double HyperElastic_Triangle::cost_function( Eigen::VectorXd& S, Eigen::VectorXd& grad) {
//	Eigen::MatrixXd S0 = S;
//	mechanical_object.density_enery_term(S0);
//	grad = mechanical_object.gradient*area;
//	return mechanical_object.energy*area;
//}

cost_function::cost_function(const Vec3i& Element_Node_index, Eigen::MatrixXd& Global_Node_vector, MaterialModel& mechanical_prop, Vec2& S_ref0)
	:HTRI(Element_Node_index, Global_Node_vector, mechanical_prop) {
	S_ref = S_ref0;
}

double cost_function::operator ()(const Eigen::VectorXd& S, Eigen::VectorXd& grad, Eigen::MatrixXd& Hessian) {
		Eigen::MatrixXd S0 = S;
		HTRI.mechanical_object.density_enery_term(S0);

		double kw = HTRI.mechanical_object.bulk_modulus;
		double c1 = HTRI.mechanical_object.energy * HTRI.area;
		double c2 = 0.5 *kw* (S - S_ref).squaredNorm();  // need to multiplied bt weight 
		grad = HTRI.mechanical_object.gradient*HTRI.area + kw * (S - S_ref);
		Hessian = HTRI.mechanical_object.Hesssian * HTRI.area ;   // a term related to kw*(S-S_ref) should be added to hessian
		// Eigen::Ref<Eigen::Matrix2d> Hesssian0(Hessian);
		// Hessian = Hesssian0.inverse();
		//std::cout << Hessian << std::endl;
       return c1+c2;

}