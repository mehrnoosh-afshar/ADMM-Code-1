#include "admm_solver.h"
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include "Timer.h"
#include <omp.h>
using namespace ADMM;
using namespace Eigen;

typedef std::chrono::high_resolution_clock Clock;

void MultiplyVectorByScalar(std::vector<int>& v, int k) {
	for_each(v.begin(), v.end(), [k](int& c) { c *= k; });
}

Solver::Solver() : initialized(false) {
	//m_constraints = std::make_shared<ConstraintSet>(ConstraintSet());
}
void Solver::shape_global_vectors(Eigen::MatrixXd& Global_Node_vector, std::vector<float>& m) {

	int node_number = Global_Node_vector.rows();
	int size_Dof = node_number * problem_dimension;
	Eigen::MatrixXd x = Global_Node_vector.transpose();
	x.resize(size_Dof, 1);
	

	m_masses.resize(size_Dof, 1);
	m_x = x;
	m_v = m_x * 0;

	for (int i = 0; i < node_number; i++)
	{
		for (int j = 0; j < problem_dimension; j++)
		{
			m_masses[i * problem_dimension + j] = m[i];
		}
	}


}

void Solver::shape_external_force() {

	int n_force = controlpoints_inds.size();
	std::vector<int> initials;
	initials = controlpoints_inds;
	MultiplyVectorByScalar(initials, problem_dimension);
	ext_force = VecX::Zero(m_masses.size());
	for (int i = 0; i < n_force; ++i)
	{
			for (int j = 0; j < problem_dimension; j++)
			{
				ext_force[initials[i] + j] = force_value(i,j)/ m_masses[initials[i] + j];
			}
	}
	
}


void Solver::solver_step() {
	if (m_settings.verbose > 0) {
		std::cout << "\nSimulating with dt: " <<
			m_settings.timestep_s << "s..." << std::endl;
	}


	// Other const-variable short names and runtime data
	const int dof = m_x.rows();
	const int n_nodes = dof / problem_dimension;
	const double dt = m_settings.timestep_s;
	const int n_energyterms = energyterms.size();
	const int n_threads = std::min(n_energyterms, omp_get_max_threads());
	m_runtime = RuntimeData(); // reset so for each step we rewrite timer data
	Timer timeclock;

	// Add gravity
	if (std::abs(m_settings.gravity) > 0) {
		for (int i = 0; i < n_nodes; ++i) { m_v[i * 3 + 1] += dt * m_settings.gravity; }
	}

	// Position without elasticity/constraints 
	VecX x_bar = m_x + dt * m_v+ ext_force*std::pow(dt,2);
	VecX M_xbar = m_masses.asDiagonal() * x_bar;
	VecX curr_x = x_bar; // Temperorary x used in optimization

	// Initialize ADMM vars
	VecX curr_z = m_D * m_x;  // this  curr_z = m_D * m_x is upon convergence
	VecX curr_u = VecX::Zero(curr_z.rows());
	VecX solver_termB = VecX::Zero(dof); // it is a time varying matrix in run time
	VecX solver_termB_prime = VecX::Zero(dof+solver_term_D.size());
	// Passive collisions are detected each GS iteration,
	// so we'll skip them in the global collision detection loop.
	//bool detect_passive = m_settings.linsolver != 1;

	// Run a timestep
	int s_i = 0;
	for (; s_i < m_settings.admm_iters; ++s_i) {

		// Local step
		timeclock.Start(); // strat a time here
        #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
		for (int i = 0; i < n_energyterms; ++i) {
			energyterms[i]->update(m_D, curr_x, curr_z, curr_u);  
		}
		m_runtime.local_ms += timeclock.GetDuration();

		// Collision detection, which also updates the BVHs and stuff
		//t.reset();
		//m_constraints->collider->clear_hits();
		//m_constraints->collider->detect(surface_inds, curr_x, detect_passive);
		//m_runtime.collision_ms += t.elapsed_ms();

		// Global step
		timeclock.Start();
		solver_termB.noalias() = M_xbar + solver_Dt_Wt_W * (curr_z - curr_u);
		solver_termB_prime << solver_termB, solver_term_D; 
		
		// if using augmented A give curr_x0 to solver then extract curr_x from it
		VecX curr_x0;
		m_runtime.inner_iters += m_linsolver->solve(curr_x0, solver_termB_prime);
		curr_x = curr_x0.segment(0, x_bar.size());
		m_runtime.global_ms += timeclock.GetDuration(); // add a time here

	} // end solver loop

	// Computing new velocity and setting the new state
	m_v.noalias() = (curr_x - m_x) * (1.0 / dt);
	m_x = curr_x;

	// Output run time
	// if (m_settings.verbose > 0) { m_runtime.print(m_settings); } // check this for run time data 
} // end timestep iteration

bool Solver::initialize(const Settings& settings_) {
	
	m_settings = settings_;

	const int dof = m_x.rows();
	if (m_settings.verbose > 0) { std::cout << "Solver::initialize: " << std::endl; }

	if (m_settings.timestep_s <= 0.0) {
		std::cerr << "\n**Solver Error: timestep set to " << m_settings.timestep_s <<
			"s, changing to 1/24s." << std::endl;
		m_settings.timestep_s = 1.0 / 24.0;
	}
	if (!(m_masses.rows() == dof && dof >= 3)) {
		std::cerr << "\n**Solver Error: Problem with node data!" << std::endl;
		return false;
	}
	if (m_v.rows() != dof) { m_v.resize(dof); }

	// Clear previous runtime stuff settings
	m_v.setZero();

	// If we want energy-based constraints, set them up now.
	//if (m_settings.linsolver == 0 || m_settings.linsolver == 2) {
	//	std::unordered_map<int, Vec3>::iterator pinIter = m_constraints->pins.begin();
	//	for (; pinIter != m_constraints->pins.end(); ++pinIter) {
	//		m_pin_energies[pinIter->first] = std::make_shared<SpringPin>(SpringPin(pinIter->first, pinIter->second));
	//		energyterms.emplace_back(m_pin_energies[pinIter->first]);
	//	}
	// } // end create energy based hard constraints

	// Set up the selector matrix (D) and weight (W) matrix
	std::vector<Eigen::Triplet<double> > Dmatrix_triplets;
	std::vector<double> weights;
	Eigen::MatrixXd D_local_reduction;
	int n_energyterms = energyterms.size(); // this should be equal to the number of elements 
	for (int i = 0; i < n_energyterms; ++i) {
		//std::cout << energyterms[i]->Element_Node_index << std::endl;
		//std::cout << energyterms[i]->get_weight() << std::endl;
		energyterms[i]->get_local_Dmatrix(Dmatrix_triplets, D_local_reduction,weights);		
	}

	// Create the Selector+Reduction matrix
	m_W_diag = Eigen::Map<VecX>(&weights[0], weights.size());
	int n_D_rows = weights.size();
	m_D.resize(n_D_rows, dof);
	m_D.setZero();
	m_D.setFromTriplets(Dmatrix_triplets.begin(), Dmatrix_triplets.end());
    m_Dt = m_D.transpose();

	// Compute mass matrix
	SparseMat M(dof, dof);
	Eigen::VectorXi nnz = Eigen::VectorXi::Ones(dof); // non zeros per column
	M.reserve(nnz);
	for (int i = 0; i < dof; ++i) { M.coeffRef(i, i) = m_masses[i]; }


	// Set global matrices
	SparseMat W(n_D_rows, n_D_rows);
	W.reserve(n_D_rows);
	for (int i = 0; i < n_D_rows; ++i) { W.coeffRef(i, i) = m_W_diag[i]; }
	const double dt2 = (m_settings.timestep_s * m_settings.timestep_s);  
	solver_Dt_Wt_W = dt2 * m_Dt * W * W;
	solver_termA = M + SparseMat(solver_Dt_Wt_W * m_D);

	// Set global matrices for fixed boundary conditions  
	// this form is just duitable for completely fixed points 
	// if you want to extend this put it in another class 
    
	// return solver_termA to triplet to make the augmented triplet 
	// solver_agment_A =[A,CT;C,0]

	std::vector<Eigen::Triplet<double>> augm_teriplet;
	for (int i = 0; i < solver_termA.outerSize(); i++) {
		for (typename SparseMat::InnerIterator it(solver_termA, i); it; ++it)
			augm_teriplet.emplace_back(it.row(), it.col(), it.value());
	}


	int size_bc = BCpoints_inds.size() * problem_dimension;
	solver_term_D = VecX::Zero(size_bc);

	int j = 0;
	for (int i = 0; i < BCpoints_inds.size(); i++)
	{
		for (int k = 0; k < problem_dimension; k++)
		{
			solver_term_D(j) = m_x(BCpoints_inds[i] * problem_dimension + k);
			augm_teriplet.emplace_back(dof + j, BCpoints_inds[i] * problem_dimension + k, 1); // this adds values of C
			augm_teriplet.emplace_back(BCpoints_inds[i] * problem_dimension + k,dof + j , 1); // this adds values of CT
            //C.coeffRef(j, BCpoints_inds[i] * problem_dimension + k) = 1;
			++j;
		}
	}
	// Make augmented A and C matrix 
	solver_augment_termA.resize(dof+ size_bc, dof+ size_bc);
	solver_augment_termA.setZero();
	solver_augment_termA.setFromTriplets(augm_teriplet.begin(), augm_teriplet.end());

	// Set up the linear solver
     
	m_linsolver = std::make_shared<LDLTSolver>(LDLTSolver());
	

	// If we haven't set a global solver, make one:
	if (!m_linsolver) { throw std::runtime_error("What happened to the global solver?"); }
	// if (m_settings.constraint_w > 0.0) { m_constraints->constraint_w = m_settings.constraint_w; } // add constriants if you are doing collision detection 
	m_linsolver->update_system(solver_augment_termA);

	// Make sure they don't have any collision obstacles
	//if (m_settings.linsolver == 0) {
    //	if (m_constraints->collider->passive_objs.size() > 0 ||
	//		m_constraints->collider->dynamic_objs.size() > 0) {
	//		throw std::runtime_error("**Solver::add_obstacle Error: No collisions with LDLT solver");
	//	}
	//}

	// All done
	if (m_settings.verbose >= 1) { printf("%d nodes, %d energy terms\n", (int)m_x.size() / 3, (int)energyterms.size()); }
	initialized = true;
	return true;

} // end init

template<typename T> void myclamp(T& val, T min, T max) { if (val < min) { val = min; } if (val > max) { val = max; } }
bool Solver::Settings::parse_args(int argc, char** argv) {

	// Check args with params
	for (int i = 1; i < argc - 1; ++i) {
		std::string arg(argv[i]);
		std::stringstream val(argv[i + 1]);
		if (arg == "-help" || arg == "--help" || arg == "-h") { help(); return true; }
		else if (arg == "-dt") { val >> timestep_s; }
		else if (arg == "-v") { val >> verbose; }
		else if (arg == "-it") { val >> admm_iters; }
		else if (arg == "-g") { val >> gravity; }
		else if (arg == "-ls") { val >> linsolver; }
		else if (arg == "-ck") { val >> constraint_w; }
	}

	// Check if last arg is one of our no-param args
	std::string arg(argv[argc - 1]);
	if (arg == "-help" || arg == "--help" || arg == "-h") { help(); return true; }

	return false;

} // end parse settings args

void Solver::Settings::help() {
	std::stringstream ss;
	ss << "\n==========================================\nArgs:\n" <<
		"\t-dt: time step (s)\n" <<
		"\t-v: verbosity (higher -> show more)\n" <<
		"\t-it: # admm iters\n" <<
		"\t-g: gravity (m/s^2)\n" <<
		"\t-ls: linear solver (0=LDLT, 1=NCMCGS, 2=UzawaCG) \n" <<
		"\t-ck: constraint weights (-1 = auto) \n" <<
		"==========================================\n";
	printf("%s", ss.str().c_str());
}

void Solver::RuntimeData::print(const Settings& settings) {
	std::cout << "\nTotal global step: " << global_ms << "ms";;
	std::cout << "\nTotal local step: " << local_ms << "ms";
	std::cout << "\nTotal collision update: " << collision_ms << "ms";
	std::cout << "\nAvg global step: " << global_ms / double(settings.admm_iters) << "ms";;
	std::cout << "\nAvg local step: " << local_ms / double(settings.admm_iters) << "ms";
	std::cout << "\nAvg collision update: " << collision_ms / double(settings.admm_iters) << "ms";
	std::cout << "\nADMM Iters: " << settings.admm_iters;
	std::cout << "\nAvg Inner Iters: " << float(inner_iters) / float(settings.admm_iters);
	std::cout << std::endl;
}