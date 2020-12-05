#pragma once
#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H 1

#include <Eigen/Sparse>
#include "Triangle_Energy_Term.h"
#include "LinearSolver.h"





namespace ADMM {

	// The main solver
	class Solver {
	public:

		Solver();  // default constructor 

		typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
		typedef Eigen::Matrix<double, 3, 1> Vec3;
		typedef Eigen::Matrix<double, 2, 1> Vec2;
		typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;

		// Solver settings
		struct Settings {
			bool parse_args(int argc, char** argv); // parse from terminal args. Returns true if help()
			void help();		// -help	print details, parse_args returns true if used
			double timestep_s;	// -dt <flt>	timestep in seconds
			int verbose;		// -v <int>	terminal output level (higher=more)
			int admm_iters;		// -it <int>	number of admm solver iterations
			double gravity;		// -g <flt>	force of (-y) gravity . if you do not want to have gravity put it equal to 0 
			int linsolver;		// -ls <int>	0=LDLT, 1=NCMCGS, 2=UzawaCG
			double constraint_w;	// -ck <flt>	constraint weights (-1 = auto)
			Settings() : timestep_s(1.0 / 24.0), verbose(1), admm_iters(10),
				gravity(-9.8), linsolver(0), constraint_w(-1) {}
		};

		// RuntimeData struct used for logging.
		// Add timings are per time step.
		struct RuntimeData {
			double global_ms; // total ms for global solver
			double local_ms; // total ms for local solver
			double collision_ms; // total ms for collision update/detection
			int inner_iters; // total global step iterations
			RuntimeData() : global_ms(0), local_ms(0), collision_ms(0), inner_iters(0) {}
			void print(const Settings& settings);
		};

		

		// Per-node (x3) data (for x, y, and z)
		int problem_dimension;
		VecX m_x; // node positions, scaled x3 for 3D case and scaled x2 for 2D case
		VecX m_v; // node velocities, scaled x3 for 3D case and scaled x2 for 2D case
		VecX m_masses; // node masses, scaled x3 for 3D case and scaled x2 for 2D case
		VecX ext_force;  // external forces applied to Dofs
		std::vector<int> surface_inds; // indices of surface vertices
		std::vector<int> controlpoints_inds; // indices of control points vertices _ the indices should start from zeo
		// std::vector<float>  force_value; // the value of each row (force sets) in controlpoints_inds-_ the indices should start from zeo
		Eigen::MatrixXd force_value;;
		std::vector<int> BCpoints_inds;  //indices of BC points  vertices_ the indices should start from zeo

		std::vector< std::shared_ptr< HyperElastic_Triangle> > energyterms; // minimized (implicit) this vector of pointer is where the physical data and
		                                                                    // mesh data is integrated into Solver. 

		// Adds nodes to the Solver.
		// Returns the current total number of nodes after insert.
		// Assumes m is scaled x3 (i.e. 3 values per node).

		void Solver::shape_global_vectors(Eigen::MatrixXd& Global_Node_vector, std::vector<float>& m);
		void Solver::shape_external_force();
		// Returns true on success.
		virtual bool initialize(const Settings& settings_ = Settings());

		// Performs a Solver step // both local step and global step is in there
		virtual void solver_step();

		// Returns the runtime data from the last time step.
		virtual const RuntimeData& runtime_data() { return m_runtime; }

		// Outputs solver_termA to a file (for debugging/analysis).
		//virtual void save_matrix(const std::string& filename);

		// Returns the current settings
		const Settings& settings() { return m_settings; }

	protected:

		Settings m_settings; // copied from init
		RuntimeData m_runtime; // reset each iteration
		bool initialized; // has init been called?

		// Solver used in the global step
		std::shared_ptr<LinearSolver> m_linsolver; // linear solver and constraints are for global step
		//std::shared_ptr<ConstraintSet> m_constraints;

		// Global matrices
		SparseMat m_D, m_Dt, C; // reduction matrix
		VecX m_W_diag , solver_term_D; // diagonal of the weight matrix

		// Solver variables computed in initialize
		SparseMat solver_termA;
		SparseMat solver_Dt_Wt_W;
		SparseMat solver_augment_termA;
	}; // end class Solver

} // end namespace admm

#endif