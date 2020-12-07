#include <Eigen/Sparse>
#include <Eigen/Dense>

#include"Triangle_Mesh.h"
#include"Triangle_Energy_Term.h"
#include <LBFGS.h>
//#include <chrono>
#include <omp.h>
#include "Timer.h"
//#include "optim.hpp"
#include"signedSVD.h"
#include"admm_solver.h"
#include "mex.h"


#define WIN32_LEAN_AND_MEAN
//#include <Windows.

using namespace Eigen;
bool myfunction(int i, int j) { return (i < j); }

Eigen::Matrix<double, Eigen::Dynamic, 1> AdmmWithTriangleMesh(Triangle_Mesh& trimesh) {
     MatrixXd Global_Node_vector;
   Eigen::MatrixXi Element_Node_index;

   Global_Node_vector = trimesh.Node;
   Element_Node_index = trimesh.Element;
   int Number_elemnets = Element_Node_index.rows();

   MaterialModel hyperelastic_materialmodel;
   const char* Type = "neo_hookean";
   hyperelastic_materialmodel.matrial_name_category(Type);

   std::vector<HyperElastic_Triangle> energyterms_list;

    ADMM::Solver solver_m;
   // solver setting  
    ADMM::Solver::Settings settings;
    settings.admm_iters = 100;
    settings.gravity = 0.0;
    settings.timestep_s = 0.001;
    std::vector<int> f (&trimesh.force_index(0), trimesh.force_index.data() + trimesh.force_index.cols() * trimesh.force_index.rows());
    std::vector<int> bc(&trimesh.Bc_index(0), trimesh.Bc_index.data() + trimesh.Bc_index.cols() * trimesh.Bc_index.rows());
    //std::vector<float> fvalue(&trimesh.force_value(0), trimesh.force_value.data() + trimesh.force_value.cols() * trimesh.force_value.rows());

    solver_m.controlpoints_inds = f;
    solver_m.BCpoints_inds = bc;
    solver_m.force_value = trimesh.force_value;
    solver_m.problem_dimension = 2;

    Vec3i Element_Node_index_sorted;

    for (int i = 0; i < Number_elemnets; i++)
    {
        Element_Node_index_sorted = Element_Node_index.row(i);
        //std::cout << Element_Node_index_sorted << std::endl;
        std::sort(Element_Node_index_sorted.data(), Element_Node_index_sorted.data() + Element_Node_index_sorted.size(), myfunction);
        HyperElastic_Triangle dummy(Element_Node_index_sorted, Global_Node_vector, hyperelastic_materialmodel);
        std::shared_ptr<HyperElastic_Triangle> pShared = std::make_shared<HyperElastic_Triangle>(dummy);
        // energyterms_list.push_back(HyperElastic_Triangle(Element_Node_index.row(i), Global_Node_vector, hyperelastic_materialmodel));
        solver_m.energyterms.push_back(pShared);
    }


    solver_m.shape_global_vectors(Global_Node_vector, trimesh.m);
    solver_m.shape_external_force();
    solver_m.initialize(settings);

    auto begin = std::chrono::high_resolution_clock::now();
    solver_m.solver_step();           // Hi this is the step that I need to be faster 
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << solver_m.m_x << std::endl;
    std::cout << "" << "time" << std::endl;
    std::cout << elapsed.count() * 1e-9 << "time" << std::endl;
    ADMM::Solver::RuntimeData R1;
    R1 = solver_m.runtime_data();
    return solver_m.m_x;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Triangle_Mesh trimesh;
    trimesh.load_mesh_data(prhs[3], prhs[4]);
    trimesh.load_point_indeces(prhs[2], prhs[0]);
    trimesh.load_forcevalue(prhs[1]);
    double density = 1000; // Kg/m
    trimesh.weighted_masses(density);
    auto mxResult = AdmmWithTriangleMesh(trimesh);
    std::cout << mxResult << std::endl;
    plhs[0] = mxCreateDoubleMatrix(mxResult.rows(),mxResult.cols(),mxREAL);
    Eigen::Map<Eigen::MatrixXd> map(mxGetPr(plhs[0]), mxResult.rows(), mxResult.cols());
    map = mxResult;
}

void admm(mxArray* forceIndex, mxArray* forceValue, mxArray* Bcindex, mxArray*Node_corr, mxArray* Element_index, int E, double nu) {
    Triangle_Mesh trimesh;
    trimesh.load_mesh_data(Node_corr, Element_index);
    trimesh.load_point_indeces(Bcindex, forceIndex);
    trimesh.load_forcevalue(forceValue);
    double density = 1000; // Kg/m
    trimesh.weighted_masses(density);
    AdmmWithTriangleMesh(trimesh);
}

int main()
{
    
    // 2D exmaple 
    Triangle_Mesh trimesh;
    trimesh.read_mesh_data("data/data.mat", "data/data2.mat");
    trimesh.read_point_indeces("data/Bcindex.mat", "data/Forceindex.mat");
    trimesh.read_forcevalue("data/ForceValue.mat");
    double density = 1000; // Kg/m
    trimesh.weighted_masses(density);
    AdmmWithTriangleMesh(trimesh);
    // std::cout << trimesh.Element << std::endl; 
   // std::cout << trimesh.force_value << std::endl;
   // for (size_t i = 0; i < trimesh.m.size(); i++)
   // {
   //     std::cout << trimesh.m[i] << std::endl;
   // }

   
  
    //std::cout << R1.local_ms << "time" << std::endl;
    //std::cout << solver_m.m_x[200]  << std::endl;
    //std::cout << solver_m.m_x[100] << std::endl;
    //std::cout << solver_m.m_x[70] << std::endl;
    //std::cout << solver_m.m_x[3] << std::endl;
    //std::cout << solver_m.m_x[4] << std::endl;
    //std::cout << solver_m.m_x[48] << std::endl;
}