#pragma once
#include "mat.h"
#include "matrix.h"
#include<iostream>
#include <vector>
#include <stdio.h>
#include <Eigen/Dense>
#include <string>

typedef Eigen::Matrix<int, 1, 3> Vec3i;
typedef Eigen::Matrix<double, 1, 2> Vec2d;

class Triangle_Mesh
{
public:
	
	Eigen::MatrixXd Node;
	Eigen::MatrixXi Element;
	Eigen::MatrixXd force_value;
	Eigen::MatrixXi Bc_index;
	Eigen::MatrixXi force_index;

	std::vector<float> m;
	void read_mesh_data(const char *Node_path , const char *Element_path);
	void weighted_masses(float density_kgm);
	void read_point_indeces(const char* Bcindeces_path, const char* Forceindeces_path);
	void read_forcevalue(const char* forcevalue_vec_path);
	Triangle_Mesh();
//	~Triangle_Mesh();

private:


};

void Triangle_Mesh::read_mesh_data(const char *Node_path, const char *Element_path)
{
	MATFile* pmat_node;
	MATFile* pmat_element;
	mxArray* pMxArray_node; // <- Check this, I think this is missing something
	mxArray* pMxArray_element;


	//OutputDebugStringW(L"Here1\n");
	pmat_node = matOpen(Node_path, "r");
	pmat_element = matOpen(Element_path, "r");

	if (pmat_node == NULL) {
		std::cout << "no file for mesh's node";
	}
	if (pmat_element == NULL) {
		std::cout << "no file for mesh's element";
	}

	pMxArray_node = matGetVariable(pmat_node, "Nodes");
	pMxArray_element = matGetVariable(pmat_element, "Elements");
	
	if (pMxArray_node == NULL)
	{
		std::cout << "Error reading existing matrix for Nodes !!!" << std::endl;
	}
	if (pMxArray_element == NULL)
	{
		std::cout << "Error reading existing matrix for Elements !!!" << std::endl;
	}

	double* pr_node = mxGetPr(pMxArray_node);
	double* pr_element = mxGetPr(pMxArray_element);

	size_t Mnode = mxGetM(pMxArray_node);
	size_t Nnode = mxGetN(pMxArray_node);

	size_t Melement = mxGetM(pMxArray_element);
	size_t Nelement = mxGetN(pMxArray_element);

	Node.resize(Mnode, Nnode);
	Element.resize(Melement, Nelement);

	for (int i = 0; i < Nnode; i++)
	{
		for (int j = 0; j < Mnode; j++)
		{
			Node(j, i) = *(pr_node + Mnode * i + j);
		}
	}

	for (int i = 0; i < Nelement; i++)
	{
		for (int j = 0; j < Melement; j++)
		{
			Element(j, i) = *(pr_element + Melement * i + j);
		}
	}

	matClose(pmat_node);
	matClose(pmat_element);
	//OutputDebugStringW(L"Here2\n");

}

void Triangle_Mesh::read_point_indeces(const char* Bcindeces_path, const char* Forceindeces_path)
{
	MATFile* pmat_bc;
	MATFile* pmat_fc;
	mxArray* pMxArray_bc; // <- Check this, I think this is missing something
	mxArray* pMxArray_fc;


	//OutputDebugStringW(L"Here1\n");
	pmat_bc = matOpen(Bcindeces_path, "r");
	pmat_fc = matOpen(Forceindeces_path, "r");

	if (pmat_bc == NULL) {
		std::cout << "no file for mesh's node";
	}
	if (pmat_fc == NULL) {
		std::cout << "no file for mesh's element";
	}

	pMxArray_bc = matGetVariable(pmat_bc, "NBC");
	pMxArray_fc = matGetVariable(pmat_fc, "NF");

	if (pMxArray_bc == NULL)
	{
		std::cout << "Error reading existing matrix for Nodes !!!" << std::endl;
	}
	if (pMxArray_fc == NULL)
	{
		std::cout << "Error reading existing matrix for Elements !!!" << std::endl;
	}

	double* pr_bc = mxGetPr(pMxArray_bc);
	double* pr_force = mxGetPr(pMxArray_fc);

	size_t MBC= mxGetM(pMxArray_bc);
	size_t NBC = mxGetN(pMxArray_bc);

	size_t MNF = mxGetM(pMxArray_fc);
	size_t NNF = mxGetN(pMxArray_fc);

	Bc_index.resize(MBC, NBC);
	force_index.resize(MNF, NNF);

	for (int i = 0; i < NBC; i++)
	{
		for (int j = 0; j < MBC; j++)
		{
			Bc_index(j, i) = *(pr_bc + MBC * i + j)-1; // this -1 is beacuse Matlab indices strat from 1
		}
	}

	for (int i = 0; i < NNF; i++)
	{
		for (int j = 0; j < MNF; j++)
		{
			force_index(j, i) = *(pr_force + MNF * i + j)-1; // this -1 is beacuse Matlab indices strat from 1
		}
	}

	matClose(pmat_bc);
	matClose(pmat_fc);
	//OutputDebugStringW(L"Here2\n");

}

void Triangle_Mesh::read_forcevalue(const char* forcevalue_vec_path)
{
	MATFile* pmat_fvalue;
	mxArray* pMxArray_fvalue; // <- Check this, I think this is missing something


	//OutputDebugStringW(L"Here1\n");
	pmat_fvalue = matOpen(forcevalue_vec_path, "r");

	if (pmat_fvalue == NULL) {
		std::cout << "no file for mesh's node";
	}
	

	pMxArray_fvalue = matGetVariable(pmat_fvalue, "Force");

	if (pMxArray_fvalue == NULL)
	{
		std::cout << "Error reading existing matrix for Nodes !!!" << std::endl;
	}

	double* pr_fvalue = mxGetPr(pMxArray_fvalue);

	size_t Mfvalue = mxGetM(pMxArray_fvalue);
	size_t Nfvalue = mxGetN(pMxArray_fvalue);



	force_value.resize(Mfvalue, Nfvalue);

	for (int i = 0; i < Nfvalue; i++)
	{
		for (int j = 0; j < Mfvalue; j++)
		{
			force_value(j, i) = *(pr_fvalue + Mfvalue * i + j);
		}
	}


	matClose(pmat_fvalue);
	//OutputDebugStringW(L"Here2\n");

}



void Triangle_Mesh::weighted_masses(float density_kgm) 
{
	m.resize(Node.rows());


	int num_nodes = Node.rows();
	int num_element = Element.rows();
	for (int f = 0; f < num_element; ++f) {
		//std::cout << trimesh.Element.row(f) << std::endl;
		Vec3i nodes_index = Element.row(f);
		Eigen::Matrix3d D;
		D.row(0) << Node.row(nodes_index[0] - 1).col(0), Node.row(nodes_index[0] - 1).col(1), 1;
		D.row(1) << Node.row(nodes_index[1] - 1).col(0), Node.row(nodes_index[1] - 1).col(1), 1;
		D.row(2) << Node.row(nodes_index[2] - 1).col(0), Node.row(nodes_index[2] - 1).col(1), 1;
		//std::cout << f << std::endl;
		//Vec2d edge1 = Node.row(nodes_index[1]) - Node.row(nodes_index[0]);
		//Vec2d edge2 = Node.row(nodes_index[2]) - Node.row(nodes_index[0]);
		//float area = 0.5f * (edge1.cross(edge2)).norm();
		float area = 0.5f * D.determinant();
		float tri_mass = density_kgm * area;

		m[nodes_index[0] - 1] += tri_mass / 3.f;
		m[nodes_index[1] - 1] += tri_mass / 3.f;
		m[nodes_index[2] - 1] += tri_mass / 3.f;

	}
}


Triangle_Mesh::Triangle_Mesh()
{
}

//Triangle_Mesh::~Triangle_Mesh()
//{
//}