#pragma once
#include "mat.h"
#include "matrix.h"
#include<iostream>
#include <vector>
#include <stdio.h>
#include <Eigen/Dense>
#include <string>

typedef Eigen::Matrix<int, 1, 4> Vec4i;


class Tet_Mesh
{
public:

	Eigen::MatrixXd Node;
	Eigen::MatrixXi Element;
	std::vector<float> m;
	void read_mesh_data(const char* Node_path, const char* Element_path);
	void weighted_masses(float density_kgm);
	Tet_Mesh();
	//	~Triangle_Mesh();

private:


};

void Tet_Mesh::read_mesh_data(const char* Node_path, const char* Element_path)
{
	MATFile* pmat_node;
	MATFile* pmat_element;
	mxArray* pMxArray_node;
	mxArray* pMxArray_element;

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

}

void Tet_Mesh::weighted_masses(float density_kgm)
{
	m.resize(Node.rows());


	int num_nodes = Node.rows();
	int num_element = Element.rows();
	for (int f = 0; f < num_element; ++f) {
		Vec4i nodes_index = Element.row(f);
		Eigen::Matrix<float, 3, 3> edges;
		edges.col(0) = (Node.row(nodes_index[1]-1) - Node.row(nodes_index[0] - 1)).transpose();
		edges.col(1) = (Node.row(nodes_index[2] - 1) - Node.row(nodes_index[0] - 1)).transpose();
		edges.col(2) = (Node.row(nodes_index[3] - 1) - Node.row(nodes_index[0] - 1)).transpose();
		float v = std::abs((edges).determinant() / 6.f);
		float tet_mass = density_kgm * v;
		m[nodes_index[0] - 1] += tet_mass / 4.f;
		m[nodes_index[1] - 1] += tet_mass / 4.f;
		m[nodes_index[2] - 1] += tet_mass / 4.f;
		m[nodes_index[3] - 1] += tet_mass / 4.f;

	}
}


Tet_Mesh::Tet_Mesh()
{
}

//Triangle_Mesh::~Triangle_Mesh()
//{
//}