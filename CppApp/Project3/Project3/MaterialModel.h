#pragma once
#include<iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <math.h> 

	static enum StringValue {
		linearelstic,
		neo_hookean,
		odegen
	};

	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;

	class MaterialModel
	{

	protected:
		double E; // module of elasticity 
		double nu; // possion ration


	public:
		double bulk_modulus; // bulk moduless
		std::string  material_name;
		double energy; // say is strain energy density term
		VecX gradient;
		Eigen::MatrixXd Hesssian; 
		inline void set_meachanicalproperties(double e, double nu0);
		inline void matrial_name_category(const char* category);
		inline void density_enery_term(Eigen::MatrixXd& S);

		MaterialModel() {
				material_name = "linear";
				E = 3000;
				nu = 0.49;
				bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));
		};
		// ~MaterialModel();





	};

	inline void MaterialModel::set_meachanicalproperties(double e, double nu0)
	{
		E = e; nu = nu0;
	}
	inline void MaterialModel::matrial_name_category(const char* category)
	{

		material_name = (char*)category;
	}
	inline StringValue hashit(std::string const& material_name) {
		if (material_name == "linearelstic") return linearelstic;
		if (material_name == "neo_hookean") return neo_hookean;
		if (material_name == "odegen") return odegen;
	}
	inline void MaterialModel::density_enery_term(Eigen::MatrixXd& S)
	{
		double k = E / (3.0 * (1.0 - 2.0 * nu));
		bulk_modulus = k;
		double miu = E / (2.0 * (1.0 + nu));

		double I, J, I_bar;
		I = 0; J = 0; I_bar = 0;
		if (S.size() == 2)
		{
			 J = S(0) * S(1);
			 I = std::pow(S(0), 2) + std::pow(S(1), 2)+1; // the one is to compensate for S(2)=1 for plane strain case
			 I_bar = std::pow(J, -(2.0 / 3.0)) * I;
		}
		else
		{
			 J = S(0) * S(1) * S(2);
			 I = std::pow(S(0), 2) + std::pow(S(1), 2) + std::pow(S(2), 2);
			 I_bar = std::pow(J, -2.0 / 3.0) * I;
		}
		switch (hashit(material_name))
		{
		case linearelstic:
			energy = 1;
			break;
		case neo_hookean:
			energy = 0.5 * miu * (I_bar - 3) + 0.5 * k * std::pow(J - 1, 2);
			if (S.size() == 2)
			{
				gradient.resize(2, 1);
				gradient << 0.5f * miu * (2 * S(0) * std::pow(J, -2.0f / 3.0f) - (2.0f / 3.0f) * S(1) * std::pow(J, -5.0f / 3.0f) * I) + k * S(1) * (J - 1),
					0.5f * miu * (2 * S(1) * std::pow(J, -2.0f / 3.0f) - (2.0f / 3.0f) * S(0) * std::pow(J, -5.0f / 3.0f) * I) +  k * S(0) * (J - 1);
				Hesssian.resize(2, 2);
				Hesssian << 0.5f * miu * (-0.6667f * pow(J, -2.0 / 3.2) + (10.0 / 9.0) * I * std::pow(J, -8.0 / 3.0) * S(1) * S(1)) + k * S(1) * S(1), 2 * k * J - k - (4.0 / 9.0) * miu * pow(J, -5.0 / 3.2) * I,
					2 * k * J - k - (4.0 / 9.0) * miu * pow(J, -5.0 / 3.2) * I, 0.5f * miu * (-0.6667f * pow(J, -2.0 / 3.2) + (10.0 / 9.0) * I * std::pow(J, -8.0 / 3.0) * S(0) * S(0)) + k * S(0) * S(0);
				
			}
			else
			{
				gradient.resize(3, 1);
				gradient << 0.5 * miu * (2 * S(0) * std::pow(J, -2 / 3) - (2 / 3) * S(1) *S(2)* std::pow(J, -5 / 3) * I) +  k * S(1) *S(2)* (J - 1),
					0.5 * miu * (2 * S(1) * std::pow(J, -2 / 3) - (2 / 3) * S(0) *S(2)* std::pow(J, -5 / 3) * I) +  k * S(0) *S(2)* (J - 1),
					0.5 * miu * (2 * S(2) * std::pow(J, -2 / 3) - (2 / 3) * S(0) * S(1) * std::pow(J, -5 / 3) * I) +  k * S(0) * S(1) * (J - 1);
				Hesssian.resize(3, 3);
				Hesssian << 0.5f * miu * (-0.6667f * pow(J, -2.0 / 3.2) + (10.0 / 9.0) * I * std::pow(J, -8.0 / 3.0) * S(1) * S(1) * S(2) * S(2)) + k * S(1) * S(1) * S(2) * S(2),
					2 * k * J * S(2) - k * S(2) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(0) * S(0) * S(2) + (4.0 / 3.0) * S(1) * S(1) * S(2) + (4.0 / 3.0) * I * S(2) - (10.0 / 9.0) * I * S(2)),
					2 * k * J * S(1) - k * S(1) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(0) * S(0) * S(1) + (4.0 / 3.0) * S(1) * S(2) * S(2) + (4.0 / 3.0) * I * S(1) - (10.0 / 9.0) * I * S(1)),
					2 * k * J * S(2) - k * S(2) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(0) * S(0) * S(2) + (4.0 / 3.0) * S(1) * S(1) * S(2) + (4.0 / 3.0) * I * S(2) - (10.0 / 9.0) * I * S(2)),
					0.5f * miu * (-0.6667f * pow(J, -2.0 / 3.2) + (10.0 / 9.0) * I * std::pow(J, -8.0 / 3.0) * S(0) * S(0) * S(2) * S(2)) + k * S(0) * S(0) * S(2) * S(2),
					2 * k * J * S(0) - k * S(0) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(1) * S(1) * S(0) + (4.0 / 3.0) * S(0) * S(2) * S(2) + (4.0 / 3.0) * I * S(0) - (10.0 / 9.0) * I * S(0)),
					2 * k * J * S(1) - k * S(1) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(0) * S(0) * S(1) + (4.0 / 3.0) * S(1) * S(2) * S(2) + (4.0 / 3.0) * I * S(1) - (10.0 / 9.0) * I * S(1)),
					2 * k * J * S(0) - k * S(0) - 0.5 * miu * std::pow(J, -5.0 / 3.0) * ((4.0 / 3.0) * S(1) * S(1) * S(0) + (4.0 / 3.0) * S(0) * S(2) * S(2) + (4.0 / 3.0) * I * S(0) - (10.0 / 9.0) * I * S(0)),
					0.5f * miu * (-0.6667f * pow(J, -2.0 / 3.2) + (10.0 / 9.0) * I * std::pow(J, -8.0 / 3.0) * S(0) * S(0) * S(1) * S(1)) + k * S(0) * S(0) * S(1) * S(1);

			}
			
			
				      
			break;
		case odegen:
			energy = 1;
			break;

		}
	}

//	 MaterialModel::MaterialModel()
//	{
//		material_name = "linear";
//		E = 3000;
//		nu = 0.49;
//		bulk_modulus = E / (3 * (1 - 2 * nu));
//
//	}


//MaterialModel::~MaterialModel()
//{
//}

