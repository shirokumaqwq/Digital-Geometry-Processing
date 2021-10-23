#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Core>

#include "TutteParameterization.h"

typedef Eigen::Triplet<double> T;

class ARAP :public Parameterization {
public:
	ARAP(Mesh* mesh);
	
	~ARAP() {
		mesh_copy.clear();
		local_x.clear();
		local_cot.clear();
	}
	void Recover();	
	void InitLocal();
	std::vector<Eigen::Matrix2d> Local();
	void Global(std::vector<Eigen::Matrix2d> local_rotation);
protected:
	Mesh mesh_copy;
	
	std::vector<std::vector<Eigen::Vector2d>> local_x;
	std::vector<std::vector<double>> local_cot;

};