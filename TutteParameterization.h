#pragma once

#include <vector>
#include "MeshDefinition.h"
#include<Eigen/Sparse>

typedef Eigen::Triplet<double> T;

class Parameterization {
public:
	Parameterization(Mesh* mesh);
	~Parameterization();

	void TutteParameterization();

private:
	void FindBoundary();
	std::vector<OpenMesh::Vec2d> CircleBoundaryCoordinate();

protected:
	std::vector<T> uniform_weight(std::vector<T>& tripletList);
	std::vector<T> cot_weight(std::vector<T>& tripletList);
	std::vector<int> Boundary_idx;
	Mesh* mesh;
};