#pragma once
#include <vector>
#include "MeshDefinition.h"

class Curvature {
public:
	Curvature();
	~Curvature();

	void MeanCurvature();
	void AbsoluteMeanCurvature();
	void GaussianCurvature();

	std::vector<OpenMesh::Vec3d> GetMeanCurvatureColor();
	std::vector<OpenMesh::Vec3d> GetAbsoluteMeanCurvatureColor();
	std::vector<OpenMesh::Vec3d> GetGaussianCurvatureColor();

	bool ifready();

	OpenMesh::Vec3d GetVertexColor(int n);
	void SetMesh(Mesh *imesh);

	OpenMesh::Vec3d Color(double V, double maxV, double minV);
private:
	int pnum;
	Mesh* mesh;

	std::vector<double> amean_curvature;
	std::vector<double> mean_curvature;
	std::vector<double> gaussian_curvature;

	std::vector<OpenMesh::Vec3d> ccolor;

	std::vector<double> LAR;  //local averaing region
};