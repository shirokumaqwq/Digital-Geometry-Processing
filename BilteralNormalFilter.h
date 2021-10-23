#pragma once
#include <vector>
#include "MeshDefinition.h"

class BilateralNormalFilter {
public:
	BilateralNormalFilter(Mesh* mesh, double sigma_r = 1.0,size_t n = 10);
	BilateralNormalFilter(Mesh* mesh, double sigma_r = 1.0);
	~BilateralNormalFilter();
	
	void Filtering();
	void SetN(size_t n);

private:
	double ComputeSigmac();
	std::vector<OpenMesh::Vec3d> ComputeNormal(std::vector<OpenMesh::Vec3d> last_normal);
	void UpdateVerticesPoistion(std::vector<OpenMesh::Vec3d> new_normal);
	double ComputeVolume();
	bool IsClosed();

	Mesh* mesh;
	double sigma_r2;
	double sigma_c2;
	size_t n=20;
	bool is_self_adaptive = false;
	double avg_edge_length;
	double qwq = 0;
};