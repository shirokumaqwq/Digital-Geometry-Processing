#include "BilteralNormalFilter.h"
#include <iostream>
#include <time.h>

BilateralNormalFilter::BilateralNormalFilter(Mesh* mesh,  double sigma_r,size_t n) {
	this->mesh = mesh;
	this->n = n;
	this->is_self_adaptive = false;
	this->sigma_r2 = sigma_r * sigma_r;
	
	ComputeSigmac();
	
}
BilateralNormalFilter::BilateralNormalFilter(Mesh* mesh, double sigma_r)
{
	this->mesh = mesh;
	this->sigma_r2 = sigma_r * sigma_r;
	this->is_self_adaptive = true;
	ComputeSigmac();
}


BilateralNormalFilter::~BilateralNormalFilter() {
	mesh = nullptr;
	is_self_adaptive = false;
	qwq = 0;
}

void BilateralNormalFilter::SetN(size_t n)
{
	this->n = n;
}

double BilateralNormalFilter::ComputeSigmac()
{
	int n = 0;
	double d = 0;
	for (auto& fh : mesh->faces())
	{
		for (auto& ffh : fh.faces())
		{
			d += (mesh->calc_face_centroid(fh) - mesh->calc_face_centroid(ffh)).norm();
			n++;
		}
	}
	d /= n;
	this->sigma_c2 = d * d;
	this->avg_edge_length = d;
	std::cout << "sigmac^2:" << sigma_c2 << std::endl;

	return d;
}

std::vector<OpenMesh::Vec3d> BilateralNormalFilter::ComputeNormal(std::vector<OpenMesh::Vec3d> last_normal)
{
	double Kp = 0;
	std::vector<OpenMesh::Vec3d> new_normal(mesh->n_faces(), { 0,0,0 });
	//std::cout << "normalsize:" << new_normal.size() << std::endl;

	for (auto& fh : mesh->faces())
	{
		OpenMesh::Vec3d facenormal = { 0,0,0};
		for (auto& ffh : fh.faces())
		{
			double dc = (mesh->calc_face_centroid(fh) - mesh->calc_face_centroid(ffh)).norm();
			double dn = (last_normal[fh.idx()] - last_normal[ffh.idx()]).norm();
			double Ws = exp(-dc * dc / 2 / sigma_c2);
			double Wr = exp(-dn * dn / 2 / sigma_r2);
			double Area = mesh->calc_face_area(ffh);
			//Kp += Area * Ws * Wr;
			facenormal += Area * Ws * Wr * last_normal[ffh.idx()];
		}
		//facenormal /= Kp;
		new_normal[fh.idx()] = facenormal.normalize();
	}


	return new_normal;
}

void BilateralNormalFilter::UpdateVerticesPoistion(std::vector<OpenMesh::Vec3d> new_normal)
{
	qwq = 0;
	for (auto& vh : mesh->vertices())
	{
		OpenMesh::Vec3d dx = { 0,0,0 };
		int k = 0;
		for (auto& fvh : vh.faces())
		{
			OpenMesh::Vec3d tnormal = new_normal[fvh.idx()];
			dx += tnormal * (mesh->calc_face_centroid(fvh) - mesh->point(vh)).dot(tnormal);
			k++;
		}
		qwq = (dx / k).norm() > qwq ? (dx / k).norm() : qwq;
		//std::cout << mesh->point(vh) << " + " << dx /k << std::endl;
		mesh->set_point(vh, mesh->point(vh) + dx /k);
	}
	std::cout << "顶点移动比0.1平均边长:" << qwq / (0.1 * avg_edge_length) << std::endl;
}

void BilateralNormalFilter::Filtering()
{
	if (mesh == nullptr)
		throw"mesh is null!";

	clock_t start, end;
	start = clock();
	//初始化面法向
	std::vector<OpenMesh::Vec3d> last_normal(mesh->n_faces());
	for (auto& fh : mesh->faces())
		last_normal[fh.idx()] = mesh->calc_face_normal(fh);
	bool isclosed = IsClosed();
	double last_volume = 0.1;
	
	int k = 0;
	for (; k < n; k++)
	{
		if (isclosed)
			last_volume = ComputeVolume();
		last_normal = ComputeNormal(last_normal);
		UpdateVerticesPoistion(last_normal);
		if (qwq / (0.1 * avg_edge_length) < 1)
		{
			k++;
			break;
		}
		if (isclosed)
		{
			double pwp = ComputeVolume() / last_volume;
			//std::cout << k << "次体积比:" << qwq << std::endl;
			for (auto& vh : mesh->vertices())
				mesh->set_point(vh, mesh->point(vh) * pwp);
		}

	}
	
	mesh->update_face_normals();
	end = clock();
	std::cout << "循环" << k<< "次滤波，共用时" << end - start << "ms" << std::endl;
}

double BilateralNormalFilter::ComputeVolume()
{
	if (!IsClosed())
		throw"mesh is not closed!";
	OpenMesh::Vec3d origin = { 0,0,0 };
	double mesh_volume = 0;
	for (auto& fh : mesh->faces())
	{
		//实际是三倍体积
		mesh_volume += -mesh->calc_face_normal(fh).dot(origin - mesh->calc_face_centroid(fh)) * mesh->calc_face_area(fh);
	}
	return mesh_volume;
}

bool BilateralNormalFilter::IsClosed()
{
	if (mesh == nullptr)
		throw"mesh is null!";
	for (auto& vh : mesh->vertices())
	{
		if (vh.is_boundary())
			return false;
	}
	return true;
}