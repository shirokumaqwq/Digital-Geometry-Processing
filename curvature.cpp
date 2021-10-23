#include "curvature.h"
#include <iostream>

Curvature::Curvature()
{

}

Curvature::~Curvature()
{

}

void Curvature::MeanCurvature()
{
	if (mesh == nullptr || LAR .empty())
		throw"mesh is null! or LAR is empty!";

	std::vector<OpenMesh::Vec3d> meanC;
	meanC.resize(pnum);
	for (size_t n = 0; n < meanC.size(); n++)
		meanC[n] = { 0,0,0 };
	std::cout << "初始化meanC数组大小:" << meanC.size() << std::endl;

	for (auto& vh : mesh->vertices())
	{
		if (vh.is_boundary())
			meanC[vh.idx()] = { 0,0,0 };
		else {
			for (auto& eh = vh.halfedge(); eh.opp().next() != vh.halfedge(); eh = eh.opp().next())
			{
				auto tempv1 = mesh->calc_edge_vector(eh.prev());
				auto tempv2 = mesh->calc_edge_vector(eh.next().opp());
				double cota = tempv1.dot(tempv2) / tempv1.cross(tempv2).norm();

				tempv1 = mesh->calc_edge_vector(eh.opp().prev());
				tempv2 = mesh->calc_edge_vector(eh.opp().next().opp());
				double cotb = tempv1.dot(tempv2) / tempv1.cross(tempv2).norm();

				meanC[vh.idx()] += (cota + cotb) * mesh->calc_edge_vector(eh.opp()) / 4 / LAR[vh.idx()];

			}
		}
	}
	
	mean_curvature.resize(pnum);
	for (size_t n = 0; n < pnum; n++)
	{
		
		mean_curvature[n] = mesh->calc_vertex_normal(mesh->vertex_handle(n)).dot(meanC[n]) > 0 ? meanC[n].norm() : -meanC[n].norm();
	}


}

void Curvature::AbsoluteMeanCurvature()
{
	amean_curvature.resize(pnum);
	for (size_t n = 0; n < pnum; n++)
	{
		amean_curvature[n] = abs(mean_curvature[n]);
	}

}

void Curvature::GaussianCurvature()
{
	if (mesh == nullptr || LAR.empty())
		throw"mesh is null!";

	gaussian_curvature.resize(pnum);
	for (size_t n = 0; n < pnum; n++)
	{
		gaussian_curvature[n] = 2 * 3.1415926;
		//std::cout << n << "'s gaussc:" << gaussian_curvature[n] << std::endl;

	}
	for (auto& vh : mesh->vertices())
	{
		if (vh.is_boundary())
			gaussian_curvature[vh.idx()] = 0;
		else {
			for (auto& eh : vh.outgoing_halfedges())
			{
				double theta = acos(mesh->calc_edge_vector(eh).normalize().dot(mesh->calc_edge_vector(eh.prev().opp()).normalize()));
				gaussian_curvature[vh.idx()] -= theta;
				
				//std::cout << theta << std::endl;
			}
			if (abs(gaussian_curvature[vh.idx()]) < 10e-6)
				gaussian_curvature[vh.idx()] = 0;
			gaussian_curvature[vh.idx()] = abs(gaussian_curvature[vh.idx()]) < 10e-4 ? 0 : gaussian_curvature[vh.idx()] / LAR[vh.idx()];

		}
		//std::cout << vh.idx() << "'s gaussc:" << gaussian_curvature[vh.idx()]<<"  is boundary?" <<vh.is_boundary()<< std::endl;
	}
}

std::vector<OpenMesh::Vec3d> Curvature::GetMeanCurvatureColor()
{
	if (mean_curvature.empty())
		MeanCurvature();

	ccolor.clear();
	ccolor.resize(pnum);

	double maxV = -1;
	double minV = INFINITY;

	for (size_t n = 0; n < pnum; n++)
	{
		maxV = maxV > mean_curvature[n] ? maxV : mean_curvature[n];
		minV = minV < mean_curvature[n] ? minV : mean_curvature[n];
	}

	for (size_t n = 0; n < pnum; n++)
	{
		std::cout << mean_curvature[n] << std::endl;
		ccolor[n] = Color(mean_curvature[n], maxV, minV);
	}

	return ccolor;
}

std::vector<OpenMesh::Vec3d> Curvature::GetAbsoluteMeanCurvatureColor()
{
	if (amean_curvature.empty())
		AbsoluteMeanCurvature();

	ccolor.clear();
	ccolor.resize(pnum);

	double maxV = -1;
	double minV = INFINITY;

	for (size_t n = 0; n < pnum; n++)
	{
		maxV = maxV > amean_curvature[n] ? maxV : amean_curvature[n];
		minV = minV < amean_curvature[n] ? minV : amean_curvature[n];
	}

	for (size_t n = 0; n < pnum; n++)
	{
		//std::cout << amean_curvature[n] << std::endl;
		ccolor[n] = Color(amean_curvature[n], maxV, minV);
	}
	//std::cout << "max:" << maxV << ", min:" << minV << std::endl;
	return ccolor;
}

std::vector<OpenMesh::Vec3d> Curvature::GetGaussianCurvatureColor()
{
	if (gaussian_curvature.empty())
		GaussianCurvature();

	ccolor.clear();
	ccolor.resize(pnum);

	double maxV = -1;
	double minV = INFINITY;

	for (size_t n = 0; n < pnum; n++)
	{
		
		maxV = maxV > gaussian_curvature[n] ? maxV : gaussian_curvature[n];
		minV = minV < gaussian_curvature[n] ? minV : gaussian_curvature[n];
	}

	for (size_t n = 0; n < pnum; n++)
	{
		//std::cout << "gaussc:" << gaussian_curvature[n] << std::endl;
		ccolor[n] = Color(gaussian_curvature[n], maxV, minV);
	}
	//std::cout << "max:" << maxV << ", min:" << minV << std::endl;
	return ccolor;
}

OpenMesh::Vec3d Curvature::GetVertexColor(int n)
{
	if (ccolor.size() < n || n < 0)
	{
		std::cout << "color out of range" << std::endl;
		return{ 0,0,0 };
	}

	return ccolor[n];
}

void Curvature::SetMesh(Mesh *imesh)
{
	mesh = imesh;
	pnum = imesh->n_vertices();
	LAR.resize(pnum);

	for (size_t n = 0; n < pnum; n++)
		LAR[n] = 0;
	//compute barycentric cell form LAR
	/*for (auto& fh : mesh->faces())
	{
		auto gc = mesh->calc_face_centroid(fh);
		for (auto& eh : fh.edges())
		{
			LAR[eh.v0().idx()] += (mesh->point(eh.v0()) - mesh->point(eh.v1())).cross(gc - mesh->point(eh.v0())).norm() / 4;
			LAR[eh.v1().idx()] += (mesh->point(eh.v0()) - mesh->point(eh.v1())).cross(gc - mesh->point(eh.v1())).norm() / 4;
		}
	}*/
	//mix方法
	for (auto& fh : mesh->faces())
	{
		//检查是否为钝角三角形
		bool isobtuse = false;
		int tempidx;
		for (auto& heh : fh.halfedges())
		{
			if (mesh->calc_edge_vector(heh).dot(mesh->calc_edge_vector(heh.opp().next())) <= 0)
			{
				isobtuse = true;
				tempidx = heh.from().idx();
				break;
			}
		}

		if (isobtuse)
		{
			std::cout << "obtuse:" << tempidx << std::endl;
			double area = mesh->calc_face_area(fh) / 4;
			for (auto& vh : fh.vertices())
			{
				if (vh.idx() == tempidx)
					LAR[vh.idx()] += area * 2;
				else
					LAR[vh.idx()] += area;
			}
		}
		else {
			//计算外心
			std::vector<unsigned int> vidx;
			std::vector<OpenMesh::Vec3d> vlist;
			for (auto& vh : fh.vertices())
			{
				vlist.push_back(mesh->point(vh));
				vidx.push_back(vh.idx());
			}double c0 = (vlist[1] - vlist[0]).dot(vlist[2] - vlist[0]);
			double c1 = (vlist[0] - vlist[1]).dot(vlist[2] - vlist[1]);
			double c2 = (vlist[0] - vlist[2]).dot(vlist[1] - vlist[2]);
			double c = (c0 + c1 + c2) * 2;
			OpenMesh::Vec3d ccenter = { (c1 + c2) / c,(c0 + c2) / c,(c0 + c1) / c };

			
			double s0 = (ccenter - vlist[0]).cross(vlist[1] - vlist[0]).norm() / 4;
			double s1 = (ccenter - vlist[1]).cross(vlist[2] - vlist[1]).norm() / 4;
			double s2 = (ccenter - vlist[2]).cross(vlist[0] - vlist[2]).norm() / 4;

			LAR[vidx[0]] += s0 + s2;
			LAR[vidx[1]] += s0 + s1;
			LAR[vidx[2]] += s1 + s2;
		}
	}
	
	for (size_t n = 0; n < pnum; n++)
		std::cout << "LAR"<<n<<":"<<LAR[n] << std::endl;
}

//rgb range 
OpenMesh::Vec3d Curvature::Color(double V,double maxV,double minV)
{
	if (maxV < minV)
		throw "maxV<minV!";
	V = (V - minV) / (maxV - minV);

	OpenMesh::Vec3d color;
	if (V < 0.25)
		color = { 0,4 * V,1 };
	else if (V < 0.5)
		color = { 0,1,1 - 4 * (V - 0.25) };
	else if (V < 0.75)
		color = { 4 * (V - 0.5), 1, 0 };
	else
		color = { 1,1 - 4 * (V - 0.75),0 };
	//std::cout << V << ":" << color << std::endl;
	return color;
}

bool Curvature::ifready()
{
	if (!ccolor.empty() && ccolor.size() == pnum)
		return true;
	else
		return false;
}