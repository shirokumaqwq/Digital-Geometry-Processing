#include "TutteParameterization.h"
#include <iostream>


Parameterization::Parameterization(Mesh* mesh)
{
	this->mesh = mesh;
}

Parameterization::~Parameterization()
{
	mesh = nullptr;
	Boundary_idx.clear();
}

void Parameterization::FindBoundary()
{
	if (mesh == nullptr)
		throw"mesh is null!";
	Boundary_idx.clear();
	OpenMesh::SmartVertexHandle v_b_start ;
	for (auto& vh : mesh->vertices())
	{
		if (vh.is_boundary())
		{
			v_b_start = vh;
			Boundary_idx.push_back(vh.idx());
			//std::cout << vh.idx() << std::endl;
			break;
		}
	}

	OpenMesh::SmartVertexHandle v_b = v_b_start;
	OpenMesh::SmartVertexHandle last_v = v_b;
	do
	{
		for (auto& hh : v_b.outgoing_halfedges())
			{
				if (hh.is_boundary() && hh.to().idx() != last_v.idx())
				{
					last_v = v_b;
					v_b = hh.to();
					Boundary_idx.push_back(v_b.idx());
					//std::cout << v_b.idx() << std::endl;
					break;
				}
			}
	} while (v_b.idx() != v_b_start.idx());
	Boundary_idx.pop_back();
}

std::vector<OpenMesh::Vec2d> Parameterization::CircleBoundaryCoordinate() 
{
	if (Boundary_idx.empty())
		FindBoundary();
	
	int v_num = Boundary_idx.size();
	std::vector<OpenMesh::Vec2d> boundary_coordinate(v_num);
	double pi = 3.1415926535;
	for (int n = 0; n < v_num; n++)
	{
		boundary_coordinate[n] = { cos(2 * pi / v_num * n),sin(2 * pi / v_num * n) };
	}
	return boundary_coordinate;
}

std::vector<T> Parameterization::uniform_weight(std::vector<T>& tripletList)
{
	for (auto& vh : mesh->vertices())
	{
		if (!vh.is_boundary())
		{
			int qwq = 0;
			for (auto& vvh : vh.vertices())
			{
				tripletList.push_back(T(vh.idx(), vvh.idx(), -1));
				qwq++;
			}
			tripletList.push_back(T(vh.idx(), vh.idx(), qwq));
		}
	}
	return tripletList;
}

std::vector<T> Parameterization::cot_weight(std::vector<T>& tripletList)
{
	for (auto& vh : mesh->vertices())
	{
		if (!vh.is_boundary())
		{
			double qwq = 0;
			for (auto& heh : vh.outgoing_halfedges())
			{
				//cos/sin
				auto v11 = mesh->calc_edge_vector(heh.prev());
				auto v12 = mesh->calc_edge_vector(heh.next().opp());
				double cota = v11.dot(v12) / v11.cross(v12).norm();

				auto v21 = mesh->calc_edge_vector(heh.opp().prev());
				auto v22 = mesh->calc_edge_vector(heh.opp().next().opp());
				double cotb = v21.dot(v22) / v21.cross(v22).norm();
				tripletList.push_back(T(vh.idx(), heh.to().idx(), -cota - cotb));
				qwq += cota + cotb;
			}
			tripletList.push_back(T(vh.idx(), vh.idx(), qwq));
		}
	}
	return tripletList;
}


void Parameterization::TutteParameterization()
{
	FindBoundary();
	std::vector<OpenMesh::Vec2d> boundary_coordinate = CircleBoundaryCoordinate();

	std::vector<T> tripletList;
	int vnum = mesh->n_vertices();
	Eigen::VectorXd U = Eigen::VectorXd::Zero(vnum);
	Eigen::VectorXd V = Eigen::VectorXd::Zero(vnum);

	for (int k=0;k<Boundary_idx.size();k++)
	{
		tripletList.push_back(T(Boundary_idx[k], Boundary_idx[k], 1));
		U[Boundary_idx[k]] = boundary_coordinate[k][0];
		V[Boundary_idx[k]] = boundary_coordinate[k][1];
	}

	//uniform_weight(tripletList);
	cot_weight(tripletList);
	
	Eigen::SparseMatrix<double> L(vnum, vnum);
	L.setFromTriplets(tripletList.begin(), tripletList.end());

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(L);
	Eigen::VectorXd u_coor = solver.solve(U);
	Eigen::VectorXd v_coor = solver.solve(V);

	for (int k = 0; k < vnum; k++)
		mesh->set_point(mesh->vertex_handle(k), OpenMesh::Vec3d(u_coor[k], v_coor[k], 0));
}