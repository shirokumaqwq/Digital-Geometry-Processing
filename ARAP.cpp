#include "ARAP.h"
#include <iostream>
#include <Eigen/SVD>  
#include <Eigen/dense>

ARAP::ARAP(Mesh* mesh):Parameterization(mesh)
{
	
	//mesh_copy.assign(*mesh);
}

void ARAP::Recover()
{
	std::cout << "yuannum:" << mesh_copy.n_vertices() << "\n parenum:" << mesh->n_vertices() << std::endl;
	for (auto vh : mesh->vertices())
	{
		std::cout << mesh->point(vh) << " <- " << mesh_copy.point(mesh_copy.vertex_handle(vh.idx())) << std::endl;
		mesh->set_point(vh, mesh_copy.point(mesh_copy.vertex_handle(vh.idx())));
	}
}

void ARAP::InitLocal()
{
	if (mesh == nullptr)
		throw"mesh is null!";
	
	local_x.clear();
	local_x.resize(mesh->n_faces(),std::vector<Eigen::Vector2d>(3));
	local_cot.clear();
	local_cot.resize(mesh->n_faces(), std::vector<double>(3));

	std::vector<OpenMesh::Vec3d> temp3; temp3.resize(3);
	std::vector<Eigen::Vector2d> temp2; temp2.resize(3);
	std::vector<double> temptheta; temptheta.resize(3);
	for (auto& fh : mesh->faces())
	{
		int k = 0;
		int fidx = fh.idx();
		for (auto fhh : mesh->fh_range(fh))
		{
			temp3[k] = mesh->calc_edge_vector(fhh);
			k++;
		}
		temp2[0]={ temp3[0].norm(),0 };
		temp2[2]={ -temp3[0].dot(-temp3[2]) / temp3[0].norm(), -temp3[0].cross(-temp3[2]).norm() / temp3[0].norm() };
		temp2[1] = -temp2[2] - temp2[0];
		local_x[fidx] = temp2;
		//for (int k = 0; k < 3; k++)
			//std::cout << local_x[fidx][k] << " ---- " << temp3[k] << std::endl;
		//std::cout << std::endl;
		temptheta[0] = temp3[2].dot(-temp3[1]) / temp3[2].cross(temp3[1]).norm();
		temptheta[1] = temp3[0].dot(-temp3[2]) / temp3[0].cross(temp3[2]).norm();
		temptheta[2] = temp3[1].dot(-temp3[0]) / temp3[1].cross(temp3[0]).norm();
		
		local_cot[fidx] = temptheta;
		 //sd::cout << "fnum:" << fh.idx() << std::endl;
		//for (int k = 0; k < 3; k++)
			//std::cout <<"num:"<< fidx <<"   --"<< atan(1.0f / local_cot[fidx][k]) << std::endl;
		//std::cout << std::endl;
	}
	std::cout << "Init Local triangle coordinate\n" << std::endl;
}


std::vector<Eigen::Matrix2d> ARAP::Local()
{
	std::vector<Eigen::Vector2d> u_hedge;
	u_hedge.resize(3);


	std::vector<Eigen::Matrix2d> local_rotation;
	for (auto& xfh : mesh->faces())
	{
		Eigen::Matrix2d tempM = Eigen::Matrix2d::Zero();

		int fidx = xfh.idx();
		auto ufh = mesh->face_handle(fidx);

		int k = 0;
		for (auto ufhh : mesh->fh_range(ufh))
		{
			auto qwq = mesh->calc_edge_vector(ufhh);
			u_hedge[k] = { qwq[0],qwq[1] };
			k++;
		}
		for (int k = 0; k < 3; k++)
		{
			tempM += local_cot[fidx][k] * u_hedge[k] * local_x[fidx][k].transpose();
			//std::cout << local_cot[fidx][k] <<"\n" << u_hedge[k] << "\n" << local_x[fidx][k].transpose() << std::endl;
		}
		//Eigen::Matrix2d tempq, tempw;

		//tempq << -u_hedge[0].x(), u_hedge[1].x(),
			//	 -u_hedge[0].y(), u_hedge[1].y();
		//tempw << -local_x[fidx][0].x(), local_x[fidx][1].x(),
			//	 -local_x[fidx][0].y(), local_x[fidx][1].y();
		//std::cout << "Method1:\n"<<tempq * tempw.inverse() << std::endl;
		//Eigen::JacobiSVD<Eigen::MatrixXd> svd(tempq * tempw.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(tempM, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix2d U = svd.matrixU();
		Eigen::Matrix2d U1 = U;
		Eigen::Matrix2d V = svd.matrixV();
		if ((U * V.transpose()).determinant() < 0)
		{
			U1(0, 0) = U(1, 0); U1(0, 1) = U(1, 1);
			U1(1, 0) = U(0, 0); U1(1, 1) = U(0, 1);
		}
		local_rotation.push_back(U1*V.transpose());
		//std::cout << "Method2:\n"<<tempM << std::endl << U << std::endl << V << std::endl;;
		//std::cout << local_rotation.back() << std::endl << std::endl;
	}
	return  local_rotation;
}

void ARAP::Global(std::vector<Eigen::Matrix2d> local_rotation)
{
	int vnum = mesh->n_vertices();
	Eigen::SparseMatrix<double> L(vnum, vnum);
	Eigen::VectorXd Fu = Eigen::VectorXd::Zero(vnum);
	Eigen::VectorXd Fv = Eigen::VectorXd::Zero(vnum);


	for (auto& fh : mesh->faces())
	{
		int fidx = fh.idx();
		std::vector<int> vidx;
		for (auto vh : fh.vertices())
		{
			vidx.push_back(vh.idx());
			//std::cout << mesh->point(vh) << std::endl;
		}
		//std::cout << std::endl;

		for (int k = 0; k < 3; k++)
		{
			int kk = (k + 1) % 3;
			int k_ = (k +2) % 3;
			L.coeffRef(vidx[k], vidx[k]) += local_cot[fidx][kk];
			L.coeffRef(vidx[kk], vidx[kk]) += local_cot[fidx][kk];
			L.coeffRef(vidx[k], vidx[kk]) -= local_cot[fidx][kk];
			L.coeffRef(vidx[kk], vidx[k]) -= local_cot[fidx][kk];

			Fu(vidx[k]) += local_cot[fidx][k] * (local_rotation[fidx] * local_x[fidx][k])[0];
			Fv(vidx[k]) += local_cot[fidx][k] * (local_rotation[fidx] * local_x[fidx][k])[1];
			//std::cout << local_cot[fidx][k] << "\n" << local_rotation[fidx] << "\n" << local_x[fidx][k] << "\n\n";
			Fu(vidx[k_]) -= local_cot[fidx][k] * (local_rotation[fidx] * local_x[fidx][k])[0];
			Fv(vidx[k_]) -= local_cot[fidx][k] * (local_rotation[fidx] * local_x[fidx][k])[1];
		}

	}
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//std::cout <<"L:\n"<< L << std::endl;
	//std::cout << "Fu:\n" << Fu << std::endl;
	//std::cout << "Fv:\n" << Fv << std::endl;
	solver.compute(L);
	Eigen::VectorXd u_coor = solver.solve(Fu);
	Eigen::VectorXd v_coor = solver.solve(Fv);
	int k = 0;
	for (auto& vh : mesh->vertices())
	{		
		//std::cout << u_coor[vh.idx()] << "," << u_coor[vh.idx()] << " <- " << mesh->point(vh) << std::endl;
		mesh->set_point(vh, OpenMesh::Vec3d(u_coor[k], v_coor[k], 0));
		k++;
	}
}