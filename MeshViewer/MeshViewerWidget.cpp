#include <QtCore>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewerWidget.h"

#include <vector>
#include <time.h>

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();
	bool read_OK = MeshTools::ReadMesh(mesh, filename);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

void MeshViewerWidget::Clear(void)
{
	mesh.clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	mesh.update_normals();
	if (mesh.vertices_empty())
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : mesh.edges())
	{
		double len = mesh.calc_edge_length(eh);
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}
	
	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return MeshTools::WriteMesh(mesh, filename, DBL_DECIMAL_DIG);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (!mesh.vertices_empty())
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
}

void MeshViewerWidget::Dijkstra(std::vector<double>& dis,int n)
{
	//int n = v.idx();
	//保存已经求得到n的最小路径的顶点idx
	std::vector<int> S;
	S.push_back(n);
	//保存是否求过该顶点，是则1，否则0
	std::vector<int> ifp(mesh.n_vertices(), 0);
	ifp[n] = 1;
	dis[n] = 0;

	while (S.size()!=mesh.n_vertices())
	{

		double mindis = DBL_MAX;
		int minidx = -1;
		//遍历S求离S最近的点
		
		for (int i = 0; i < S.size(); i++)
		{
			if (S[i] >= 0)
			{
				Mesh::VertexHandle v0 = mesh.vertex_handle(S[i]);
				auto vh = mesh.vv_begin(mesh.vertex_handle(S[i]));

				//对S中的每一个点，求其到邻接点的距离
				for (; vh != mesh.vv_end(mesh.vertex_handle(S[i])); vh++)
				{
					Mesh::VertexHandle v1 = *vh;
					//std::cout <<S[i]<<"-"<< v1.idx() << "   ";
					if (ifp[v1.idx()] == 0)
					{
						
						double tempd = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
						//std::cout << S[i] << "的邻接点：" << v1.idx() << "-距离" << tempd;
						if (tempd < mindis)
						{
							mindis = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
							minidx = v1.idx();
							//std::cout << "<min,为 " << mindis << std::endl;
						}
					
							//std::cout << ">min"<<std::endl;
					}
				}
			}
		}
		S.push_back(minidx);
		//std::cout << "  " << S.back();
		if (minidx >= 0)
		{
			dis[minidx] = mindis;
			ifp[minidx] = 1;
		}
		//std::cout << std::endl << std::endl;
	}
	//std::cout << std::endl << std::endl << std::endl;
}

/*state:0为求最短路径，1为求n起始的最小生成树*/
void MeshViewerWidget::DijkstaBetwenTwoP(int n, int m)
{
	//保存已经求得到n的最小路径的顶点idx
	std::vector<int> S;
	S.push_back(n);
	//保存是否求过该顶点，是则1，否则0
	std::vector<int> ifp(mesh.n_vertices(), 0);
	ifp[n] = 1;
	std::vector<double> dis(mesh.n_vertices(), DBL_MAX);
	std::vector<int> preidx(mesh.n_vertices(), -1);
	dis[n] = 0;
	
	if (!SPath.empty())
		SPath.clear();
	isDrawSPath = true;
	isDrawMinTree = false;

	while (S.size() != mesh.n_vertices())
	{

		double mindis = DBL_MAX;
		int minidx = -1;
		int prei = -1;

		//遍历S求最近的点
		for (int i = 0; i < S.size(); i++)
		{
			if (S[i] >= 0)
			{
				Mesh::VertexHandle v0 = mesh.vertex_handle(S[i]);
				auto vh = mesh.vv_begin(mesh.vertex_handle(S[i]));

				//对S中的每一个点，求其到邻接点的距离
				for (; vh != mesh.vv_end(mesh.vertex_handle(S[i])); vh++)
				{
					Mesh::VertexHandle v1 = *vh;
					if (ifp[v1.idx()] == 0)
					{
						double tempd = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
						if (tempd < mindis)
						{
							mindis = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
							minidx = v1.idx();
							prei = v0.idx();
						}
					}
				}
			}
		}
		S.push_back(minidx);
		//std::cout << "  " << S.back();
		if (minidx >= 0)
		{
			dis[minidx] = mindis;
			ifp[minidx] = 1;
			preidx[minidx] = prei;
		}
		if (minidx == m)
		{
			//std::cout << m << "到" << n << "的最短路径为：" << std::endl << m;

			int i = m;
			while (preidx[i] != -1)
			{
				//std::cout << " - " << preidx[i];
				//寻找连接的边
				for (auto eh = mesh.voh_begin(mesh.vertex_handle(i)); eh != mesh.voh_end(mesh.vertex_handle(i)); eh++)
				{

					auto etemp = *eh;

					if (etemp.to().idx() == preidx[i])
					{
						SPath.push_back(etemp);
						break;
					}
				}
				i = preidx[i];
			}
			//std::cout<<std::endl;
			break;
		}
		
	}
}

void MeshViewerWidget::ShortestPath()
{
	if (mesh.vertices_empty())
	{
		std::cout << "Mesh is empty！" << std::endl;
		return;
	}
	std::cout << "Compute Shortest path..." << std::endl;
	//生成那个完全图
	//std::vector<std::vector<double>> edgeGraph(mesh.n_vertices(), std::vector<double>(mesh.n_vertices(), DBL_MAX));
	//保存该顶点是否被遍历过
	//std::vector<int> ipoints(mesh.n_vertices());
	/*size_t n = 0;
	for (const auto& v : mesh.vertices())
	{
		for (size_t i = 0; i < mesh.n_vertices(); i++)
			ipoints[i] = 0;
		
		Dijkstra(edgeGraph[n], n);

		//std::cout << n << ": " << v.idx() << " " << mesh.point(v)[0] << "," << mesh.point(v)[1] << "," << mesh.point(v)[2] << std::endl;
		n++;
	}
	end1 = clock();
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		for (int j = 0; j < mesh.n_vertices(); j++)
			std::cout << edgeGraph[i][j] << " ";
		std::cout << std::endl;
	}*/
	clock_t start, end;
	start = clock();
	
	DijkstaBetwenTwoP(0, 22);
	end = clock();
	std::cout << "生成两点距离用时：" << (end - start) << "ms" << std::endl << std::endl;
}

void MeshViewerWidget::MinSpanTree()
{
	clock_t t1, t2;
	t1 = clock();
	if (mesh.vertices_empty())
	{
		std::cout << "Mesh is empty！" << std::endl;
		return;
	}
	std::cout << "Compute Mininum span tree..." << std::endl;

	if (!MinTree.empty())
		MinTree.clear();
	isDrawMinTree = true;
	isDrawSPath = false;

	int n = 0;//求从n出发的最小生成树
	//保存已经求得到n的最小路径的顶点idx
	std::vector<int> S;
	S.push_back(n);
	//保存是否求过该顶点，是则1，否则0
	std::vector<int> ifp(mesh.n_vertices(), 0);
	ifp[n] = 1;
	std::vector<double> dis(mesh.n_vertices(), DBL_MAX);
	std::vector<int> preidx(mesh.n_vertices(), -1);
	dis[n] = 0;


	while (S.size() != mesh.n_vertices())
	{

		double mindis = DBL_MAX;
		int minidx = -1;
		int prei = -1;

		//遍历S求最近的点
		for (int i = 0; i < S.size(); i++)
		{
			if (S[i] >= 0)
			{
				Mesh::VertexHandle v0 = mesh.vertex_handle(S[i]);
				auto vh = mesh.vv_begin(mesh.vertex_handle(S[i]));

				//对S中的每一个点，求其到邻接点的距离
				for (; vh != mesh.vv_end(mesh.vertex_handle(S[i])); vh++)
				{
					Mesh::VertexHandle v1 = *vh;
					if (ifp[v1.idx()] == 0)
					{
						double tempd = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
						if (tempd < mindis)
						{
							mindis = (mesh.point(v1) - mesh.point(v0)).norm() + dis[v0.idx()];
							minidx = v1.idx();
							prei = v0.idx();
						}
					}
				}
			}
		}
		S.push_back(minidx);
		//std::cout << "  " << S.back();
		if (minidx >= 0)
		{
			dis[minidx] = mindis;
			ifp[minidx] = 1;
			preidx[minidx] = prei;
			
			MinTree.push_back(mesh.vertex_handle(prei));
			MinTree.push_back(mesh.vertex_handle(minidx));
			//std::cout << prei<<" - " << minidx << std::endl;
		}
	}
	t2 = clock();
	std::cout << "求最小生成树用时：" << t2 - t1 << "ms" << std::endl;
}

void MeshViewerWidget::DrawMeanCurvature()
{
	curvatureState = MEAN;
	curvature.SetMesh(&mesh);
	curvature.MeanCurvature();
	curvature.GetMeanCurvatureColor();
	//curvature.GetAbsoluteMeanCurvatureColor();
}

void MeshViewerWidget::DrawGaussianCurvature()
{
	curvatureState = GAUSSIAN;
	curvature.SetMesh(&mesh);
	curvature.GaussianCurvature();
	curvature.GetGaussianCurvatureColor();
}

void MeshViewerWidget::AddNoise()
{
	std::cout << "Add noise to mesh successfully!" << std::endl;

	for (auto& vh : mesh.vertices())
	{
		//仅选择30%的顶点添加一个沿顶点法向长度与相邻边平均边长相关的随机噪声
		if (!vh.is_boundary()&&rand()%10>6)
		{
			int k = 0; double avgl = 0;
			for (auto& eh : vh.edges())
			{
				avgl += mesh.calc_edge_length(eh);
				k++;
			}
			double noise = (rand() % 20 - 10) / 7;
			mesh.set_point(vh, mesh.point(vh) + mesh.calc_vertex_normal(vh) * avgl / k * noise);
		}
	}
	update();
}

void MeshViewerWidget::DoBilateralNormalFilter()
{
	std::cout << "Start filtering:" << std::endl;
	size_t n = 10;  //n为迭代次数
	BilateralNormalFilter filter(&mesh, 1, 20);
	filter.Filtering();
	update();
}

void MeshViewerWidget::DoParameterization()
{
	std::cout << "Tutte Parameterization" << std::endl;
	Parameterization parameter(&mesh);
	parameter.TutteParameterization();
	UpdateMesh();
	update();
}

void MeshViewerWidget::DoARAP()
{
	std::cout << "ARAP" << std::endl;
	ARAP parameter(&mesh);
	//parameter.TutteParameterization();
	parameter.InitLocal();
	parameter.TutteParameterization();
	for(int k=0;k<5;k++)
		parameter.Global(parameter.Local());

	UpdateMesh();
	update();
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (mesh.n_vertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	case CURVATURE:
		DrawCurvature();
	default:
		break;
	}

	if(isDrawSPath)
		DrawShortPath();
	if (isDrawMinTree)
		DrawMinTree();
}

void MeshViewerWidget::DrawPoints(void) const
{
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		glNormal3dv(mesh.normal(vh).data());
		glVertex3dv(mesh.point(vh).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : mesh.faces())
	{
		glNormal3dv(mesh.normal(fh).data());
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glVertex3dv(mesh.point(fvh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawCurvature(void) 
{
	//std::cout << "draw curvature" << std::endl;
	if (!curvature.ifready())
		return;
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	
	for (const auto& fh : mesh.faces())
	{
		//auto fnormal = mesh.calc_edge_vector(fh.halfedge().next()).cross(mesh.calc_edge_vector(fh.halfedge()));
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			
			auto tcolor = curvature.GetVertexColor(fvh.idx());
			glColor3d(tcolor[0], tcolor[1], tcolor[2]);
		
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::DrawSmooth(void) const
{
	glColor3d(0.8, 0.8, 0.8);
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	for (const auto& fh : mesh.faces())
	{
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			auto vh0 = mesh.from_vertex_handle(heh);
			auto vh1 = mesh.to_vertex_handle(heh);
			glNormal3dv(mesh.normal(vh0).data());
			glVertex3dv(mesh.point(vh0).data());
			glNormal3dv(mesh.normal(vh1).data());
			glVertex3dv(mesh.point(vh1).data());
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawShortPath()
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(5.0f);
	glColor3d(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < SPath.size(); i++)
	{
		OpenMesh::SmartHalfedgeHandle heh = SPath[i];
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		//std::cout << "画线：" << vh0.idx() << "-" << vh1.idx() << std::endl;
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawMinTree()
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(5.0f);
	glColor3d(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 0; i < MinTree.size(); i+= 2)
	{
		OpenMesh::VertexHandle vh0 = MinTree[i];
		OpenMesh::VertexHandle vh1 = MinTree[i + 1];
		//std::cout << "数组大小："<<MinTree.size()<<"画线：" << vh0.idx() << "-" << vh1.idx() << std::endl;
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}