#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "MeshDefinition.h"
#include "curvature.h"
#include "BilteralNormalFilter.h"
#include "TutteParameterization.h"
#include "ARAP.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
	void ShortestPath(void);
	void MinSpanTree(void);
	void DrawMeanCurvature(void);
	void DrawGaussianCurvature(void);
	void DoBilateralNormalFilter(void);
	void AddNoise(void);
	void DoParameterization(void);
	void DoARAP(void);

protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawSmooth(void) const;
	void DrawCurvature(void) ;
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
	void DrawShortPath(void);
	void DrawMinTree(void);
	void Dijkstra(std::vector<double> &dis, int n);
	void DijkstaBetwenTwoP(int n, int m);

protected:
	Mesh mesh;
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
	bool isDrawSPath = false;
	bool isDrawMinTree = false;
	enum CURVATURESTATE {DEFAULT, MEAN, GAUSSIAN };
	CURVATURESTATE curvatureState = DEFAULT;

	Curvature curvature;
	std::vector<OpenMesh::SmartHalfedgeHandle> SPath;
	std::vector<OpenMesh::VertexHandle> MinTree;
};
