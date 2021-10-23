#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>

class MeshParamWidget : public QWidget
{
	Q_OBJECT

public:
	MeshParamWidget(QWidget *parent = 0);
	~MeshParamWidget(void);
private:
	void CreateTabWidget(void);
	void CreateLayout(void);


signals:
	void PrintInfoSignal();
	void sShortestPath();
	void sMinTree();
	void sGaussianCurvature();
	void sMeanCurvature();
	void sBilateralNormalFilter();
	void sAddNoise();
	void sParameterization();
	void sARAP();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QPushButton* pbShortestPath;	
	QPushButton* pbMinTree;
	QPushButton* pbMeanCurvature;
	QPushButton* pbGaussianCurvature;
	QPushButton* pBilateralNormalFilter;
	QPushButton* pAddNoise;
	QPushButton* pParameterization;
	QPushButton* pARAP;
};
