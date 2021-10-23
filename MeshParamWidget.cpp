#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));
	pbShortestPath = new QPushButton(tr("Find Shortest Path"));
	connect(pbShortestPath, SIGNAL(clicked()), SIGNAL(sShortestPath()));
	pbMinTree = new QPushButton(tr("Find Minimum Spanning Tree"));
	connect(pbMinTree, SIGNAL(clicked()), SIGNAL(sMinTree()));
	pbMeanCurvature = new QPushButton(tr("Compute Mean Curvature"));
	connect(pbMeanCurvature, SIGNAL(clicked()), SIGNAL(sMeanCurvature()));
	pbGaussianCurvature= new QPushButton(tr("Compute Gaussian Curvature"));
	connect(pbGaussianCurvature, SIGNAL(clicked()), SIGNAL(sGaussianCurvature()));
	pBilateralNormalFilter = new QPushButton(tr("Do Blateral Normal Filtering"));
	connect(pBilateralNormalFilter, SIGNAL(clicked()), SIGNAL(sBilateralNormalFilter()));
	pAddNoise = new QPushButton(tr("Add Noise"));
	connect(pAddNoise, SIGNAL(clicked()), SIGNAL(sAddNoise()));
	pParameterization = new QPushButton(tr("Do Tutte Parameterization"));
	connect(pParameterization, SIGNAL(clicked()), SIGNAL(sParameterization()));
	pARAP = new QPushButton(tr("Do ARAP Parameterization"));
	connect(pARAP, SIGNAL(clicked()), SIGNAL(sARAP()));


	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(pbShortestPath);
	layout->addWidget(pbMinTree);
	layout->addWidget(pbMeanCurvature);
	layout->addWidget(pbGaussianCurvature);
	layout->addWidget(pBilateralNormalFilter);
	layout->addWidget(pAddNoise);
	layout->addWidget(pParameterization);
	layout->addWidget(pARAP);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}
