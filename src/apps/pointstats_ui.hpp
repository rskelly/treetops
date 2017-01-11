#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <cmath>

#include <QWidget>
#include <QDir>
#include <QMessageBox>
#include <QtCore>

#include "util.hpp"
#include "pointstats.hpp"
#include "ui_pointstats.h"
#include "filelist.hpp"

namespace geotools {

namespace ui {

class PointStatsCallbacks: public QObject, public geotools::util::Callbacks {
	Q_OBJECT
public:
	void stepCallback(float status) const;
	void overallCallback(float status) const;
	void statusCallback(const std::string &msg) const;
signals:
	void stepProgress(int) const;
	void overallProgress(int) const;
	void statusUpdate(QString) const;
};

class PSWorkerThread;

class PointStatsForm: public QWidget, public Ui::PointStatsForm {
	friend class PSWorkerThread;
	Q_OBJECT
private:
	QWidget *m_form;
	QDir m_last;
	PSWorkerThread *m_workerThread;
	geotools::util::Callbacks *m_callbacks;
	geotools::point::PointStatsConfig m_config;
	FileList m_fileList;
	QString m_filter;
	bool m_cancel;

	void updateFileList();
	void updateFileButtons();
	void checkRun();
	void updateSnapUi();
	void updateAreaUi();
	void updateTypeUi();

public:
	PointStatsForm(QWidget *p = Q_NULLPTR);
	void setupUi(QWidget *Form);
	~PointStatsForm();
	void show();

public slots:
	void destFileClicked();
	void destFileChanged(QString);
	void snapModeChanged(int);
	void cancelClicked();
	void runClicked();
	void exitClicked();
	void crsConfigClicked();
	void typeSelected(int);
	void threadsChanged(int);
	void quantileChanged(int);
	void quantilesChanged(int);
	void attributeSelected(int);
	void resolutionXChanged(double);
	void resolutionYChanged(double);
	void gapFunctionSelected(int);
	void originXChanged(double);
	void originYChanged(double);
	void quantileFilterFromChanged(int);
	void quantileFilterToChanged(int);
	void quantileFilterChanged(int);
	void maxAngleChanged(int);
	void classItemClicked(QListWidgetItem*);
	void fileListChanged();
	void gapThresholdChanged(double);
	void areaModeChanged(int);
	void areaSizeChanged(double);
	void done();
};

class PSWorkerThread: public QThread {
private:
	geotools::util::Bounds m_bounds;
	PointStatsForm *m_parent;
	std::string m_error;
	void run();
public:
	void init(PointStatsForm *parent, const geotools::util::Bounds &bounds);
	bool hasError();
	std::string getError();
};

}

}

#endif

