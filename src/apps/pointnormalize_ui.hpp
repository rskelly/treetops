/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pointnormalize_ui.hpp
 * Author: rob
 *
 * Created on November 27, 2016, 12:59 PM
 */

#ifndef POINTNORMALIZE_UI_HPP
#define POINTNORMALIZE_UI_HPP

#include <QDir>
#include <QWidget>
#include <QThread>
#include <QButtonGroup>

#include "ui_pointnormalize.h"
#include "pointnormalize.hpp"
#include "filelist.hpp"
#include "util.hpp"

using namespace geotools::point;

namespace geotools {

namespace ui {

class PointNormalizeCallbacks: public QObject, public geotools::util::Callbacks {
	Q_OBJECT
public:
	void stepCallback(float status) const;
	void overallCallback(float status) const;signals:
	void stepProgress(int) const;
	void overallProgress(int) const;
};

class WorkerThread;

class PointNormalizeForm: public QWidget, public Ui::PointNormalizeForm {
	Q_OBJECT
private:
	QWidget *m_form;
	FileList m_fileList;
	QDir m_last;
	QString m_filter;
	geotools::util::Callbacks *m_callbacks;
	WorkerThread *m_workerThread;
	PointNormalizeConfig m_config;
	bool m_cancel;
	QButtonGroup *m_bgrpGround;
	QButtonGroup *m_bgrpNegative;

public slots:
	void exitClicked();
	void cancelClicked();
	void helpClicked();
	void runClicked();
	void outputFolderClicked();
	void done();
	void checkRun();
	void fileListChanged();
	void keepNegativeToggled(bool);
	void keepGroundToggled(bool);
	void threadsChanged(int);
	void overwriteChanged(bool);
	void bufferChanged(double);

public:
	PointNormalizeForm(QWidget *p = Q_NULLPTR);
	void setupUi(QWidget *Form);
	const PointNormalizeConfig& config();
	geotools::util::Callbacks* callbacks();
	bool* cancel();
	~PointNormalizeForm();
};

class WorkerThread: public QThread {
private:
	PointNormalizeForm *m_parent;
	std::string m_error;
	void run();
public:
	void init(PointNormalizeForm *parent);
	bool hasError();
	std::string getError();
};

}
}

#endif /* POINTNORMALIZE_UI_HPP */

