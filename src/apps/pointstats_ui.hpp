#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <cmath>

#include <QtWidgets/QWidget>
#include <QDir>
#include <QMessageBox>
#include <QtCore>

#include "util.hpp"
#include "pointstats.hpp"
#include "ui_pointstats.h"

namespace geotools {

    namespace ui {

        class PointStatsCallbacks : public QObject, public geotools::util::Callbacks {
            Q_OBJECT
        public:
            void stepCallback(float status) const;
            void overallCallback(float status) const;
        signals:
            void stepProgress(int) const;
            void overallProgress(int) const;
        };

        class WorkerThread;

        class PointStatsForm : public QWidget, public Ui::PointStatsForm {
            friend class WorkerThread;
            Q_OBJECT
        private:
            QWidget *m_form;
            QDir m_last;
            WorkerThread *m_workerThread;
            geotools::util::Callbacks *m_callbacks;
            geotools::point::PointStatsConfig m_config;
            
            void updateFileList();
            void updateFileButtons();
            void updateTypeUi();
            void checkRun();

        public:
            PointStatsForm(QWidget *p = Q_NULLPTR);
            void setupUi(QWidget *Form);
            ~PointStatsForm();

        public slots:
            void fileListSelectionChanged();
            void selectFilesClicked();
            void removeFilesClicked();
            void clearFilesClicked();
            void destFileClicked();
            void snapToGridChanged(bool);
            void cancelClicked();
            void runClicked();
            void crsConfigClicked();
            void typeSelected(int);
            void threadsChanged(int);
            void quantileChanged(int);
            void quantilesChanged(int);
            void attributeSelected(int);
            void resolutionChanged(double);
            void gapFunctionSelected(int);
            void originXChanged(double);
            void originYChanged(double);
            void quantileFilterFromChanged(int);
            void quantileFilterToChanged(int);
            void quantileFilterChanged(int);
            void maxAngleChanged(int);
            void classItemClicked(QListWidgetItem*);
            void done();
        };

        class WorkerThread : public QThread {
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

