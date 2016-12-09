#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <vector>
#include <string>
#include <math.h>

#include <QWidget>
#include <QDir>
#include <QtCore>
#include <QMessageBox>

#include "mosaic.hpp"
#include "util.hpp"
#include "ui_mosaic.h"

namespace geotools {

	namespace raster {

		class MosaicCallbacks: public QObject, public geotools::util::Callbacks {
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

	}

	namespace ui {

		class WorkerThread;

		class MosaicForm: public QWidget, public Ui::MosaicForm {
			Q_OBJECT
			friend class WorkerThread;
		private:
			QWidget *m_form;
			int m_tileSize;
			double m_distance;
			int m_threads;
			geotools::raster::MosaicCallbacks *m_callbacks;
			WorkerThread *m_workerThread;
			std::string m_destFile;
			std::vector<std::string> m_tifFiles;
			QDir m_last;

			void updateFileList();
			void updateFileButtons();
			void checkRun();

		public:
			MosaicForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			virtual ~MosaicForm();

		public slots:
			void fileListSelectionChanged();
			void selectFilesClicked();
			void removeFilesClicked();
			void clearFilesClicked();
			void cancelClicked();
			void runClicked();
			void destFileClicked();
			void upClicked();
			void downClicked();
			void distanceChanged(QString);
			void threadsChanged(int);
			void tileSizeChanged(QString);
			void done();
		};

		class WorkerThread: public QThread {
		private:
			MosaicForm *m_parent;

		public:

			void init(MosaicForm *parent) {
				m_parent = parent;
			}

			void run() {
				geotools::raster::Mosaic m;
				m.setCallbacks(m_parent->m_callbacks);
				try {
					m.mosaic(m_parent->m_tifFiles, m_parent->m_destFile,
							m_parent->m_distance, m_parent->m_tileSize,
							m_parent->m_threads);
				} catch (const std::exception &e) {
					QMessageBox err((QWidget *) m_parent);
					err.setText("Error");
					err.setInformativeText(QString(e.what()));
					err.exec();
				}
			}

			virtual ~WorkerThread() {}

		};

	}

}

#endif

