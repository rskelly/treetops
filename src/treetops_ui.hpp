#ifndef __TREETOPS_UI_HPP__
#define __TREETOPS_UI_HPP__

#include <QWidget>
#include <QMessageBox>
#include <QtCore>
#include <QDir>

#include "util.hpp"
#include "treetops.hpp"
#include "ui_treetops.h"

using namespace geo::treetops::config;

namespace geo {

	namespace treetops {

		class TreetopsCallbacks: public QObject, public geo::util::Callbacks {
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

		class TTWorkerThread;

		class TreetopsForm: public QObject, public Ui::TreetopsForm {
			friend class TTWorkerThread;
			Q_OBJECT
		private:
			bool m_cancel;
			QWidget *m_form;
			geo::util::Callbacks *m_callbacks;
			TTWorkerThread *m_workerThread;
			QDir m_last;
			TreetopsConfig m_config;

			// Check if the program is runnable; set buttons accordingly.
			void checkRun();

			// Update the view.
			void updateView();

			// Reset the progress bars and status message.
			void resetProgress();

		public:
			TreetopsForm();
			void setupUi(QWidget *parent);
			void show();
			virtual ~TreetopsForm();

		public slots:
			void doSmoothChanged(bool);
			void doTopsChanged(bool);
			void doCrownsChanged(bool);

			void smoothWindowSizeChanged(int);
			void smoothSigmaChanged(double);
			void smoothOriginalCHMChanged(QString);
			void smoothSmoothedCHMChanged(QString);
			void smoothSmoothedCHMDriverChanged(QString);

			void topsThresholdsChanged(QString);
			void topsSmoothedCHMChanged(QString);
			void topsTreetopsDatabaseChanged(QString);
			void topsTreetopsDatabaseDriverChanged(QString);
			void topsTreetopsSRIDClicked();
			void topsTreetopsSRIDChanged(int);
			void topsMaxNullsChanged(double);

			void crownsThresholdsChanged(QString);
			void crownsTreetopsDatabaseChanged(QString);
			void crownsSmoothedCHMChanged(QString);
			void crownsOriginalCHMChanged(QString);
			void crownsCrownsRasterChanged(QString);
			void crownsCrownsRasterDriverChanged(QString);
			void crownsCrownsDatabaseChanged(QString);
			void crownsCrownsDatabaseDriverChanged(QString);
			void crownsDoDatabaseChanged(bool);
			void crownsUpdateHeightsChanged(bool);
			void crownsRemoveHolesChanged(bool);
			void crownsRemoveDanglesChanged(bool);

			void exitClicked();
			void runClicked();
			void cancelClicked();
			void helpClicked();

			void smoothOriginalCHMClicked();
			void smoothSmoothedCHMClicked();

			void topsThresholdsClicked();
			void topsSmoothedCHMClicked();
			void topsTreetopsDatabaseClicked();

			void crownsCrownsDatabaseClicked();
			void crownsCrownsRasterClicked();
			void crownsTreetopsDatabaseClicked();
			void crownsSmoothedCHMClicked();
			void crownsOriginalCHMClicked();
			void crownsThresholdsClicked();

			void done();
		};

		class TTWorkerThread: public QThread {
		private:
			TreetopsForm *m_parent;
			std::string m_message;
			bool m_isError;
			void run();

		public:
			void init(TreetopsForm *parent);
			virtual ~TTWorkerThread();
			std::string message() const;
			bool isError() const;
		};

	}

}

#endif

