#ifndef __TREETOPS_UI_HPP__
#define __TREETOPS_UI_HPP__

#include <QWidget>
#include <QMessageBox>
#include <QtCore>
#include <QDir>

#include "util.hpp"
#include "treetops.hpp"
#include "ui_treetops.h"

using namespace geotools::treetops::config;

namespace geotools {

	namespace treetops {

		class TreetopsCallbacks: public QObject, public geotools::util::Callbacks {
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

		class TreetopsForm: public QWidget, public Ui::TreetopsForm {
			friend class TTWorkerThread;
			Q_OBJECT
		private:
			bool m_cancel;
			QWidget *m_form;
			geotools::util::Callbacks *m_callbacks;
			TTWorkerThread *m_workerThread;
			QDir m_last;
			TreetopsConfig m_config;

			std::vector<QWidget*> m_smoothGroup;
			std::vector<QWidget*> m_topsGroup;
			std::vector<QWidget*> m_crownsGroup;

			// Check if the program is runnable; set buttons accordingly.
			void checkRun();

			// Update the view.
			void updateView();

			// Reset the progress bars and status message.
			void resetProgress();

			// (Dis|En)able a group of components.
			void enableGroup(const std::vector<QWidget*> &grp, bool enable);

		public:
			TreetopsForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *parent);
			~TreetopsForm();

		public slots:
			void doSmoothChanged(bool);
			void doTopsChanged(bool);
			void doCrownsChanged(bool);

			void smoothWindowSizeChanged(int);
			void smoothSigmaChanged(double);
			void smoothOriginalCHMChanged(QString);
			void smoothSmoothedCHMChanged(QString);

			void topsMinHeightChanged(double);
			void topsWindowSizeChanged(int);
			void topsOriginalCHMChanged(QString);
			void topsSmoothedCHMChanged(QString);
			void topsTreetopsDatabaseChanged(QString);
			void topsTreetopsSRIDClicked();
			void topsTreetopsSRIDChanged(int);

			void crownsRadiusChanged(double);
			void crownsHeightFractionChanged(double);
			void crownsMinHeightChanged(double);
			void crownsTreetopsDatabaseChanged(QString);
			void crownsSmoothedCHMChanged(QString);
			void crownsCrownsRasterChanged(QString);
			void crownsCrownsDatabaseChanged(QString);

			void exitClicked();
			void runClicked();
			void cancelClicked();
			void helpClicked();

			void smoothOriginalCHMClicked();
			void smoothSmoothedCHMClicked();

			void topsOriginalCHMClicked();
			void topsSmoothedCHMClicked();
			void topsTreetopsDatabaseClicked();

			void crownsCrownsDatabaseClicked();
			void crownsCrownsRasterClicked();
			void crownsTreetopsDatabaseClicked();
			void crownsSmoothedCHMClicked();

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

