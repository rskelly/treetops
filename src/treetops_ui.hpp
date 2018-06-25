#ifndef __TREETOPS_UI_HPP__
#define __TREETOPS_UI_HPP__

#include <QtWidgets/QWidget>
#include <QtWidgets/QMessageBox>
#include <QtCore/QDir>
#include <QtCore/QThread>

#include "util.hpp"
#include "treetops.hpp"
#include "settings.hpp"
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
		class TTClockThread;

		class TreetopsForm: public QDialog, public Ui::TreetopsForm {
			friend class TTWorkerThread;
			Q_OBJECT
		private:
			bool m_cancel;
			QWidget *m_form;
			geo::util::Callbacks *m_callbacks;
			TTWorkerThread *m_workerThread;
			TTClockThread *m_clockThread;
			TreetopsConfig m_config;

			Settings m_settings;
			std::string m_settingsFile;

			// Check if the program is runnable; set buttons accordingly.
			void checkRun();

			// Update the view.
			void updateView();

			// Reset the progress bars and status message.
			void resetProgress();

		public:
			TreetopsForm();
			void setupUi(QWidget *parent);
			void showForm();
			void setRunTime(const std::string& time);
			void loadSettings();
			virtual ~TreetopsForm();

		public slots:
			void settingsFileClicked();
			void settingsFileChanged(QString);

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

		// A thread for performing the work of the application.
		class TTWorkerThread: public QThread {
		private:
			// The form that owns this thread.
			TreetopsForm *m_parent;
			// The error message from the last error.
			std::string m_message;
			// True if the thread is running or exiting in an error state.
			bool m_isError;

			// Run the thread.
			void run();

			// Reset the error message, etc.
			void reset();

		public:

			// Initialize the thread with a pointer to its parent form.
			void init(TreetopsForm *parent);

			// Destroy the thread.
			virtual ~TTWorkerThread();

			// Return the last error message/
			std::string message() const;

			// Return true if an error has occurred.
			bool isError() const;
		};

		// A thread for handling the UI clock.
		class TTClockThread: public QThread {
		private:
			TreetopsForm *m_parent;
			geo::util::Stopwatch m_sw;
			bool m_running;

			void run();

		public:
			void init(TreetopsForm *parent);
			virtual ~TTClockThread();
			std::string time();
			void start();
			void stop();
		};

	}

}

#endif

