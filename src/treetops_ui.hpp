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

		class TreetopsForm: public QDialog, public Ui::TreetopsForm, public TreetopsConfigListener {
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

			// Handle updates from the config.
			void configUpdate(TreetopsConfig& config, long field);

		public:
			TreetopsForm();
			void setupUi(QWidget *parent);
			void showForm();
			void setRunTime(const std::string& time);
			void loadSettings();
			virtual ~TreetopsForm();

			/**
			 * Called when Treetops fails to convert the internal database to the chosen format.
			 * Will cause the app to change the filename and driver to SQLite and move the file into
			 * the appropriate place. Modifies the tops and crowns DB.
			 */
			void topsConvertFix();

			/**
			 * Called when Treetops fails to convert the internal database to the chosen format.
			 * Will cause the app to change the filename and driver to SQLite and move the file into
			 * the appropriate place. Modifies the crowns DB.
			 */
			void crownsConvertFix();

		public slots:
			void settingsFileClicked();
			void settingsFileChanged(QString);

			void doSmoothChanged(bool);
			void doTopsChanged(bool);
			void doCrownsChanged(bool);

			void originalCHMChanged(QString);
			void originalCHMBandChanged(int);
			void smoothedCHMChanged(QString);
			void smoothedCHMDriverChanged(QString);
			void treetopsDatabaseChanged(QString);
			void treetopsDatabaseDriverChanged(QString);
			void crownsRasterChanged(QString);
			void crownsRasterDriverChanged(QString);
			void crownsDatabaseChanged(QString);
			void crownsDatabaseDriverChanged(QString);

			void crownsDatabaseClicked();
			void crownsRasterClicked();
			void treetopsDatabaseClicked();
			void smoothedCHMClicked();
			void originalCHMClicked();

			void smoothWindowSizeChanged(int);
			void smoothSigmaChanged(double);

			void topsThresholdsChanged(QString);
			void topsThresholdsClicked();
			void topsTreetopsSRIDClicked();
			void topsTreetopsSRIDChanged(int);
			void topsMaxNullsChanged(double);

			void crownsThresholdsChanged(QString);
			void crownsDoDatabaseChanged(bool);
			void crownsUpdateHeightsChanged(bool);
			void crownsRemoveHolesChanged(bool);
			void crownsRemoveDanglesChanged(bool);

			void exitClicked();
			void runClicked();
			void cancelClicked();
			void helpClicked();

			void crownsThresholdsClicked();

			void done();
		};

		/**
		 * A thread for performing the work of the application.
		 */
		class TTWorkerThread: public QThread {
		private:
			TreetopsForm *m_parent;		///<! The form that owns this thread.
			std::string m_message;		///<! The error message from the last error.
			bool m_isError;				///<! True if the thread is running or exiting in an error state.

			/**
			 * Run the thread.
			 */
			void run();

			/**
			 * Reset the error message, etc.
			 */
			void reset();

		public:

			/**
			 * Initialize the thread with a pointer to its parent form.
			 *
			 * @param parent The parent widget.
			 */
			void init(TreetopsForm *parent);

			/**
			 * Destroy the thread.
			 */
			virtual ~TTWorkerThread();

			/**
			 * Return the last error message.
			 *
			 * @return The last error message.
			 */
			std::string message() const;

			/**
			 * Return true if an error has occurred.
			 *
			 * @return True if an error has occurred.
			 */
			bool isError() const;
		};

		/**
		 * A thread for handling the UI clock.
		 */
		class TTClockThread: public QThread {
		private:
			TreetopsForm *m_parent;		///<! The parent widget.
			geo::util::Stopwatch m_sw;	///<! A stopwatch instance.
			bool m_running;				///<! True when the application is running.

			/**
			 * Run the clock thread.
			 */
			void run();

		public:
			/**
			 * Initialize with a pointer to the parent.
			 *
			 * @param parent A pointer to the parent widget.
			 */
			void init(TreetopsForm *parent);

			/**
			 * Destroy the clock thread.
			 */
			virtual ~TTClockThread();

			/**
			 * Return the time.
			 *
			 * @return The time.
			 */
			std::string time() const;

			/**
			 * Start the clock.
			 */
			void start();

			/**
			 * Stop the clock.
			 */
			void stop();
		};

	}

}

#endif

