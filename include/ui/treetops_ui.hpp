#ifndef __TREETOPS_UI_HPP__
#define __TREETOPS_UI_HPP__

#include <QtWidgets/QWidget>
#include <QtWidgets/QMessageBox>
#include <QtCore/QDir>
#include <QtCore/QThread>

#include "treetops.hpp"
#include "settings.hpp"

#include "ui/ui_treetops.h"

using namespace tt;
using namespace tt::config;

namespace tt {
namespace ui {

class TreetopsForm: public QDialog, Ui::TreetopsForm {
	Q_OBJECT
private:

	Settings m_settings;

	// Check if the program is runnable; set buttons accordingly.
	void checkRun();

	// Update the view.
	void updateView();

	// Reset the progress bars and status message.
	void resetProgress();

	// Handle updates from the config.
	void configUpdate(long); //TreetopsConfig& config, long field);

public:
	//TreetopsForm(QWidget* parent = nullptr);
	void setupUi(QWidget* parent);
	void showForm();
	void setRunTime(const std::string& time);
	void loadSettings();
	~TreetopsForm() = default;

signals:
	void configUpdateReceived(long);
	void updateReceived(long);
	
public slots:
	void settingsFileClicked();
	void settingsFileChanged(QString);

	void doSmoothChanged(bool);
	void doTopsChanged(bool);
	void doCrownsChanged(bool);

	void originalCHMChanged(QString);
	void originalCHMBandChanged(int);
	void originalCHMClicked();

	void smoothedCHMChanged(QString);
	void smoothedCHMDriverChanged(QString);
	void smoothedCHMClicked();
	void smoothWindowSizeChanged(int);
	void smoothSigmaChanged(double);

	void treetopsDatabaseChanged(QString);
	void treetopsDatabaseDriverChanged(QString);
	void treetopsDatabaseClicked();
	void topsThresholdsChanged(QString);
	void topsThresholdsEditingFinished(); // TODO: Temporary see #113.
	void topsThresholdsClicked();
	void topsMaxNullsChanged(double);

	void crownsRasterChanged(QString);
	void crownsRasterDriverChanged(QString);
	void crownsDatabaseChanged(QString);
	void crownsDatabaseDriverChanged(QString);
	void crownsDatabaseClicked();
	void crownsRasterClicked();
	void crownsThresholdsChanged(QString);
	void crownsThresholdsEditingFinished(); // TODO: Temporary see #113.
	void crownsDoDatabaseChanged(bool);
	void crownsUpdateHeightsChanged(bool);
	void crownsRemoveHolesChanged(bool);
	void crownsRemoveDanglesChanged(bool);
	void crownsThresholdsClicked();
	void crownsKeepSmoothedChanged(bool);

	void exitClicked();
	void runClicked();
	void cancelClicked();
	void helpClicked();

	void handleConfigUpdate(long);

	/**
	 * \brief Called when the worker thread stops.
	 */
	void stopped();

	/**
	 * \brief Called when the worker thread starts.
	 */
	void started();
};

} // ui
} // tt

#endif

