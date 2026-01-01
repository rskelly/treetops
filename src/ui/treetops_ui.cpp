#ifndef GIT_REV
#define GIT_REV 0000000
#endif
#define stringyx(x) stringy(x)
#define stringy(GIT_REV) #GIT_REV


#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <QtCore/QDir>
#include <QtCore/QUrl>
#include <QtCore/QString>
#include <QtGui/QDesktopServices>

#include "constants.hpp"
#include "treetops.hpp"
#include "settings.hpp"

#include "ui/ui_treetops.h"
#include "ui/ui_util.hpp"
#include "ui/treetops_ui.hpp"

using namespace tt;
using namespace tt::ui;
using namespace tt::ui::util;
using namespace tt::config;


// TreetopsForm implementation

void TreetopsForm::setRunTime(const std::string& time) {
	lblRunTime->setText(QString(time.c_str()));
}

/*
TreetopsForm::~TreetopsForm() {
}
*/

void TreetopsForm::showForm() {
	setupUi(this);
	loadSettings();
	show();
}

void TreetopsForm::loadSettings() {

	txtSettingsFile->setText(qstr(m_settings.settingsFile()));

	txtOriginalCHM->setText(qstr(m_settings.get("originalCHM", "")));
	spnOriginalCHMBand->setValue(m_settings.get("originalCHMBand", 0));

	// -- smoothing
	grpSmoothing->setChecked(m_settings.get("doSmoothing", false));
	spnSmoothWindow->setValue(m_settings.get("smoothWindowSize", 0));
	spnSmoothSigma->setValue(m_settings.get("smoothSigma", 0.0));
	txtSmoothedCHM->setText(qstr(m_settings.get("smoothedCHM", "")));
	cboSmoothedCHMDriver->setCurrentText(qstr(m_settings.get("smoothedCHMDriver", RASTER_DRIVERS[0])));

	// -- tops
	grpTops->setChecked(m_settings.get("doTops", false));
	txtTopsThresholds->setText(qstr(m_settings.get("topsThresholdsList", "")));
	spnTopsMaxNulls->setValue(m_settings.get("topsMaxNulls", 0));
	txtTreetopsDatabase->setText(qstr(m_settings.get("treetopsDatabase", "")));
	cboTreetopsDatabaseDriver->setCurrentText(qstr(m_settings.get("treetopsDatabaseDriver", VECTOR_DRIVERS[0])));

	// -- crowns
	grpCrowns->setChecked(m_settings.get("doCrowns", false));
	txtCrownsThresholds->setText(qstr(m_settings.get("crownsThresholdsList", "")));
	chkCrownsDoDatabase->setChecked(m_settings.get("crownsDoDatabase", false));
	chkCrownsUpdateHeights->setChecked(m_settings.get("crownsUpdateHeights", false));
	chkCrownsRemoveHoles->setEnabled(m_settings.get("doCrowns", false) && m_settings.get("crownsDoDatabase", false));
	chkCrownsRemoveDangles->setEnabled(m_settings.get("crownsDoDatabase", false) && m_settings.get("doCrowns", false));
	chkCrownsRemoveHoles->setChecked(m_settings.get("crownsRemoveHoles", false));
	chkCrownsRemoveDangles->setChecked(m_settings.get("crownsRemoveDangles", false));
	chkCrownsKeepSmoothed->setChecked(m_settings.get("crownsKeepSmoothed", false));
	txtCrownsRaster->setText(qstr(m_settings.get("crownsRaster", "")));
	cboCrownsRasterDriver->setCurrentText(qstr(m_settings.get("crownsRasterDriver", "GTiff")));
	txtCrownsDatabase->setText(qstr(m_settings.get("crownsDatabase", "")));
	cboCrownsDatabaseDriver->setCurrentText(qstr(m_settings.get("crownsDatabaseDriver", "Spatialite")));
}

void TreetopsForm::setupUi(QWidget *form) {
	Ui::TreetopsForm::setupUi(form);

	QString title = form->windowTitle();
	form->setWindowTitle(title + " <Rev: " + stringyx(GIT_REV) + ">");

	// Create callbacks and worker thread
	/*
	m_settings.get("setMonitor(new tt::TreetopsMonitor"));
	m_workerThread = new TTWorkerThread();
	m_workerThread->init(this, &m_settings.;
	m_clockThread = new TTClockThread();
	m_clockThread->init(this);
	*/

	// Populate combos.
	QStringList rasterDrivers;
	for(int i = 0; i < sizeof(RASTER_DRIVERS) / sizeof(char*); ++i)
		rasterDrivers << RASTER_DRIVERS[i];

	QStringList vectorDrivers;
	for(int i = 0; i < sizeof(VECTOR_DRIVERS) / sizeof(char*); ++i)
		vectorDrivers << "ESRI Shapefile" << "SQLite"; //qstr(it.first);

	cboSmoothedCHMDriver->addItems(rasterDrivers);
	cboTreetopsDatabaseDriver->addItems(vectorDrivers);
	cboCrownsRasterDriver->addItems(rasterDrivers);
	cboCrownsDatabaseDriver->addItems(vectorDrivers);

	// Connect events
	connect(btnSettingsFile, SIGNAL(clicked()), this, SLOT(settingsFileClicked()));
	connect(txtSettingsFile, SIGNAL(textChanged(QString)), this, SLOT(settingsFileChanged(QString)));

	connect(txtOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(originalCHMChanged(QString)));
	connect(spnOriginalCHMBand, SIGNAL(valueChanged(int)), this, SLOT(originalCHMBandChanged(int)));
	connect(btnOriginalCHM, SIGNAL(clicked()), this, SLOT(originalCHMClicked()));

	// -- smoothing
	connect(spnSmoothWindow, SIGNAL(valueChanged(int)), this, SLOT(smoothWindowSizeChanged(int)));
	connect(spnSmoothSigma, SIGNAL(valueChanged(double)), this, SLOT(smoothSigmaChanged(double)));
	connect(txtSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothedCHMChanged(QString)));
	connect(cboSmoothedCHMDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(smoothedCHMDriverChanged(QString)));
	connect(btnSmoothedCHM, SIGNAL(clicked()), this, SLOT(smoothedCHMClicked()));

	// -- tops
	// TODO: Needs validator, see #113. connect(txtTopsThresholds, SIGNAL(textEdited(QString)), this, SLOT(topsThresholdsChanged(QString)));
	connect(txtTopsThresholds, SIGNAL(editingFinished()), this, SLOT(topsThresholdsEditingFinished()));
	connect(spnTopsMaxNulls, SIGNAL(valueChanged(double)), this, SLOT(topsMaxNullsChanged(double)));
	connect(btnTopsThresholds, SIGNAL(clicked()), this, SLOT(topsThresholdsClicked()));

	connect(txtTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(treetopsDatabaseChanged(QString)));
	connect(cboTreetopsDatabaseDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(treetopsDatabaseDriverChanged(QString)));
	connect(btnTreetopsDatabase, SIGNAL(clicked()), this, SLOT(treetopsDatabaseClicked()));

	// -- crowns
	// TODO: Needs validator, see #113. connect(txtCrownsThresholds, SIGNAL(textEdited(QString)), this, SLOT(crownsThresholdsChanged(QString)));
	connect(txtCrownsThresholds, SIGNAL(editingFinished()), this, SLOT(crownsThresholdsEditingFinished()));
	connect(chkCrownsDoDatabase, SIGNAL(toggled(bool)), this, SLOT(crownsDoDatabaseChanged(bool)));
	connect(chkCrownsUpdateHeights, SIGNAL(toggled(bool)), this, SLOT(crownsUpdateHeightsChanged(bool)));
	connect(btnCrownsThresholds, SIGNAL(clicked()), this, SLOT(crownsThresholdsClicked()));
	connect(chkCrownsRemoveHoles, SIGNAL(toggled(bool)), this, SLOT(crownsRemoveHolesChanged(bool)));
	connect(chkCrownsRemoveDangles, SIGNAL(toggled(bool)), this, SLOT(crownsRemoveDanglesChanged(bool)));

	connect(txtCrownsRaster, SIGNAL(textChanged(QString)), this, SLOT(crownsRasterChanged(QString)));
	connect(txtCrownsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsDatabaseChanged(QString)));
	connect(cboCrownsRasterDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(crownsRasterDriverChanged(QString)));
	connect(cboCrownsDatabaseDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(crownsDatabaseDriverChanged(QString)));
	connect(btnCrownsDatabase, SIGNAL(clicked()), this, SLOT(crownsDatabaseClicked()));
	connect(btnCrownsRaster, SIGNAL(clicked()), this, SLOT(crownsRasterClicked()));
	connect(chkCrownsKeepSmoothed, SIGNAL(toggled(bool)), this, SLOT(crownsKeepSmoothedChanged(bool)));

	// -- section toggles
	connect(grpSmoothing, SIGNAL(toggled(bool)), this, SLOT(doSmoothChanged(bool)));
	connect(grpTops, SIGNAL(toggled(bool)), this, SLOT(doTopsChanged(bool)));
	connect(grpCrowns, SIGNAL(toggled(bool)), this, SLOT(doCrownsChanged(bool)));

	// -- program buttons
	connect(btnExit, SIGNAL(clicked()), this, SLOT(exitClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(helpClicked()));

	// -- handle m_settings.updates in the ui thread.
	connect(this, SIGNAL(updateReceived(long)), this, SLOT(handleConfigUpdate(long)));

	// -- callbacks
	//connect(dynamic_cast<TreetopsMonitor*>(m_settings.monitor()), SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
	//connect(dynamic_cast<TreetopsMonitor*>(m_settings.monitor()), SIGNAL(statusUpdate(QString)), lblStatus, SLOT(setText(QString)));

	// -- worker thread.
	//connect(m_workerThread, SIGNAL(finished()), this, SLOT(stopped()));
	//connect(m_workerThread, SIGNAL(started()), this, SLOT(started()));

	//m_settings.setListener(this);
	//m_settings.setActive(true);
	//m_settings.update(TopsThresholds|CrownsThresholds);
	checkRun();
}

void TreetopsForm::resetProgress() {
	prgStep->setValue(0);
	lblStatus->setText("[Not Started]");
}

void TreetopsForm::settingsFileClicked() {
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getOutputFile(this, "Settings File", lastDir, ALL_PATTERN, path, false);
	m_settings.lastDir(path);
	m_settings.settingsFile(path);
}

void TreetopsForm::settingsFileChanged(QString text) {
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.settingsFile(path);
}

void TreetopsForm::topsMaxNullsChanged(double maxNulls) {
	m_settings.set("topsMaxNulls", maxNulls);
}

void TreetopsForm::crownsRemoveHolesChanged(bool on) {
	m_settings.set("crownsRemoveHoles", on);
}

void TreetopsForm::crownsRemoveDanglesChanged(bool on) {
	m_settings.set("crownsRemoveDangles", on);
}

void TreetopsForm::crownsKeepSmoothedChanged(bool on) {
	m_settings.set("crownsKeepSmoothed", on);
}

void TreetopsForm::updateView() {
	bool enable = true; //!(m_workerThread && m_workerThread->isRunning());
	grpFiles->setEnabled(enable);
	grpSmoothing->setEnabled(enable);
	grpTops->setEnabled(enable);
	grpCrowns->setEnabled(enable);
}

void TreetopsForm::originalCHMClicked() {
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getInputFile(this, "CHM for Smoothing", lastDir, ALL_PATTERN, path);
	bool active = m_settings.get("active", false);
	m_settings.lastDir(path);
	m_settings.set("smoothedCHMDriver", sstr(cboSmoothedCHMDriver->currentText()));
	m_settings.set("treetopsDatabaseDriver", sstr(cboTreetopsDatabaseDriver->currentText()));
	m_settings.set("crownsRasterDriver", sstr(cboCrownsRasterDriver->currentText()));
	m_settings.set("crownsDatabaseDriver", sstr(cboCrownsDatabaseDriver->currentText()));
	m_settings.set("active", active);
	m_settings.set("originalCHM", path);
	txtOriginalCHM->setText(QString(path.c_str()));
}

void TreetopsForm::originalCHMBandChanged(int band) {
	m_settings.set("originalCHMBand", band);
}

void TreetopsForm::smoothedCHMClicked() {
	std::string oldExt = tt::util::extension(m_settings.get("smoothedCHM", ""));
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getOutputFile(this, "Smoothed CHM", lastDir, ALL_PATTERN, path);
	m_settings.lastDir(path);
	m_settings.set("smoothedCHM", path);
	txtSmoothedCHM->setText(QString(path.c_str()));
}

void TreetopsForm::smoothedCHMDriverChanged(QString text) {
	m_settings.set("smoothedCHMDriver", sstr(text));
}

void TreetopsForm::originalCHMChanged(QString text) {
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.set("originalCHM", path);
}

void TreetopsForm::smoothedCHMChanged(QString text) {
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.set("smoothedCHM", path);
}

void TreetopsForm::treetopsDatabaseChanged(QString text) {
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.set("treetopsDatabase", path);
}

void TreetopsForm::treetopsDatabaseClicked() {
	std::string oldExt = tt::util::extension(m_settings.get("treetopsDatabase", ""));
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getOutputFile(this, "Treetops Database", lastDir, ALL_PATTERN, path);
	m_settings.lastDir(path);
	m_settings.set("treetopsDatabase", path);
	txtTreetopsDatabase->setText(QString(path.c_str()));
}

void TreetopsForm::treetopsDatabaseDriverChanged(QString text) {
	m_settings.set("treetopsDatabaseDriver", sstr(text));
}

void TreetopsForm::topsThresholdsClicked() {
	std::vector<TopThreshold> thresholds = m_settings.topThresholds();
	getTopsThresholds(this, thresholds);
	m_settings.topThresholds(thresholds);
}

void TreetopsForm::crownsThresholdsClicked() {
	std::vector<CrownThreshold> thresholds = m_settings.crownThresholds();
	getCrownsThresholds(this, thresholds);
	m_settings.crownThresholds(thresholds);
}

void TreetopsForm::crownsRasterClicked() {
	// Get the extension of the stored file.
	std::string ext = tt::util::extension(m_settings.get("crownsRaster", ""));
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getOutputFile(this, "Crowns Raster", lastDir, ALL_PATTERN, path);
	m_settings.lastDir(path);
	m_settings.set("crownsRaster", path);
	// If the extension doesn't match unselect the driver. TODO: Automate this.
	if(ext != tt::util::extension(path))
		cboCrownsDatabaseDriver->setCurrentText("");
}

void TreetopsForm::crownsDatabaseClicked() {
	// Get the extension of the stored file.
	std::string ext = tt::util::extension(m_settings.get("crownsDatabase", ""));
	std::string path;
	std::string lastDir = m_settings.lastDir();
	getOutputFile(this, "Crowns Database", lastDir, ALL_PATTERN, path);
	m_settings.lastDir(path);
	m_settings.set("crownsDatabase", path);
	// If the extension doesn't match unselect the driver. TODO: Automate this.
	if(ext != tt::util::extension(path))
		cboCrownsDatabaseDriver->setCurrentText("");
	txtCrownsDatabase->setText(QString(path.c_str()));
}

void TreetopsForm::crownsDatabaseDriverChanged(QString text) {
	m_settings.set("crownsDatabaseDriver", sstr(text));
}

void TreetopsForm::crownsDoDatabaseChanged(bool state) {
	m_settings.set("crownsDoDatabase", state);
}

void TreetopsForm::crownsUpdateHeightsChanged(bool state) {
	m_settings.set("crownsUpdateHeights", state);
}

void TreetopsForm::doSmoothChanged(bool v) {
	m_settings.set("doSmoothing", v);
}

void TreetopsForm::doTopsChanged(bool v) {
	m_settings.set("doTops", v);
	if(!m_settings.get("doTops", true))
		grpCrowns->setChecked(false);
}

void TreetopsForm::doCrownsChanged(bool doCrowns) {
	if(doCrowns && !m_settings.get("doTops", true))
		grpTops->setChecked(true);
	m_settings.set("doCrowns", doCrowns);
}

void TreetopsForm::crownsRasterChanged(QString text) {
	// Get the extension of the file in the settings.
	std::string ext = tt::util::extension(m_settings.get("crownsRaster", ""));
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.set("crownsRaster", path);
	// If the extension has changed, unselect the driver. TODO: Automate the selection.
	if(ext != tt::util::extension(path))
		cboCrownsRasterDriver->setCurrentText("");
}

void TreetopsForm::crownsRasterDriverChanged(QString text) {
	m_settings.set("crownsRasterDriver", sstr(text));
}

void TreetopsForm::crownsDatabaseChanged(QString text) {
	// Get the extension of the file in the settings.
	std::string ext = tt::util::extension(m_settings.get("crownsDatabase", ""));
	std::string path = text.toStdString();
	m_settings.lastDir(path);
	m_settings.set("crownsDatabase", path);
	// If the extension has changed, unselect the driver. TODO: Automate the selection.
	if(ext != tt::util::extension(m_settings.get("crownsDatabase", "")))
		cboCrownsDatabaseDriver->setCurrentText("");
}

void TreetopsForm::topsThresholdsChanged(QString thresh) {
	m_settings.parseTopThresholds(sstr(thresh));
}

// TODO: Temporary, see #113.
void TreetopsForm::topsThresholdsEditingFinished() {
	QString thresh = txtTopsThresholds->text();
	m_settings.parseTopThresholds(sstr(thresh));
}

void TreetopsForm::crownsThresholdsChanged(QString thresh) {
	m_settings.parseCrownThresholds(sstr(thresh));
}

// TODO: Temporary, see #113.
void TreetopsForm::crownsThresholdsEditingFinished() {
	QString thresh = txtCrownsThresholds->text();
	m_settings.parseCrownThresholds(sstr(thresh));
}

void TreetopsForm::smoothWindowSizeChanged(int size) {
	m_settings.set("smoothWindowSize", size);
}

void TreetopsForm::smoothSigmaChanged(double sigma) {
	m_settings.set("smoothSigma", sigma);
}

void TreetopsForm::runClicked() {
	//if (m_workerThread->isRunning())
	//	return;
	//m_settings.monitor()->setCanceled(false);
	//m_workerThread->start();
}

void TreetopsForm::started() {
	btnRun->setEnabled(false);
	btnCancel->setEnabled(true);
	btnExit->setEnabled(false);
	//m_clockThread->start();
	checkRun();
	updateView();
}

void TreetopsForm::stopped() {
	//m_clockThread->stop();
	//m_clockThread->wait();
	/*
	if (m_workerThread->isError()) {
		errorDialog(this, "Error", m_workerThread->message());
		resetProgress();
	} else if(!m_workerThread->message().empty()) {
		infoDialog(this, "Notice", m_workerThread->message());
	}
		*/
	checkRun();
	updateView();
}

void TreetopsForm::exitClicked() {
	_debug("quit");
	close();
}

void TreetopsForm::cancelClicked() {
	_debug("cancel");
	//m_settings.monitor()->cancel();
	checkRun();
}

void TreetopsForm::helpClicked() {
	_debug("help");
	QDesktopServices::openUrl(QUrl("https://github.com/rskelly/treetops/wiki/Tree-Tops-and-Crowns", QUrl::TolerantMode));
}

void TreetopsForm::checkRun() {
	/*
	if(m_workerThread) {
		btnRun->setEnabled(m_settings.canRun() && !m_workerThread->isRunning());
		btnCancel->setEnabled(m_workerThread->isRunning());
		btnExit->setEnabled(!m_workerThread->isRunning());
	} else {
		btnRun->setEnabled(false);
		btnCancel->setEnabled(false);
		btnExit->setEnabled(true);
	}
		*/
}

void TreetopsForm::configUpdate(long field) {
	// Emit the event that will trigger handleConfig update in the UI thread.
	emit configUpdateReceived(field);
}

void TreetopsForm::handleConfigUpdate(long field) {
	/*
	if(field & SettingsFile) {
		std::string path = m_settings.settingsFile();
		if(isfile(path)) {
			QMessageBox::StandardButton reply = QMessageBox::question(this, "Settings",
					"A settings file exists with this name. Click 'Open' to use the saved settings or 'Reset' to overwrite the file with the new settings.",
					QMessageBox::Open|QMessageBox::Reset);
			if(reply != QMessageBox::Reset) {
				m_settings.set("active", false);
				m_settings.load(m_settings. path);
				loadSettings();
				m_settings.set("active", true);
				return;
			}
		}
	}


	if(field & SettingsFile)
		txtSettingsFile->setText(qstr(m_settings.settingsFile()));

	if(field & DoSmoothing)
		grpSmoothing->layout()->setEnabled(m_settings.doSmoothing());

	if(field & DoTops)
		grpTops->layout()->setEnabled(m_settings.doTops());

	if(field & DoCrowns)
		grpCrowns->layout()->setEnabled(m_settings.doCrowns());

	if(field & OriginalCHM)
		txtOriginalCHM->setText(qstr(m_settings.originalCHM()));

	if(field & SmoothedCHM)
		txtSmoothedCHM->setText(qstr(m_settings.smoothedCHM()));

	if(field & TopsThresholds)
		txtTopsThresholds->setText(qstr(m_settings.topsThresholdsList()));

	if(field & TreetopsDatabase)
		txtTreetopsDatabase->setText(qstr(m_settings.treetopsDatabase()));

	if(field & TreetopsDatabaseDriver)
		cboTreetopsDatabaseDriver->setCurrentText(qstr(m_settings.treetopsDatabaseDriver()));

	if((field & DoCrowns) || (field & CrownsDoDatabase)) {
		bool doCrownsAndDb = m_settings.doCrowns() && m_settings.crownsDoDatabase();
		chkCrownsRemoveHoles->setEnabled(doCrownsAndDb);
		chkCrownsRemoveDangles->setEnabled(doCrownsAndDb);
	}

	if(field & CrownsDatabase)
		txtCrownsDatabase->setText(qstr(m_settings.crownsDatabase()));

	if(field & CrownsDatabaseDriver)
		cboCrownsDatabaseDriver->setCurrentText(qstr(m_settings.crownsDatabaseDriver()));

	if(field & CrownsRaster)
		txtCrownsRaster->setText(qstr(m_settings.crownsRaster()));

	if(field & CrownsThresholds)
		txtCrownsThresholds->setText(qstr(m_settings.crownThresholdsList()));

	if(field & CrownsKeepSmoothed)
		chkCrownsKeepSmoothed->setChecked(m_settings.crownsKeepSmoothed());

	if(field & TopsDBFormatChanged)
		infoDialog(this, "Settings Changed", "The tops database is too large for the Shapefile format. Changed to SQLite.");

	if(field & CrownsDBFormatChanged)
		infoDialog(this, "Settings Changed", "The crowns database is too large for the Shapefile format. Changed to SQLite.");

	m_settings.save(m_settings.;
	*/

	checkRun();
}
