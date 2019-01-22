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

#include "geo.hpp"
#include "raster.hpp"
#include "treetops.hpp"
#include "treetops_ui.hpp"
#include "ui_util.hpp"
#include "settings.hpp"

using namespace geo::ui;
using namespace geo::ui::util;
using namespace geo::treetops;
using namespace geo::treetops::config;

std::string __lastDir(Settings& settings, const std::string& filename) {
	if (Util::isFile(filename)) {
		settings.lastDir = Util::parent(filename);
	}
	else {
		settings.lastDir = filename;
	}
	return filename;
}

// TreetopsCallbacks implementation

void TreetopsCallbacks::stepCallback(float status) const {
	emit stepProgress((int) std::round(status * 100));
}

void TreetopsCallbacks::overallCallback(float status) const {
	emit overallProgress((int) std::round(status * 100));
}

void TreetopsCallbacks::statusCallback(const std::string &msg) const {
	emit statusUpdate(qstr(msg));
}


void TTClockThread::init(TreetopsForm *parent) {
	m_parent = parent;
}

void TTClockThread::start() {
	m_running = true;
	QThread::start();
}

void TTClockThread::stop() {
	m_running = false;
}

void TTClockThread::run() {
	m_sw.reset();
	m_sw.start();
	while(m_running) {
		m_parent->setRunTime(m_sw.time());
		sleep(1);
	}
}

TTClockThread::~TTClockThread() {
}

// WorkerThread implementation

void TTWorkerThread::run() {

	using namespace geo::treetops;
	using namespace geo::treetops::config;

	// Clear the message to indicate no issues.
	m_message.clear();
	m_toFix = FixNone;

	Treetops t;

	try {
		// Reset error state.
		reset();

		// Set the pointer to the callback object.
		t.setCallbacks(m_parent->m_callbacks);

		const TreetopsConfig &config = m_parent->m_config;
		const TreetopsCallbacks *cb = (TreetopsCallbacks *) m_parent->m_callbacks;

		// Calculate the number of steps to complete the job(s).
		int steps = (((int)config.doSmoothing()) + ((int)config.doTops()) + ((int)config.doCrowns())) * 2;
		int step = 0;

		if (cb)
			cb->overallCallback(0.01f);

		if (config.doSmoothing()) {
			cb->overallCallback((float) ++step / steps);
			t.smooth(config, m_parent->m_cancel);
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doTops()) {
			cb->overallCallback((float) ++step / steps);
			try {
				t.treetops(config, m_parent->m_cancel);
			} catch(const geo::treetops::util::DBConvertException& ex) {
				m_toFix |= FixTops;
				m_message = "Saving to the selected database format has failed. The output has been converted to SQLite.";
			}
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doCrowns()) {
			cb->overallCallback((float) ++step / steps);
			try {
				t.treecrowns(config, m_parent->m_cancel);
			} catch(const geo::treetops::util::DBConvertException& ex) {
				m_toFix |= FixCrowns;
				m_message = "Saving to the selected database format has failed. The output has been converted to SQLite.";
			}
			cb->overallCallback((float) ++step / steps);
		}

		cb->statusCallback("Done.");
		cb->overallCallback(1.0f);

	} catch (const std::exception& e) {
		m_message = stripBoost(e.what());
		m_isError = true;
		try {
			t.cleanup();
		} catch(std::exception& ex) {
			g_warn(ex.what());
		}
	}
}

void TTWorkerThread::init(TreetopsForm *parent) {
	m_parent = parent;
	reset();
}

void TTWorkerThread::reset() {
	m_isError = false;
	m_message.clear();
}

std::string TTWorkerThread::message() const {
	return m_message;
}

bool TTWorkerThread::isError() const {
	return m_isError;
}

int TTWorkerThread::toFix() const {
	return m_toFix;
}

TTWorkerThread::~TTWorkerThread(){}


// TreetopsForm implementation

TreetopsForm::TreetopsForm() :
	Ui::TreetopsForm(),
	m_cancel(false),
	m_form(nullptr),
	m_callbacks(nullptr),
	m_workerThread(nullptr),
	m_clockThread(nullptr) {
}

void TreetopsForm::topsConvertFix() {
	std::string text = Util::extension(m_config.treetopsDatabase());
	std::string tops = m_config.treetopsDatabase().substr(0, m_config.treetopsDatabase().find(text)) + ".sqlite";
	m_config.setActive(false);
	m_config.setTreetopsDatabase(tops);
	m_config.setTreetopsDatabaseDriver("SQLite");
	if(m_config.crownsDatabaseDriver() != "SQLite") {
		std::string cext = Util::extension(m_config.crownsDatabase());
		std::string crowns = m_config.crownsDatabase().substr(0, m_config.crownsDatabase().find(cext)) + ".sqlite";
		m_config.setCrownsDatabase(crowns);
		m_config.setCrownsDatabaseDriver("SQLite");
	}
	m_config.setActive(true);
	m_config.update(TreetopsDatabase|TreetopsDatabaseDriver|CrownsDatabase|CrownsDatabaseDriver);
}

void TreetopsForm::crownsConvertFix() {
	std::string ext = Util::extension(m_config.crownsDatabase());
	std::string crowns = m_config.crownsDatabase().substr(0, m_config.crownsDatabase().find(ext)) + ".sqlite";
	m_config.setActive(false);
	m_config.setCrownsDatabase(crowns);
	m_config.setCrownsDatabaseDriver("SQLite");
	m_config.setActive(true);
	m_config.update(CrownsDatabase|CrownsDatabaseDriver);
}

void TreetopsForm::setRunTime(const std::string& time) {
	lblRunTime->setText(QString(time.c_str()));
}

TreetopsForm::~TreetopsForm() {
	// Save the settings.
	m_settings.save(m_config);
	delete m_form;
	delete m_callbacks;
	if(m_clockThread) {
		m_clockThread->stop();
		m_clockThread->wait();
		delete m_clockThread;
	}
	if (m_workerThread) {
		m_workerThread->wait();
		delete m_workerThread;
	}
}

void TreetopsForm::showForm() {
	if(!m_form) {
		m_form = new QWidget();
		this->setupUi(m_form);
	}
	m_form->show();
}

void TreetopsForm::loadSettings() {

	// Populate fields with saved or default values.
	// -- smoothing
	grpSmoothing->setChecked(m_config.doSmoothing());
	spnSmoothWindow->setValue(m_config.smoothWindowSize());
	spnSmoothSigma->setValue(m_config.smoothSigma());

	txtSettingsFile->setText(qstr(m_config.settings()));

	txtOriginalCHM->setText(qstr(m_config.originalCHM()));
	spnOriginalCHMBand->setValue(m_config.originalCHMBand());
	txtSmoothedCHM->setText(qstr(m_config.smoothedCHM()));
	if(!m_config.smoothedCHMDriver().empty())
		cboSmoothedCHMDriver->setCurrentText(qstr(m_config.smoothedCHMDriver()));

	txtTreetopsDatabase->setText(qstr(m_config.treetopsDatabase()));
	if(!m_config.treetopsDatabaseDriver().empty())
		cboTreetopsDatabaseDriver->setCurrentText(qstr(m_config.treetopsDatabaseDriver()));

	txtCrownsRaster->setText(qstr(m_config.crownsRaster()));
	cboCrownsRasterDriver->setCurrentText(qstr(m_config.crownsRasterDriver()));
	txtCrownsDatabase->setText(qstr(m_config.crownsDatabase()));
	if(!m_config.crownsDatabaseDriver().empty())
		cboCrownsDatabaseDriver->setCurrentText(qstr(m_config.crownsDatabaseDriver()));

	// -- tops
	grpTops->setChecked(m_config.doTops());
	txtTopsThresholds->setText(qstr(m_config.topsThresholdsList()));
	spnTopsMaxNulls->setValue(m_config.topsMaxNulls());

	// -- crowns
	grpCrowns->setChecked(m_config.doCrowns());
	txtCrownsThresholds->setText(qstr(m_config.crownsThresholdsList()));
	chkCrownsDoDatabase->setChecked(m_config.crownsDoDatabase());
	chkCrownsUpdateHeights->setChecked(m_config.crownsUpdateHeights());
	chkCrownsRemoveHoles->setEnabled(m_config.crownsDoDatabase() && m_config.doCrowns());
	chkCrownsRemoveDangles->setEnabled(m_config.crownsDoDatabase() && m_config.doCrowns());
	chkCrownsRemoveHoles->setChecked(m_config.crownsRemoveHoles());
	chkCrownsRemoveDangles->setChecked(m_config.crownsRemoveDangles());

}

void TreetopsForm::setupUi(QWidget *form) {
	Ui::TreetopsForm::setupUi(form);
	m_form = form;

	QString title = form->windowTitle();
	form->setWindowTitle(title + " <Rev: " + stringyx(GIT_REV) + ">");

	// Create callbacks and worker thread
	m_callbacks = new TreetopsCallbacks();
	m_workerThread = new TTWorkerThread();
	m_workerThread->init(this);
	m_clockThread = new TTClockThread();
	m_clockThread->init(this);

	// Populate combos.
	QStringList rasterDrivers;
	for(const auto &it : geo::raster::Raster::drivers({"GTiff", "HFA"}))
		rasterDrivers << qstr(it.first);

	QStringList vectorDrivers;
	for(const auto &it : geo::db::DB::drivers({"ESRI Shapefile", "SQLite"}))
		vectorDrivers << qstr(it.first);

	cboSmoothedCHMDriver->addItems(rasterDrivers);
	cboTreetopsDatabaseDriver->addItems(vectorDrivers);
	cboCrownsRasterDriver->addItems(rasterDrivers);
	cboCrownsDatabaseDriver->addItems(vectorDrivers);

	// Connect events
	connect(btnSettingsFile, SIGNAL(clicked()), this, SLOT(settingsFileClicked()));
	connect(txtSettingsFile, SIGNAL(textChanged(QString)), this, SLOT(settingsFileChanged(QString)));

	connect(txtOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(originalCHMChanged(QString)));
	connect(spnOriginalCHMBand, SIGNAL(valueChanged(int)), this, SLOT(originalCHMBandChanged(int)));
	connect(txtSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothedCHMChanged(QString)));
	connect(btnOriginalCHM, SIGNAL(clicked()), this, SLOT(originalCHMClicked()));
	connect(btnSmoothedCHM, SIGNAL(clicked()), this, SLOT(smoothedCHMClicked()));
	connect(cboSmoothedCHMDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(smoothedCHMDriverChanged(QString)));

	connect(txtTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(treetopsDatabaseChanged(QString)));
	connect(btnTreetopsDatabase, SIGNAL(clicked()), this, SLOT(treetopsDatabaseClicked()));
	connect(cboTreetopsDatabaseDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(treetopsDatabaseDriverChanged(QString)));

	connect(txtCrownsRaster, SIGNAL(textChanged(QString)), this, SLOT(crownsRasterChanged(QString)));
	connect(txtCrownsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsDatabaseChanged(QString)));
	connect(btnCrownsRaster, SIGNAL(clicked()), this, SLOT(crownsRasterClicked()));
	connect(cboCrownsRasterDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(crownsRasterDriverChanged(QString)));
	connect(btnCrownsDatabase, SIGNAL(clicked()), this, SLOT(crownsDatabaseClicked()));
	connect(cboCrownsDatabaseDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(crownsDatabaseDriverChanged(QString)));

	// -- section toggles
	connect(grpSmoothing, SIGNAL(toggled(bool)), this, SLOT(doSmoothChanged(bool)));
	connect(grpTops, SIGNAL(toggled(bool)), this, SLOT(doTopsChanged(bool)));
	connect(grpCrowns, SIGNAL(toggled(bool)), this, SLOT(doCrownsChanged(bool)));
	// -- smoothing
	connect(spnSmoothWindow, SIGNAL(valueChanged(int)), this, SLOT(smoothWindowSizeChanged(int)));
	connect(spnSmoothSigma, SIGNAL(valueChanged(double)), this, SLOT(smoothSigmaChanged(double)));

	// -- tops
	// TODO: Needs validator, see #113. connect(txtTopsThresholds, SIGNAL(textEdited(QString)), this, SLOT(topsThresholdsChanged(QString)));
	connect(txtTopsThresholds, SIGNAL(editingFinished()), this, SLOT(topsThresholdsEditingFinished()));
	connect(spnTopsMaxNulls, SIGNAL(valueChanged(double)), this, SLOT(topsMaxNullsChanged(double)));
	connect(btnTopsThresholds, SIGNAL(clicked()), this, SLOT(topsThresholdsClicked()));
	// -- crowns
	// TODO: Needs validator, see #113. connect(txtCrownsThresholds, SIGNAL(textEdited(QString)), this, SLOT(crownsThresholdsChanged(QString)));
	connect(txtCrownsThresholds, SIGNAL(editingFinished()), this, SLOT(crownsThresholdsEditingFinished()));
	connect(chkCrownsDoDatabase, SIGNAL(toggled(bool)), this, SLOT(crownsDoDatabaseChanged(bool)));
	connect(chkCrownsUpdateHeights, SIGNAL(toggled(bool)), this, SLOT(crownsUpdateHeightsChanged(bool)));
	connect(btnCrownsThresholds, SIGNAL(clicked()), this, SLOT(crownsThresholdsClicked()));

	connect(chkCrownsRemoveHoles, SIGNAL(toggled(bool)), this, SLOT(crownsRemoveHolesChanged(bool)));
	connect(chkCrownsRemoveDangles, SIGNAL(toggled(bool)), this, SLOT(crownsRemoveDanglesChanged(bool)));

	// -- program buttons
	connect(btnExit, SIGNAL(clicked()), this, SLOT(exitClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(helpClicked()));

	// -- callbacks
	if (m_callbacks) {
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(statusUpdate(QString)), lblStatus, SLOT(setText(QString)));
	}
	// -- worker thread.
	connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

	m_config.setListener(this);
	m_config.setActive(true);
	m_config.update(TopsThresholds|CrownsThresholds);
	checkRun();
}

void TreetopsForm::resetProgress() {
	prgStep->setValue(0);
	lblStatus->setText("[Not Started]");
}

void TreetopsForm::settingsFileClicked() {
	std::string filename;
	getOutputFile(m_form, "Settings File", m_settings.lastDir, ALL_PATTERN, filename, false);
	txtSettingsFile->setText(QString(filename.c_str()));
}

void TreetopsForm::settingsFileChanged(QString filename) {
	m_config.setSettings(filename.toStdString());
}

void TreetopsForm::topsMaxNullsChanged(double maxNulls) {
	m_config.setTopsMaxNulls(maxNulls);
}

void TreetopsForm::crownsRemoveHolesChanged(bool on) {
	m_config.setCrownsRemoveHoles(on);
}

void TreetopsForm::crownsRemoveDanglesChanged(bool on) {
	m_config.setCrownsRemoveDangles(on);
}

void TreetopsForm::updateView() {
}

void TreetopsForm::originalCHMClicked() {
	std::string filename;
	getInputFile(m_form, "CHM for Smoothing", m_settings.lastDir, ALL_PATTERN, filename);
	bool active = m_config.setActive(false);
	if(m_config.smoothedCHMDriver().empty())
		m_config.setSmoothedCHMDriver(cboSmoothedCHMDriver->currentText().toStdString());
	if(m_config.treetopsDatabaseDriver().empty())
		m_config.setTreetopsDatabaseDriver(cboTreetopsDatabaseDriver->currentText().toStdString());
	if(m_config.crownsRasterDriver().empty())
		m_config.setCrownsRasterDriver(cboCrownsRasterDriver->currentText().toStdString());
	if(m_config.crownsDatabaseDriver().empty())
		m_config.setCrownsDatabaseDriver(cboCrownsDatabaseDriver->currentText().toStdString());
	m_config.setActive(active);
	m_config.setOriginalCHM(filename, true);
}

void TreetopsForm::originalCHMBandChanged(int band) {
	m_config.setOriginalCHMBand(band);
}

void TreetopsForm::smoothedCHMClicked() {
	std::string oldExt = Util::extension(m_config.smoothedCHM());
	std::string filename;
	getOutputFile(m_form, "Smoothed CHM", m_settings.lastDir, ALL_PATTERN, filename);
	m_config.setSmoothedCHM(filename);
}

void TreetopsForm::smoothedCHMDriverChanged(QString text) {
	m_config.setSmoothedCHMDriver(text.toStdString());
}

void TreetopsForm::originalCHMChanged(QString text) {
	m_config.setOriginalCHM(__lastDir(m_settings, text.toStdString()));
}

void TreetopsForm::smoothedCHMChanged(QString text) {
	m_config.setSmoothedCHM(__lastDir(m_settings, text.toStdString()));
}

void TreetopsForm::treetopsDatabaseChanged(QString text) {
	m_config.setTreetopsDatabase(__lastDir(m_settings, text.toStdString()));
}

void TreetopsForm::treetopsDatabaseClicked() {
	std::string oldExt = Util::extension(m_config.treetopsDatabase());
	std::string filename;
	getOutputFile(m_form, "Treetops Database", m_settings.lastDir, ALL_PATTERN, filename);
	m_config.setTreetopsDatabase(filename);
}

void TreetopsForm::treetopsDatabaseDriverChanged(QString text) {
	m_config.setTreetopsDatabaseDriver(text.toStdString());
}

void TreetopsForm::topsThresholdsClicked() {
	std::vector<TopThreshold> thresholds = m_config.topsThresholds();
	getTopsThresholds(m_form, thresholds);
	m_config.setTopsThresholds(thresholds);
}

void TreetopsForm::crownsThresholdsClicked() {
	std::vector<CrownThreshold> thresholds = m_config.crownsThresholds();
	getCrownsThresholds(m_form, thresholds);
	m_config.setCrownsThresholds(thresholds);
}

void TreetopsForm::crownsRasterClicked() {
	std::string oldExt = Util::extension(m_config.crownsRaster());
	std::string filename;
	getOutputFile(m_form, "Crowns Raster", m_settings.lastDir, ALL_PATTERN, filename);
	m_config.setCrownsRaster(filename);
}

void TreetopsForm::crownsRasterDriverChanged(QString text) {
	m_config.setCrownsRasterDriver(text.toStdString());
}

void TreetopsForm::crownsDatabaseClicked() {
	std::string oldExt = Util::extension(m_config.crownsDatabase());
	std::string filename;
	getOutputFile(m_form, "Crowns Database", m_settings.lastDir, ALL_PATTERN, filename);
	m_config.setCrownsDatabase(filename);
}

void TreetopsForm::crownsDatabaseDriverChanged(QString text) {
	m_config.setCrownsDatabaseDriver(text.toStdString());
}

void TreetopsForm::crownsDoDatabaseChanged(bool state) {
	m_config.setCrownsDoDatabase(state);
}

void TreetopsForm::crownsUpdateHeightsChanged(bool state) {
	m_config.setCrownsUpdateHeights(state);
}

void TreetopsForm::doSmoothChanged(bool v) {
	m_config.setDoSmoothing(v);
}

void TreetopsForm::doTopsChanged(bool v) {
	m_config.setDoTops(v);
	if(!m_config.doTops())
		grpCrowns->setChecked(false);
}

void TreetopsForm::doCrownsChanged(bool doCrowns) {
	if(doCrowns && !m_config.doTops())
		grpTops->setChecked(true);
	m_config.setDoCrowns(doCrowns);
}

void TreetopsForm::crownsRasterChanged(QString text) {
	std::string oldExt = Util::extension(m_config.crownsRaster());
	m_config.setCrownsRaster(__lastDir(m_settings, text.toStdString()));
	if(oldExt != Util::extension(m_config.crownsRaster()))
		cboCrownsRasterDriver->setCurrentText("");
}

void TreetopsForm::crownsDatabaseChanged(QString text) {
	std::string oldExt = Util::extension(m_config.crownsDatabase());
	m_config.setCrownsDatabase(__lastDir(m_settings, text.toStdString()));
	if(oldExt != Util::extension(m_config.crownsDatabase()))
		cboCrownsDatabaseDriver->setCurrentText("");
}

void TreetopsForm::topsThresholdsChanged(QString thresh) {
	m_config.parseTopsThresholds(thresh.toStdString());
}

// TODO: Temporary, see #113.
void TreetopsForm::topsThresholdsEditingFinished() {
	QString thresh = txtTopsThresholds->text();
	m_config.parseTopsThresholds(thresh.toStdString());
}

void TreetopsForm::crownsThresholdsChanged(QString thresh) {
	m_config.parseCrownsThresholds(thresh.toStdString());
}

// TODO: Temporary, see #113.
void TreetopsForm::crownsThresholdsEditingFinished() {
	QString thresh = txtCrownsThresholds->text();
	m_config.parseCrownsThresholds(thresh.toStdString());
}

void TreetopsForm::smoothWindowSizeChanged(int size) {
	m_config.setSmoothWindowSize(size);
}

void TreetopsForm::smoothSigmaChanged(double sigma) {
	m_config.setSmoothSigma(sigma);
}

void TreetopsForm::runClicked() {
	if (m_workerThread->isRunning())
		return;
	m_cancel = false;
	btnRun->setEnabled(false);
	btnCancel->setEnabled(true);
	btnExit->setEnabled(false);
	m_clockThread->start();
	m_workerThread->start();
	resetProgress();
	checkRun();
}

void TreetopsForm::done() {
	m_clockThread->stop();
	m_clockThread->wait();
	if(m_workerThread->toFix() & FixTops)
		topsConvertFix();
	if(m_workerThread->toFix() & FixCrowns)
		crownsConvertFix();
	if (m_workerThread->isError()) {
		errorDialog(m_form, "Error", m_workerThread->message());
		resetProgress();
	} else if(!m_workerThread->message().empty()) {
		infoDialog(m_form, "Notice", m_workerThread->message());
	}
	checkRun();
}

void TreetopsForm::exitClicked() {
	g_debug("quit");
	m_settings.save(m_config);
	m_form->close();
}

void TreetopsForm::cancelClicked() {
	g_debug("cancel");
	m_cancel = true;
	checkRun();
}

void TreetopsForm::helpClicked() {
	g_debug("help");
	QDesktopServices::openUrl(QUrl("https://github.com/rskelly/treetops/wiki/Tree-Tops-and-Crowns", QUrl::TolerantMode));
}

void TreetopsForm::checkRun() {
	if(m_workerThread) {
		btnRun->setEnabled(m_config.canRun() && !m_workerThread->isRunning());
		btnCancel->setEnabled(m_workerThread->isRunning());
		btnExit->setEnabled(!m_workerThread->isRunning());
	} else {
		btnRun->setEnabled(false);
		btnCancel->setEnabled(false);
		btnExit->setEnabled(true);
	}
}

void TreetopsForm::configUpdate(TreetopsConfig& config, long field) {

	if(field & SettingsFile) {
		std::string filename = m_config.settings();
		if(Util::exists(filename)) {
			QMessageBox::StandardButton reply = QMessageBox::question(this, "Settings",
					"A settings file exists with this name. Click 'Open' to use the saved settings or 'Reset' to overwrite the file with the new settings.",
					QMessageBox::Open|QMessageBox::Reset);
			if(reply != QMessageBox::Reset) {
				m_config.setActive(false);
				m_settings.load(m_config, filename);
				loadSettings();
				m_config.setActive(true);
				return;
			}
		}
	}

	config.lock();

	if(field & SettingsFile)
		txtSettingsFile->setText(qstr(m_config.settings()));

	if(field & DoSmoothing)
		grpSmoothing->layout()->setEnabled(m_config.doSmoothing());

	if(field & DoTops)
		grpTops->layout()->setEnabled(m_config.doTops());

	if(field & DoCrowns)
		grpCrowns->layout()->setEnabled(m_config.doCrowns());

	if(field & OriginalCHM)
		txtOriginalCHM->setText(qstr(m_config.originalCHM()));

	if(field & SmoothedCHM)
		txtSmoothedCHM->setText(qstr(m_config.smoothedCHM()));

	if(field & TopsThresholds)
		txtTopsThresholds->setText(qstr(m_config.topsThresholdsList()));

	if(field & TreetopsDatabase)
		txtTreetopsDatabase->setText(qstr(m_config.treetopsDatabase()));

	if(field & TreetopsDatabaseDriver)
		cboTreetopsDatabaseDriver->setCurrentText(qstr(m_config.treetopsDatabaseDriver()));

	if((field & DoCrowns) || (field & CrownsDoDatabase)) {
		bool doCrownsAndDb = m_config.doCrowns() && m_config.crownsDoDatabase();
		chkCrownsRemoveHoles->setEnabled(doCrownsAndDb);
		chkCrownsRemoveDangles->setEnabled(doCrownsAndDb);
	}

	if(field & CrownsDatabase)
		txtCrownsDatabase->setText(qstr(m_config.crownsDatabase()));

	if(field & CrownsDatabaseDriver)
		cboCrownsDatabaseDriver->setCurrentText(qstr(m_config.crownsDatabaseDriver()));

	if(field & CrownsRaster)
		txtCrownsRaster->setText(qstr(m_config.crownsRaster()));

	if(field & CrownsThresholds)
		txtCrownsThresholds->setText(qstr(m_config.crownsThresholdsList()));

	m_settings.save(config);

	config.unlock();

	checkRun();
}
