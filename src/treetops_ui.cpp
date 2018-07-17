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
#include "crs_selector_ui.hpp"
#include "ui_util.hpp"
#include "settings.hpp"

using namespace geo::ui;
using namespace geo::ui::util;
using namespace geo::treetops;
using namespace geo::treetops::config;

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

	Treetops t;
	try {
		// Reset error state.
		reset();

		// Set the pointer to the callback object.
		t.setCallbacks(m_parent->m_callbacks);

		const TreetopsConfig &config = m_parent->m_config;
		const TreetopsCallbacks *cb = (TreetopsCallbacks *) m_parent->m_callbacks;

		// Calculate the number of steps to complete the job(s).
		int steps = (((int)config.doSmoothing) + ((int)config.doTops) + ((int)config.doCrowns)) * 2;
		int step = 0;

		if (cb)
			cb->overallCallback(0.01f);

		if (config.doSmoothing) {
			cb->overallCallback((float) ++step / steps);
			t.smooth(config, m_parent->m_cancel);
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doTops) {
			cb->overallCallback((float) ++step / steps);
			t.treetops(config, m_parent->m_cancel);
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doCrowns) {
			cb->overallCallback((float) ++step / steps);
			t.treecrowns(config, m_parent->m_cancel);
			cb->overallCallback((float) ++step / steps);
		}

		cb->statusCallback("Done.");
		cb->overallCallback(1.0f);

	} catch (const std::exception &e) {
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

void TreetopsForm::setRunTime(const std::string& time) {
	lblRunTime->setText(QString(time.c_str()));
}

TreetopsForm::~TreetopsForm() {
	// Save the settings.
	m_settings.save(m_config);

	delete m_form;
	delete m_callbacks;
	if(m_clockThread) {
		m_clockThread->exit(0);
		delete m_clockThread;
	}
	if (m_workerThread) {
		m_workerThread->exit(0);
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

	QStringList rasterDrivers;
	rasterDrivers << "";
	for(const auto &it : geo::raster::Raster::drivers())
		rasterDrivers << qstr(it.first);

	QStringList vectorDrivers;
	vectorDrivers << "";
	for(const auto &it : geo::db::DB::drivers())
		vectorDrivers << qstr(it.first);

	// Populate fields with saved or default values.
	// -- smoothing
	grpSmoothing->setChecked(m_config.doSmoothing);
	spnSmoothWindow->setValue(m_config.smoothWindowSize);
	spnSmoothSigma->setValue(m_config.smoothSigma);
	txtSmoothOriginalCHM->setText(qstr(m_config.smoothOriginalCHM));
	txtSmoothSmoothedCHM->setText(qstr(m_config.smoothSmoothedCHM));
	cboSmoothSmoothedCHMDriver->addItems(rasterDrivers);
	cboSmoothSmoothedCHMDriver->setCurrentText(qstr(m_config.smoothSmoothedCHMDriver));

	// -- tops
	grpTops->setChecked(m_config.doTops);
	txtTopsThresholds->setText(qstr(m_config.topsThresholdsList()));
	txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	txtTopsTreetopsDatabase->setText(qstr(m_config.topsTreetopsDatabase));
	cboTopsTreetopsDatabaseDriver->addItems(vectorDrivers);
	cboTopsTreetopsDatabaseDriver->setCurrentText(qstr(m_config.topsTreetopsDatabaseDriver));
	spnTopsTreetopsSRID->setValue(m_config.srid);
	spnTopsMaxNulls->setValue(m_config.topsMaxNulls);

	// -- crowns
	grpCrowns->setChecked(m_config.doCrowns);
	txtCrownsThresholds->setText(qstr(m_config.crownsThresholdsList()));
	txtCrownsOriginalCHM->setText(qstr(m_config.crownsOriginalCHM));
	txtCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights & m_config.doCrowns);
	btnCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights & m_config.doCrowns);
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	txtCrownsCrownsRaster->setText(qstr(m_config.crownsCrownsRaster));
	cboCrownsCrownsRasterDriver->addItems(rasterDrivers);
	cboCrownsCrownsRasterDriver->setCurrentText(qstr(m_config.crownsCrownsRasterDriver));
	txtCrownsCrownsDatabase->setText(qstr(m_config.crownsCrownsDatabase));
	cboCrownsCrownsDatabaseDriver->addItems(vectorDrivers);
	cboCrownsCrownsDatabaseDriver->setCurrentText(qstr(m_config.crownsCrownsDatabaseDriver));
	chkCrownsDoDatabase->setChecked(m_config.crownsDoDatabase);
	chkCrownsUpdateHeights->setChecked(m_config.crownsUpdateHeights);
	txtCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	cboCrownsCrownsDatabaseDriver->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	btnCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveHoles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveDangles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveHoles->setChecked(m_config.crownsRemoveHoles);
	chkCrownsRemoveDangles->setChecked(m_config.crownsRemoveDangles);
}

void TreetopsForm::setupUi(QWidget *form) {
	Ui::TreetopsForm::setupUi(form);
	m_form = form;

	if(!m_settings.load(m_config))
		g_warn("No settings file remembered. Starting fresh.")

	loadSettings();

	// Create callbacks and worker thread
	m_callbacks = new TreetopsCallbacks();
	m_workerThread = new TTWorkerThread();
	m_workerThread->init(this);
	m_clockThread = new TTClockThread();
	m_clockThread->init(this);

	// Connect events
	connect(btnSettingsFile, SIGNAL(clicked()), this, SLOT(settingsFileClicked()));
	connect(txtSettingsFile, SIGNAL(textChanged(QString)), this, SLOT(settingsFileChanged(QString)));
	// -- section toggles
	connect(grpSmoothing, SIGNAL(toggled(bool)), this, SLOT(doSmoothChanged(bool)));
	connect(grpTops, SIGNAL(toggled(bool)), this, SLOT(doTopsChanged(bool)));
	connect(grpCrowns, SIGNAL(toggled(bool)), this, SLOT(doCrownsChanged(bool)));
	// -- smoothing
	connect(spnSmoothWindow, SIGNAL(valueChanged(int)), this, SLOT(smoothWindowSizeChanged(int)));
	connect(spnSmoothSigma, SIGNAL(valueChanged(double)), this, SLOT(smoothSigmaChanged(double)));
	connect(txtSmoothOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothOriginalCHMChanged(QString)));
	connect(txtSmoothSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothSmoothedCHMChanged(QString)));
	connect(btnSmoothOriginalCHM, SIGNAL(clicked()), this, SLOT(smoothOriginalCHMClicked()));
	connect(btnSmoothSmoothedCHM, SIGNAL(clicked()), this, SLOT(smoothSmoothedCHMClicked()));
	connect(cboSmoothSmoothedCHMDriver, SIGNAL(currentTextChanged(QString)),
			this, SLOT(smoothSmoothedCHMDriverChanged(QString)));
	// -- tops
	connect(txtTopsThresholds, SIGNAL(textEdited(QString)), this, SLOT(topsThresholdsChanged(QString)));
	connect(spnTopsTreetopsSRID, SIGNAL(valueChanged(int)), this, SLOT(topsTreetopsSRIDChanged(int)));
	connect(spnTopsMaxNulls, SIGNAL(valueChanged(double)), this, SLOT(topsMaxNullsChanged(double)));
	connect(txtTopsSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(topsSmoothedCHMChanged(QString)));
	connect(txtTopsTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(topsTreetopsDatabaseChanged(QString)));
	connect(btnTopsSmoothedCHM, SIGNAL(clicked()), this, SLOT(topsSmoothedCHMClicked()));
	connect(btnTopsTreetopsDatabase, SIGNAL(clicked()), this, SLOT(topsTreetopsDatabaseClicked()));
	connect(cboTopsTreetopsDatabaseDriver, SIGNAL(currentTextChanged(QString)),
			this, SLOT(topsTreetopsDatabaseDriverChanged(QString)));
	connect(btnTopsThresholds, SIGNAL(clicked()), this, SLOT(topsThresholdsClicked()));
	// -- crowns
	connect(txtCrownsThresholds, SIGNAL(textEdited(QString)), this, SLOT(crownsThresholdsChanged(QString)));
	connect(txtCrownsOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(crownsOriginalCHMChanged(QString)));
	connect(txtCrownsSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(crownsSmoothedCHMChanged(QString)));
	connect(txtCrownsTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsTreetopsDatabaseChanged(QString)));
	connect(txtCrownsCrownsRaster, SIGNAL(textChanged(QString)), this, SLOT(crownsCrownsRasterChanged(QString)));
	connect(txtCrownsCrownsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsCrownsDatabaseChanged(QString)));
	connect(btnCrownsOriginalCHM, SIGNAL(clicked()), this, SLOT(crownsOriginalCHMClicked()));
	connect(btnCrownsSmoothedCHM, SIGNAL(clicked()), this, SLOT(crownsSmoothedCHMClicked()));
	connect(btnCrownsTreetopsDatabase, SIGNAL(clicked()), this, SLOT(crownsTreetopsDatabaseClicked()));
	connect(btnCrownsCrownsRaster, SIGNAL(clicked()), this, SLOT(crownsCrownsRasterClicked()));
	connect(cboCrownsCrownsRasterDriver, SIGNAL(currentTextChanged(QString)),
			this, SLOT(crownsCrownsRasterDriverChanged(QString)));
	connect(btnCrownsCrownsDatabase, SIGNAL(clicked()), this, SLOT(crownsCrownsDatabaseClicked()));
	connect(cboCrownsCrownsDatabaseDriver, SIGNAL(currentTextChanged(QString)),
			this, SLOT(crownsCrownsDatabaseDriverChanged(QString)));
	connect(btnTopsTreetopsSRID, SIGNAL(clicked()), this, SLOT(topsTreetopsSRIDClicked()));
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

	checkRun();
}

void TreetopsForm::resetProgress() {
	prgStep->setValue(0);
	//prgOverall->setValue(0);
	lblStatus->setText("[Not Started]");
}

void TreetopsForm::settingsFileClicked() {
	std::string filename;
	getOutputFile(m_form, "Settings File", m_settings.outputLastDir, ALL_PATTERN, filename, false);
	txtSettingsFile->setText(QString(filename.c_str()));
}

void TreetopsForm::settingsFileChanged(QString value) {
	std::string filename = txtSettingsFile->text().toStdString();
	if(Util::exists(filename)) {
		QMessageBox::StandardButton reply = QMessageBox::question(this, "Settings",
				"If you choose an extant file, existing settings will be replaced. Do you want to continue?",
				QMessageBox::Yes|QMessageBox::No);
		if(reply == QMessageBox::No)
			return;
	}
	m_settings.load(m_config, filename);
	loadSettings();
}

void TreetopsForm::topsTreetopsSRIDChanged(int srid) {
	m_config.srid = srid;
	checkRun();
}

void TreetopsForm::topsMaxNullsChanged(double maxNulls) {
	m_config.topsMaxNulls = maxNulls;
	checkRun();
}

void TreetopsForm::crownsRemoveHolesChanged(bool on) {
	m_config.crownsRemoveHoles = on;
	checkRun();
}

void TreetopsForm::crownsRemoveDanglesChanged(bool on) {
	m_config.crownsRemoveDangles = on;
	checkRun();
}

void TreetopsForm::topsTreetopsSRIDClicked() {
	CRSSelector cs(m_form);
	cs.enableVertical(false);
	cs.setHorizontalSRID(m_config.srid);
	if (cs.exec())
		spnTopsTreetopsSRID->setValue(cs.getHorizontalSRID());
}

void TreetopsForm::updateView() {
}

void TreetopsForm::smoothOriginalCHMClicked() {
	getInputFile(m_form, "CHM for Smoothing", m_settings.originalCHMLastDir, ALL_PATTERN, m_config.smoothOriginalCHM);
	txtSmoothOriginalCHM->setText(qstr(m_config.smoothOriginalCHM));
	m_config.crownsOriginalCHM = m_config.smoothOriginalCHM;
	txtCrownsOriginalCHM->setText(qstr(m_config.crownsOriginalCHM));
	checkRun();
}

void TreetopsForm::smoothSmoothedCHMClicked() {
	std::string oldExt = Util::extension(m_config.smoothSmoothedCHM);
	getOutputFile(m_form, "Smoothed CHM", m_settings.outputLastDir, ALL_PATTERN, m_config.smoothSmoothedCHM);
	txtSmoothSmoothedCHM->setText(qstr(m_config.smoothSmoothedCHM));
	m_config.topsSmoothedCHM = m_config.smoothSmoothedCHM;
	txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	m_config.crownsSmoothedCHM = m_config.smoothSmoothedCHM;
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	if(oldExt != Util::extension(m_config.smoothSmoothedCHM))
		cboSmoothSmoothedCHMDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::smoothSmoothedCHMDriverChanged(QString text) {
	m_config.smoothSmoothedCHMDriver = text.toStdString();
	checkRun();
}

void TreetopsForm::topsSmoothedCHMClicked() {
	getInputFile(m_form, "Smoothed CHM for Treetops", m_settings.smoothedCHMLastDir, ALL_PATTERN, m_config.topsSmoothedCHM);
	txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	m_config.crownsSmoothedCHM = m_config.topsSmoothedCHM;
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	checkRun();
}

void TreetopsForm::topsTreetopsDatabaseClicked() {
	std::string oldExt = Util::extension(m_config.topsTreetopsDatabase);
	getOutputFile(m_form, "Treetops Database", m_settings.topsDatabaseLastDir, ALL_PATTERN, m_config.topsTreetopsDatabase);
	txtTopsTreetopsDatabase->setText(qstr(m_config.topsTreetopsDatabase));
	m_config.crownsTreetopsDatabase = m_config.topsTreetopsDatabase;
	txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	if(oldExt != Util::extension(m_config.topsTreetopsDatabase))
		cboTopsTreetopsDatabaseDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::topsTreetopsDatabaseDriverChanged(QString text) {
	m_config.topsTreetopsDatabaseDriver = text.toStdString();
	checkRun();
}

void TreetopsForm::topsThresholdsClicked() {
	getTopsThresholds(m_form, m_config.topsThresholds);
	txtTopsThresholds->setText(qstr(m_config.topsThresholdsList()));
	checkRun();
}

void TreetopsForm::crownsThresholdsClicked() {
	getCrownsThresholds(m_form, m_config.crownsThresholds);
	txtCrownsThresholds->setText(qstr(m_config.crownsThresholdsList()));
	checkRun();
}

void TreetopsForm::crownsOriginalCHMClicked() {
	getInputFile(m_form, "Original CHM for Crown Delineation", m_settings.originalCHMLastDir, ALL_PATTERN, m_config.crownsOriginalCHM);
	txtCrownsOriginalCHM->setText(qstr(m_config.crownsOriginalCHM));
	checkRun();
}

void TreetopsForm::crownsSmoothedCHMClicked() {
	getInputFile(m_form, "Smoothed CHM for Crown Delineation", m_settings.smoothedCHMLastDir, ALL_PATTERN, m_config.crownsSmoothedCHM);
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseClicked() {
	getInputFile(m_form, "Treetops Database", m_settings.topsDatabaseLastDir, ALL_PATTERN, m_config.crownsTreetopsDatabase);
	txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	checkRun();
}

void TreetopsForm::crownsCrownsRasterClicked() {
	std::string oldExt = Util::extension(m_config.crownsCrownsRaster);
	getOutputFile(m_form, "Crowns Raster", m_settings.crownsRasterLastDir, ALL_PATTERN, m_config.crownsCrownsRaster);
	txtCrownsCrownsRaster->setText(qstr(m_config.crownsCrownsRaster));
	if(oldExt != Util::extension(m_config.crownsCrownsRaster))
		cboCrownsCrownsRasterDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::crownsCrownsRasterDriverChanged(QString text) {
	m_config.crownsCrownsRasterDriver = text.toStdString();
	checkRun();
}


void TreetopsForm::crownsCrownsDatabaseClicked() {
	std::string oldExt = Util::extension(m_config.crownsCrownsDatabase);
	getOutputFile(m_form, "Crowns Database", m_settings.crownsDatabaseLastDir, ALL_PATTERN, m_config.crownsCrownsDatabase);
	txtCrownsCrownsDatabase->setText(qstr(m_config.crownsCrownsDatabase));
	if(oldExt != Util::extension(m_config.crownsCrownsDatabase))
		cboCrownsCrownsDatabaseDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::crownsCrownsDatabaseDriverChanged(QString text) {
	m_config.crownsCrownsDatabaseDriver = text.toStdString();
	checkRun();
}

void TreetopsForm::crownsDoDatabaseChanged(bool state) {
	m_config.crownsDoDatabase = state;
	txtCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	cboCrownsCrownsDatabaseDriver->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	btnCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveHoles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveDangles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	checkRun();
}

void TreetopsForm::crownsUpdateHeightsChanged(bool state) {
	m_config.crownsUpdateHeights = state;
	txtCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights);
	btnCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights);
	checkRun();
}

void TreetopsForm::doSmoothChanged(bool v) {
	m_config.doSmoothing = v;
	grpSmoothing->layout()->setEnabled(v);
	checkRun();
}

void TreetopsForm::doTopsChanged(bool v) {
	m_config.doTops = v;
	grpTops->layout()->setEnabled(v);
	checkRun();
}

void TreetopsForm::doCrownsChanged(bool v) {
	m_config.doCrowns = v;
	grpCrowns->layout()->setEnabled(v);
	txtCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights && m_config.doCrowns);
	btnCrownsOriginalCHM->setEnabled(m_config.crownsUpdateHeights && m_config.doCrowns);
	txtCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	cboCrownsCrownsDatabaseDriver->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	btnCrownsCrownsDatabase->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveHoles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	chkCrownsRemoveDangles->setEnabled(m_config.crownsDoDatabase && m_config.doCrowns);
	checkRun();
}

void TreetopsForm::crownsOriginalCHMChanged(QString file) {
	m_config.crownsOriginalCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsSmoothedCHMChanged(QString file) {
	m_config.crownsSmoothedCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseChanged(QString file) {
	m_config.crownsTreetopsDatabase = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsCrownsRasterChanged(QString file) {
	std::string oldExt = Util::extension(m_config.crownsCrownsRaster);
	m_config.crownsCrownsRaster = file.toStdString();
	if(oldExt != Util::extension(m_config.crownsCrownsRaster))
		cboCrownsCrownsRasterDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::crownsCrownsDatabaseChanged(QString file) {
	std::string oldExt = Util::extension(m_config.crownsCrownsDatabase);
	m_config.crownsCrownsDatabase = file.toStdString();
	if(oldExt != Util::extension(m_config.crownsCrownsDatabase))
		cboCrownsCrownsDatabaseDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::topsThresholdsChanged(QString thresh) {
	m_config.parseTopsThresholds(thresh.toStdString());
	checkRun();
}

void TreetopsForm::crownsThresholdsChanged(QString thresh) {
	std::cerr << thresh.toStdString() << "\n";
	m_config.parseCrownsThresholds(thresh.toStdString());
	checkRun();
}

void TreetopsForm::topsSmoothedCHMChanged(QString file) {
	m_config.topsSmoothedCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::topsTreetopsDatabaseChanged(QString file) {
	std::string oldExt = Util::extension(m_config.topsTreetopsDatabase);
	m_config.topsTreetopsDatabase = file.toStdString();
	if(oldExt != Util::extension(m_config.topsTreetopsDatabase))
		cboTopsTreetopsDatabaseDriver->setCurrentText("");
	checkRun();
}

void TreetopsForm::smoothWindowSizeChanged(int size) {
	m_config.smoothWindowSize = size;
	checkRun();
}

void TreetopsForm::smoothSigmaChanged(double sigma) {
	m_config.smoothSigma = sigma;
	checkRun();
}

void TreetopsForm::smoothOriginalCHMChanged(QString file) {
	m_config.smoothOriginalCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::smoothSmoothedCHMChanged(QString file) {
	std::string oldExt = Util::extension(m_config.smoothSmoothedCHM);
	m_config.smoothSmoothedCHM = file.toStdString();
	if(oldExt != Util::extension(m_config.smoothSmoothedCHM))
		cboSmoothSmoothedCHMDriver->setCurrentText("");
	checkRun();
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
	if (m_workerThread->isError()) {
		errorDialog(m_form, "Error", m_workerThread->message());
		resetProgress();
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
	btnRun->setEnabled(m_config.canRun() && !m_workerThread->isRunning());
	btnCancel->setEnabled(m_workerThread->isRunning());
	btnExit->setEnabled(!m_workerThread->isRunning());
	m_settings.save(m_config);
}
