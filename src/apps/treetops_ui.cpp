#include <QWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QSettings>
#include <QDesktopServices>
#include <QUrl>

#include "geotools.hpp"
#include "treetops.hpp"
#include "treetops_ui.hpp"
#include "crs_selector_ui.hpp"
#include "ui_util.hpp"

using namespace geotools::ui;
using namespace geotools::ui::util;
using namespace geotools::treetops;
using namespace geotools::treetops::config;

QSettings _settings("Treetops", "dijital.ca");

// Load the config settings for this app.
void _loadConfig(TreetopsConfig &config) {
	QSettings qs("TreetopsConfig", "dijital.ca");
	config.srid = qs.value(QString("srid"), config.srid).toInt();
	config.buildIndex =
			qs.value(QString("buildIndex"), config.buildIndex).toBool();
	config.tableCacheSize = qs.value(QString("tableCacheSize"),
			config.tableCacheSize).toInt();
	config.rowCacheSize =
			qs.value(QString("rowCacheSize"), config.rowCacheSize).toInt();
	config.doSmoothing =
			qs.value(QString("doSmoothing"), config.doSmoothing).toBool();
	config.smoothWindowSize = qs.value(QString("smoothingWindowSize"),
			config.smoothWindowSize).toInt();
	config.smoothSigma =
			qs.value(QString("smoothSigma"), config.smoothSigma).toDouble();
	config.smoothOriginalCHM = qs.value(QString("smoothOriginalCHM"),
			qstr(config.smoothOriginalCHM)).toString().toStdString();
	config.smoothSmoothedCHM = qs.value(QString("smoothSmoothedCHM"),
			qstr(config.smoothSmoothedCHM)).toString().toStdString();
	config.doTops = qs.value(QString("doTops"), config.doTops).toBool();
	config.topsMinHeight = qs.value(QString("topsMinHeight"),
			config.topsMinHeight).toDouble();
	config.topsWindowSize = qs.value(QString("topsWindowSize"),
			config.topsWindowSize).toInt();
	config.topsOriginalCHM = qs.value(QString("topsOriginalCHM"),
			qstr(config.topsOriginalCHM)).toString().toStdString();
	config.topsSmoothedCHM = qs.value(QString("topsSmoothedCHM"),
			qstr(config.topsSmoothedCHM)).toString().toStdString();
	config.topsTreetopsDatabase = qs.value(QString("topsTreetopsDatabase"),
			qstr(config.topsTreetopsDatabase)).toString().toStdString();
	config.doCrowns = qs.value(QString("doCrowns"), config.doCrowns).toBool();
	config.crownsRadius =
			qs.value(QString("crownsRadius"), config.crownsRadius).toDouble();
	config.crownsHeightFraction = qs.value(QString("crownsHeightFraction"),
			config.crownsHeightFraction).toDouble();
	config.crownsMinHeight = qs.value(QString("crownsMinHeight"),
			config.crownsMinHeight).toDouble();
	config.crownsSmoothedCHM = qs.value(QString("crownsSmoothedCHM"),
			qstr(config.crownsSmoothedCHM)).toString().toStdString();
	config.crownsTreetopsDatabase = qs.value(QString("crownsTreetopsDatabase"),
			qstr(config.crownsTreetopsDatabase)).toString().toStdString();
	config.crownsCrownsRaster = qs.value(QString("crownsCrownsRaster"),
			qstr(config.crownsCrownsRaster)).toString().toStdString();
	config.crownsCrownsDatabase = qs.value(QString("crownsCrownsDatabase"),
			qstr(config.crownsCrownsDatabase)).toString().toStdString();
}

// Save the config settings for this app.
void _saveConfig(TreetopsConfig &config) {
	QSettings qs("TreetopsConfig", "dijital.ca");
	qs.setValue(QString("srid"), config.srid);
	qs.setValue(QString("buildIndex"), config.buildIndex);
	qs.setValue(QString("tableCacheSize"), config.tableCacheSize);
	qs.setValue(QString("rowCacheSize"), config.rowCacheSize); //(24 * 1024 * 1024),
	qs.setValue(QString("doSmoothing"), config.doSmoothing);
	qs.setValue(QString("smoothingWindowSize"), config.smoothWindowSize);
	qs.setValue(QString("smoothSigma"), config.smoothSigma);
	qs.setValue(QString("smoothOriginalCHM"), qstr(config.smoothOriginalCHM));
	qs.setValue(QString("smoothSmoothedCHM"), qstr(config.smoothSmoothedCHM));
	qs.setValue(QString("doTops"), config.doTops);
	qs.setValue(QString("topsMinHeight"), config.topsMinHeight);
	qs.setValue(QString("topsWindowSize"), config.topsWindowSize);
	qs.setValue(QString("topsOriginalCHM"), qstr(config.topsOriginalCHM));
	qs.setValue(QString("topsSmoothedCHM"), qstr(config.topsSmoothedCHM));
	qs.setValue(QString("topsTreetopsDatabase"),
			qstr(config.topsTreetopsDatabase));
	qs.setValue(QString("doCrowns"), config.doCrowns);
	qs.setValue(QString("crownsRadius"), config.crownsRadius);
	qs.setValue(QString("crownsHeightFraction"), config.crownsHeightFraction);
	qs.setValue(QString("crownsMinHeight"), config.crownsMinHeight);
	qs.setValue(QString("crownsSmoothedCHM"), qstr(config.crownsSmoothedCHM));
	qs.setValue(QString("crownsTreetopsDatabase"),
			qstr(config.crownsTreetopsDatabase));
	qs.setValue(QString("crownsCrownsRaster"), qstr(config.crownsCrownsRaster));
	qs.setValue(QString("crownsCrownsDatabase"),
			qstr(config.crownsCrownsDatabase));

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


// WorkerThread implementation

void TTWorkerThread::run() {

	using namespace geotools::treetops;
	using namespace geotools::treetops::config;

	try {
		Treetops t;
		t.setCallbacks(m_parent->m_callbacks);
		const TreetopsConfig &config = m_parent->m_config;
		const TreetopsCallbacks *cb =
			(TreetopsCallbacks *)m_parent->m_callbacks;

		int steps = (((int)config.doSmoothing) + ((int)config.doTops)
			+ ((int)config.doCrowns)) * 2;
		int step = 0;

		if (cb)
			cb->overallCallback(0.01f);

		if (config.doSmoothing) {
			cb->overallCallback((float) ++step / steps);
			t.smooth(config);
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doTops) {
			cb->overallCallback((float) ++step / steps);
			t.treetops(config);
			cb->overallCallback((float) ++step / steps);
		}

		if (config.doCrowns) {
			cb->overallCallback((float) ++step / steps);
			t.treecrowns(config);
			cb->overallCallback((float) ++step / steps);
		}

		cb->overallCallback(1.0f);

	} catch (const std::exception &e) {
		m_message = stripBoost(e.what());
		m_isError = true;
	}
}

void TTWorkerThread::init(TreetopsForm *parent) {
	m_parent = parent;
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

TreetopsForm::TreetopsForm(QWidget *parent) :
	QWidget(parent),
	m_cancel(false),
	m_form(nullptr),
	m_callbacks(nullptr),
	m_workerThread(nullptr) {
	if(!parent)
		this->setupUi(this);
}

TreetopsForm::~TreetopsForm() {
	_saveConfig(m_config);
	delete m_callbacks;
	if (m_workerThread) {
		m_workerThread->exit(0);
		delete m_workerThread;
	}
}

void TreetopsForm::enableGroup(const std::vector<QWidget*> &grp, bool enable) {
	for(QWidget *w : grp)
		w->setEnabled(enable);
}

void TreetopsForm::setupUi(QWidget *form) {
	Ui::TreetopsForm::setupUi(form);

	m_form = form;

	_loadConfig(m_config);

	// If the settings has a record for the last used directory, use it.
	if (_settings.contains(QString("last_dir"))) {
		m_last.setPath(_settings.value(QString("last_dir")).toString());
		if(!m_last.exists())
			m_last = QDir::home();
	} else {
		m_last = QDir::home();
	}

	// Create callbacks and worker thread
	m_callbacks = new TreetopsCallbacks();
	m_workerThread = new TTWorkerThread();
	m_workerThread->init(this);

	// Create virtual groups
	m_smoothGroup.push_back(spnSmoothWindow);
	m_smoothGroup.push_back(spnSmoothSigma);
	m_smoothGroup.push_back(txtSmoothOriginalCHM);
	m_smoothGroup.push_back(btnSmoothOriginalCHM);
	m_smoothGroup.push_back(txtSmoothSmoothedCHM);
	m_smoothGroup.push_back(btnSmoothSmoothedCHM);
	m_topsGroup.push_back(spnTopsMinHeight);
	m_topsGroup.push_back(spnTopsWindowSize);
	m_topsGroup.push_back(txtTopsOriginalCHM);
	m_topsGroup.push_back(btnTopsOriginalCHM);
	m_topsGroup.push_back(txtTopsSmoothedCHM);
	m_topsGroup.push_back(btnTopsSmoothedCHM);
	m_topsGroup.push_back(txtTopsTreetopsDatabase);
	m_topsGroup.push_back(btnTopsTreetopsDatabase);
	m_topsGroup.push_back(spnTopsTreetopsSRID);
	m_topsGroup.push_back(btnTopsTreetopsSRID);
	m_crownsGroup.push_back(spnCrownsHeightFraction);
	m_crownsGroup.push_back(spnCrownsRadius);
	m_crownsGroup.push_back(spnCrownsMinHeight);
	m_crownsGroup.push_back(txtCrownsSmoothedCHM);
	m_crownsGroup.push_back(btnCrownsSmoothedCHM);
	m_crownsGroup.push_back(txtCrownsTreetopsDatabase);
	m_crownsGroup.push_back(btnCrownsTreetopsDatabase);
	m_crownsGroup.push_back(txtCrownsCrownsRaster);
	m_crownsGroup.push_back(btnCrownsCrownsRaster);
	m_crownsGroup.push_back(txtCrownsCrownsDatabase);
	m_crownsGroup.push_back(btnCrownsCrownsDatabase);

	// Set enabled if so configured.
	enableGroup(m_smoothGroup, m_config.doSmoothing);
	enableGroup(m_topsGroup, m_config.doTops);
	enableGroup(m_crownsGroup, m_config.doCrowns);

	// Populate fields with saved or default values.
	// -- smoothing
	chkEnableSmoothing->setChecked(m_config.doSmoothing);
	spnSmoothWindow->setValue(m_config.smoothWindowSize);
	spnSmoothSigma->setValue(m_config.smoothSigma);
	txtSmoothOriginalCHM->setText(qstr(m_config.smoothOriginalCHM));
	txtSmoothSmoothedCHM->setText(qstr(m_config.smoothSmoothedCHM));
	// -- tops
	chkEnableTops->setChecked(m_config.doTops);
	spnTopsMinHeight->setValue(m_config.topsMinHeight);
	spnTopsWindowSize->setValue(m_config.topsWindowSize);
	txtTopsOriginalCHM->setText(qstr(m_config.topsOriginalCHM));
	txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	txtTopsTreetopsDatabase->setText(qstr(m_config.topsTreetopsDatabase));
	spnTopsTreetopsSRID->setValue(m_config.srid);
	// -- crowns
	chkEnableCrowns->setChecked(m_config.doCrowns);
	spnCrownsHeightFraction->setValue(m_config.crownsHeightFraction);
	spnCrownsRadius->setValue(m_config.crownsRadius);
	spnCrownsMinHeight->setValue(m_config.crownsMinHeight);
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	txtCrownsCrownsRaster->setText(qstr(m_config.crownsCrownsRaster));
	txtCrownsCrownsDatabase->setText(qstr(m_config.crownsCrownsDatabase));

	// Connect events
	// -- section toggles
	connect(chkEnableSmoothing, SIGNAL(toggled(bool)), this, SLOT(doSmoothChanged(bool)));
	connect(chkEnableTops, SIGNAL(toggled(bool)), this, SLOT(doTopsChanged(bool)));
	connect(chkEnableCrowns, SIGNAL(toggled(bool)), this, SLOT(doCrownsChanged(bool)));
	// -- smoothing
	connect(spnSmoothWindow, SIGNAL(valueChanged(int)), this, SLOT(smoothWindowSizeChanged(int)));
	connect(spnSmoothSigma, SIGNAL(valueChanged(double)), this, SLOT(smoothSigmaChanged(double)));
	connect(txtSmoothOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothOriginalCHMChanged(QString)));
	connect(txtSmoothSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(smoothSmoothedCHMChanged(QString)));
	connect(btnSmoothOriginalCHM, SIGNAL(clicked()), this, SLOT(smoothOriginalCHMClicked()));
	connect(btnSmoothSmoothedCHM, SIGNAL(clicked()), this, SLOT(smoothSmoothedCHMClicked()));
	// -- tops
	connect(spnTopsMinHeight, SIGNAL(valueChanged(double)), this, SLOT(topsMinHeightChanged(double)));
	connect(spnTopsWindowSize, SIGNAL(valueChanged(int)), this, SLOT(topsWindowSizeChanged(int)));
	connect(spnTopsTreetopsSRID, SIGNAL(valueChanged(int)), this, SLOT(topsTreetopsSRIDChanged(int)));
	connect(txtTopsOriginalCHM, SIGNAL(textChanged(QString)), this, SLOT(topsOriginalCHMChanged(QString)));
	connect(txtTopsSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(topsSmoothedCHMChanged(QString)));
	connect(txtTopsTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(topsTreetopsDatabaseChanged(QString)));
	connect(btnTopsOriginalCHM, SIGNAL(clicked()), this, SLOT(topsOriginalCHMClicked()));
	connect(btnTopsSmoothedCHM, SIGNAL(clicked()), this, SLOT(topsSmoothedCHMClicked()));
	connect(btnTopsTreetopsDatabase, SIGNAL(clicked()), this, SLOT(topsTreetopsDatabaseClicked()));
	// -- crowns
	connect(spnCrownsRadius, SIGNAL(valueChanged(double)), this, SLOT(crownsRadiusChanged(double)));
	connect(spnCrownsHeightFraction, SIGNAL(valueChanged(double)), this, SLOT(crownsHeightFractionChanged(double)));
	connect(spnCrownsMinHeight, SIGNAL(valueChanged(double)), this, SLOT(crownsMinHeightChanged(double)));
	connect(txtCrownsSmoothedCHM, SIGNAL(textChanged(QString)), this, SLOT(crownsSmoothedCHMChanged(QString)));
	connect(txtCrownsTreetopsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsTreetopsDatabaseChanged(QString)));
	connect(txtCrownsCrownsRaster, SIGNAL(textChanged(QString)), this, SLOT(crownsCrownsRasterChanged(QString)));
	connect(txtCrownsCrownsDatabase, SIGNAL(textChanged(QString)), this, SLOT(crownsCrownsDatabaseChanged(QString)));
	connect(btnCrownsSmoothedCHM, SIGNAL(clicked()), this, SLOT(crownsSmoothedCHMClicked()));
	connect(btnCrownsTreetopsDatabase, SIGNAL(clicked()), this, SLOT(crownsTreetopsDatabaseClicked()));
	connect(btnCrownsCrownsRaster, SIGNAL(clicked()), this, SLOT(crownsCrownsRasterClicked()));
	connect(btnCrownsCrownsDatabase, SIGNAL(clicked()), this, SLOT(crownsCrownsDatabaseClicked()));
	connect(btnTopsTreetopsSRID, SIGNAL(clicked()), this, SLOT(topsTreetopsSRIDClicked()));
	// -- program buttons
	connect(btnExit, SIGNAL(clicked()), this, SLOT(exitClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(helpClicked()));
	// -- callbacks
	if (m_callbacks) {
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(statusUpdate(QString)), lblStatus, SLOT(setText(QString)));
	}
	// -- worker thread.
	connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

	checkRun();
}

void TreetopsForm::resetProgress() {
	prgStep->setValue(0);
	prgOverall->setValue(0);
	lblStatus->setText("[Not Started]");
}

void TreetopsForm::topsTreetopsSRIDChanged(int srid) {
	m_config.srid = srid;
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
	m_config.smoothOriginalCHM = getInputFile(m_form, "CHM for Smoothing", m_last, RASTER_PATTERN);
	txtSmoothOriginalCHM->setText(qstr(m_config.smoothOriginalCHM));
	checkRun();
	if (m_config.topsOriginalCHM.empty()) {
		m_config.topsOriginalCHM = m_config.smoothOriginalCHM;
		txtTopsOriginalCHM->setText(qstr(m_config.topsOriginalCHM));
	}
	updateView();
}

void TreetopsForm::smoothSmoothedCHMClicked() {
	m_config.smoothSmoothedCHM = getOutputFile(m_form, "Smoothed CHM", m_last, RASTER_PATTERN);
	txtSmoothSmoothedCHM->setText(qstr(m_config.smoothSmoothedCHM));
	checkRun();
	if (m_config.topsSmoothedCHM.empty()) {
		m_config.topsSmoothedCHM = m_config.smoothSmoothedCHM;
		txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	}
	if (m_config.crownsSmoothedCHM.empty()) {
		m_config.crownsSmoothedCHM = m_config.smoothSmoothedCHM;
		txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	}
}

void TreetopsForm::topsSmoothedCHMClicked() {
	m_config.topsSmoothedCHM = getInputFile(m_form, "Smoothed CHM for Treetops", m_last, RASTER_PATTERN);
	txtTopsSmoothedCHM->setText(qstr(m_config.topsSmoothedCHM));
	checkRun();
	if (m_config.crownsSmoothedCHM.empty()) {
		m_config.crownsSmoothedCHM = m_config.topsSmoothedCHM;
		txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	}
}

void TreetopsForm::topsOriginalCHMClicked() {
	m_config.topsOriginalCHM = getInputFile(m_form, "Original CHM for Treetops", m_last, RASTER_PATTERN);
	txtTopsOriginalCHM->setText(qstr(m_config.topsOriginalCHM));
	checkRun();
}

void TreetopsForm::topsTreetopsDatabaseClicked() {
	m_config.topsTreetopsDatabase = getOutputFile(m_form, "Treetops Database", m_last, VECTOR_PATTERN);
	txtTopsTreetopsDatabase->setText(qstr(m_config.topsTreetopsDatabase));
	checkRun();
	if (m_config.crownsTreetopsDatabase.empty()) {
		m_config.crownsTreetopsDatabase = m_config.topsTreetopsDatabase;
		txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	}
}

void TreetopsForm::crownsSmoothedCHMClicked() {
	m_config.crownsSmoothedCHM = getInputFile(m_form, "Smoothed CHM for Crown Delineation", m_last, RASTER_PATTERN);
	txtCrownsSmoothedCHM->setText(qstr(m_config.crownsSmoothedCHM));
	checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseClicked() {
	m_config.crownsTreetopsDatabase = getInputFile(m_form, "Treetops Database", m_last, VECTOR_PATTERN);
	txtCrownsTreetopsDatabase->setText(qstr(m_config.crownsTreetopsDatabase));
	checkRun();
}

void TreetopsForm::crownsCrownsRasterClicked() {
	m_config.crownsCrownsRaster = getOutputFile(m_form, "Crowns Raster", m_last, RASTER_PATTERN);
	txtCrownsCrownsRaster->setText(qstr(m_config.crownsCrownsRaster));
	checkRun();
}

void TreetopsForm::crownsCrownsDatabaseClicked() {
	m_config.crownsCrownsDatabase = getOutputFile(m_form, "Crowns Database", m_last, VECTOR_PATTERN);
	txtCrownsCrownsDatabase->setText(qstr(m_config.crownsCrownsDatabase));
	checkRun();
}

void TreetopsForm::doSmoothChanged(bool v) {
	m_config.doSmoothing = v;
	enableGroup(m_smoothGroup, v);
	checkRun();
}

void TreetopsForm::doTopsChanged(bool v) {
	m_config.doTops = v;
	enableGroup(m_topsGroup, v);
	checkRun();
}

void TreetopsForm::doCrownsChanged(bool v) {
	m_config.doCrowns = v;
	enableGroup(m_crownsGroup, v);
	checkRun();
}

void TreetopsForm::crownsRadiusChanged(double radius) {
	m_config.crownsRadius = radius;
	checkRun();
}

void TreetopsForm::crownsHeightFractionChanged(double frac) {
	m_config.crownsHeightFraction = frac;
	checkRun();
}

void TreetopsForm::crownsMinHeightChanged(double height) {
	m_config.crownsMinHeight = height;
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
	m_config.crownsCrownsRaster = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsCrownsDatabaseChanged(QString file) {
	m_config.crownsCrownsDatabase = file.toStdString();
	checkRun();
}

void TreetopsForm::topsMinHeightChanged(double height) {
	m_config.topsMinHeight = height;
	checkRun();
}

void TreetopsForm::topsWindowSizeChanged(int size) {
	m_config.topsWindowSize = size;
	checkRun();
}

void TreetopsForm::topsSmoothedCHMChanged(QString file) {
	m_config.topsSmoothedCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::topsOriginalCHMChanged(QString file) {
	m_config.topsOriginalCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::topsTreetopsDatabaseChanged(QString file) {
	m_config.topsTreetopsDatabase = file.toStdString();
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
	m_config.smoothSmoothedCHM = file.toStdString();
	checkRun();
}

void TreetopsForm::runClicked() {
	if (m_workerThread->isRunning())
		return;
	m_cancel = false;
	btnRun->setEnabled(false);
	btnCancel->setEnabled(true);
	btnExit->setEnabled(false);
	m_workerThread->init(this);
	m_workerThread->start();
	resetProgress();
	checkRun();
}

void TreetopsForm::done() {
	if (m_workerThread->isError()) {
		errorDialog(m_form, "Error", m_workerThread->message());
		resetProgress();
	}
	checkRun();
}

void TreetopsForm::exitClicked() {
	g_debug("quit");
	m_form->close();
}

void TreetopsForm::cancelClicked() {
	g_debug("cancel");
	m_cancel = true;
	checkRun();
}

void TreetopsForm::helpClicked() {
	g_debug("help");
	QDesktopServices::openUrl(QUrl("http://www.dijital.ca/geotools/help/treetops.html", QUrl::TolerantMode));
}

void TreetopsForm::checkRun() {
	btnRun->setEnabled(m_config.canRun() && !m_workerThread->isRunning());
	btnCancel->setEnabled(m_workerThread->isRunning());
	btnExit->setEnabled(!m_workerThread->isRunning());
}
