#include <QWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QSettings>
#include <QDesktopServices>
#include <QUrl>

#include "geotools.h"
#include "treetops.hpp"
#include "treetops_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::treetops;
using namespace geotools::treetops::config;

QSettings _settings("Treetops", "GeoTools");

QString _str(const std::string &str) {
    return QString(str.c_str());
}

std::string _getInputFile(QWidget *form, const std::string &title, QDir &path, const std::string &filter) {
    QString res = QFileDialog::getOpenFileName(form, _str(title), path.path(), _str(filter));
    path.setPath(res);
    return res.toStdString();
}

std::string _getOutputFile(QWidget *form, const std::string &title, QDir &path, const std::string &filter) {
    QString res = QFileDialog::getSaveFileName(form, _str(title), path.path(), _str(filter));
    path.setPath(res);
    return res.toStdString();
}

void _loadConfig(TreetopsConfig &config) {
    QSettings qs("TreetopsConfig", "GeoTools");
    config.srid = qs.value(QString("srid"), config.srid).toInt();
    config.buildIndex = qs.value(QString("buildIndex"), config.buildIndex).toBool();
    config.tableCacheSize = qs.value(QString("tableCacheSize"), config.tableCacheSize).toInt();
    config.rowCacheSize = qs.value(QString("rowCacheSize"), config.rowCacheSize).toInt();
    config.doSmoothing = qs.value(QString("doSmoothing"), config.doSmoothing).toBool();
    config.smoothWindowSize = qs.value(QString("smoothingWindowSize"), config.smoothWindowSize).toInt();
    config.smoothSigma = qs.value(QString("smoothSigma"), config.smoothSigma).toDouble();
    config.smoothOriginalCHM = qs.value(QString("smoothOriginalCHM"), _str(config.smoothOriginalCHM)).toString().toStdString();
    config.smoothSmoothedCHM = qs.value(QString("smoothSmoothedCHM"), _str(config.smoothSmoothedCHM)).toString().toStdString();
    config.doTops = qs.value(QString("doTops"), config.doTops).toBool();
    config.topsMinHeight = qs.value(QString("topsMinHeight"), config.topsMinHeight).toDouble();
    config.topsWindowSize = qs.value(QString("topsWindowSize"), config.topsWindowSize).toInt();
    config.topsOriginalCHM = qs.value(QString("topsOriginalCHM"), _str(config.topsOriginalCHM)).toString().toStdString();
    config.topsSmoothedCHM = qs.value(QString("topsSmoothedCHM"), _str(config.topsSmoothedCHM)).toString().toStdString();
    config.topsTreetopsDatabase = qs.value(QString("topsTreetopsDatabase"), _str(config.topsTreetopsDatabase)).toString().toStdString();
    config.doCrowns = qs.value(QString("doCrowns"), config.doCrowns).toBool();
    config.crownsRadius = qs.value(QString("crownsRadius"), config.crownsRadius).toDouble();
    config.crownsHeightFraction = qs.value(QString("crownsHeightFraction"), config.crownsHeightFraction).toDouble();
    config.crownsMinHeight = qs.value(QString("crownsMinHeight"), config.crownsMinHeight).toDouble();
    config.crownsSmoothedCHM = qs.value(QString("crownsSmoothedCHM"), _str(config.crownsSmoothedCHM)).toString().toStdString();
    config.crownsTreetopsDatabase = qs.value(QString("crownsTreetopsDatabase"), _str(config.crownsTreetopsDatabase)).toString().toStdString();
    config.crownsCrownsRaster = qs.value(QString("crownsCrownsRaster"), _str(config.crownsCrownsRaster)).toString().toStdString();
    config.crownsCrownsDatabase = qs.value(QString("crownsCrownsDatabase"), _str(config.crownsCrownsDatabase)).toString().toStdString();
}

void _saveConfig(TreetopsConfig &config) {
    QSettings qs("TreetopsConfig", "GeoTools");
    qs.setValue(QString("srid"), config.srid);
    qs.setValue(QString("buildIndex"), config.buildIndex);
    qs.setValue(QString("tableCacheSize"), config.tableCacheSize);
    qs.setValue(QString("rowCacheSize"), config.rowCacheSize); //(24 * 1024 * 1024),
    qs.setValue(QString("doSmoothing"), config.doSmoothing);
    qs.setValue(QString("smoothingWindowSize"), config.smoothWindowSize);
    qs.setValue(QString("smoothSigma"), config.smoothSigma);
    qs.setValue(QString("smoothOriginalCHM"), _str(config.smoothOriginalCHM));
    qs.setValue(QString("smoothSmoothedCHM"), _str(config.smoothSmoothedCHM));
    qs.setValue(QString("doTops"), config.doTops);
    qs.setValue(QString("topsMinHeight"), config.topsMinHeight);
    qs.setValue(QString("topsWindowSize"), config.topsWindowSize);
    qs.setValue(QString("topsOriginalCHM"), _str(config.topsOriginalCHM));
    qs.setValue(QString("topsSmoothedCHM"), _str(config.topsSmoothedCHM));
    qs.setValue(QString("topsTreetopsDatabase"), _str(config.topsTreetopsDatabase));
    qs.setValue(QString("doCrowns"), config.doCrowns);
    qs.setValue(QString("crownsRadius"), config.crownsRadius);
    qs.setValue(QString("crownsHeightFraction"), config.crownsHeightFraction);
    qs.setValue(QString("crownsMinHeight"), config.crownsMinHeight);
    qs.setValue(QString("crownsSmoothedCHM"), _str(config.crownsSmoothedCHM));
    qs.setValue(QString("crownsTreetopsDatabase"), _str(config.crownsTreetopsDatabase));
    qs.setValue(QString("crownsCrownsRaster"), _str(config.crownsCrownsRaster));
    qs.setValue(QString("crownsCrownsDatabase"), _str(config.crownsCrownsDatabase));
    
}
    
TreetopsForm::TreetopsForm(QWidget *p) :
    QWidget(p),
    m_callbacks(nullptr) {
}

TreetopsForm::~TreetopsForm() {
    _saveConfig(m_config);
    delete m_callbacks;
    if (m_workerThread) {
        m_workerThread->exit(0);
        delete m_workerThread;
    }
}

void TreetopsForm::setupUi(QWidget *form) {
    Ui::TreetopsForm::setupUi(form);

    m_form = form;

    _loadConfig(m_config);
    
    if (_settings.contains(QString("last_dir"))) {
        m_last.setPath(_settings.value(QString("last_dir")).toString());
    } else {
        m_last = QDir::home();
    }
    
    m_callbacks = new TreetopsCallbacks();
    m_workerThread = new WorkerThread();
    m_workerThread->init(this);
    
    chkEnableSmoothing->setChecked(m_config.doSmoothing);
    spnSmoothWindow->setValue(m_config.smoothWindowSize);
    spnSmoothSigma->setValue(m_config.smoothSigma);
    txtSmoothOriginalCHM->setText(_str(m_config.smoothOriginalCHM));
    txtSmoothSmoothedCHM->setText(_str(m_config.smoothSmoothedCHM));

    chkEnableTops->setChecked(m_config.doTops);
    spnTopsMinHeight->setValue(m_config.topsMinHeight);
    spnTopsWindowSize->setValue(m_config.topsWindowSize);
    txtTopsOriginalCHM->setText(_str(m_config.topsOriginalCHM));
    txtTopsSmoothedCHM->setText(_str(m_config.topsSmoothedCHM));
    txtTopsTreetopsDatabase->setText(_str(m_config.topsTreetopsDatabase));
    spnTopsTreetopsSRID->setValue(m_config.srid);
    
    chkEnableCrowns->setChecked(m_config.doCrowns);
    spnCrownsHeightFraction->setValue(m_config.crownsHeightFraction);
    spnCrownsRadius->setValue(m_config.crownsRadius);
    spnCrownsMinHeight->setValue(m_config.crownsMinHeight);
    txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    txtCrownsTreetopsDatabase->setText(_str(m_config.crownsTreetopsDatabase));
    txtCrownsCrownsRaster->setText(_str(m_config.crownsCrownsRaster));
    txtCrownsCrownsDatabase->setText(_str(m_config.crownsCrownsDatabase));
    
    connect(chkEnableSmoothing, SIGNAL(toggled(bool)), SLOT(doSmoothChanged(bool)));
    connect(chkEnableTops, SIGNAL(toggled(bool)), SLOT(doTopsChanged(bool)));
    connect(chkEnableCrowns, SIGNAL(toggled(bool)), SLOT(doCrownsChanged(bool)));

    connect(spnSmoothWindow, SIGNAL(valueChanged(int)), SLOT(smoothWindowSizeChanged(int)));
    connect(spnSmoothSigma, SIGNAL(valueChanged(double)), SLOT(smoothSigmaChanged(double)));
    connect(txtSmoothOriginalCHM, SIGNAL(textChanged(QString)), SLOT(smoothOriginalCHMChanged(QString)));
    connect(txtSmoothSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(smoothSmoothedCHMChanged(QString)));
    connect(btnSmoothOriginalCHM, SIGNAL(clicked()), SLOT(smoothOriginalCHMClicked()));
    connect(btnSmoothSmoothedCHM, SIGNAL(clicked()), SLOT(smoothSmoothedCHMClicked()));

    connect(spnTopsMinHeight, SIGNAL(valueChanged(double)), SLOT(topsMinHeightChanged(double)));
    connect(spnTopsWindowSize, SIGNAL(valueChanged(int)), SLOT(topsWindowSizeChanged(int)));
    connect(spnTopsTreetopsSRID, SIGNAL(valueChanged(int)), SLOT(topsTreetopsSRIDChanged(int)));
    connect(txtTopsOriginalCHM, SIGNAL(textChanged(QString)), SLOT(topsOriginalCHMChanged(QString)));
    connect(txtTopsSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(topsSmoothedCHMChanged(QString)));
    connect(txtTopsTreetopsDatabase, SIGNAL(textChanged(QString)), SLOT(topsTreetopsDatabaseChanged(QString)));
    connect(btnTopsOriginalCHM, SIGNAL(clicked()), SLOT(topsOriginalCHMClicked()));
    connect(btnTopsSmoothedCHM, SIGNAL(clicked()), SLOT(topsSmoothedCHMClicked()));
    connect(btnTopsTreetopsDatabase, SIGNAL(clicked()), SLOT(topsTreetopsDatabaseClicked()));

    connect(spnCrownsRadius, SIGNAL(valueChanged(double)), SLOT(crownsRadiusChanged(double)));
    connect(spnCrownsHeightFraction, SIGNAL(valueChanged(double)), SLOT(crownsHeightFractionChanged(double)));
    connect(spnCrownsMinHeight, SIGNAL(valueChanged(double)), SLOT(crownsMinHeightChanged(double)));
    connect(txtCrownsSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(crownsSmoothedCHMChanged(QString)));
    connect(txtCrownsTreetopsDatabase, SIGNAL(textChanged(QString)), SLOT(crownsTreetopsDatabaseChanged(QString)));
    connect(txtCrownsCrownsRaster, SIGNAL(textChanged(QString)), SLOT(crownsCrownsRasterChanged(QString)));
    connect(txtCrownsCrownsDatabase, SIGNAL(textChanged(QString)), SLOT(crownsCrownsDatabaseChanged(QString)));
    connect(btnCrownsSmoothedCHM, SIGNAL(clicked()), SLOT(crownsSmoothedCHMClicked()));
    connect(btnCrownsTreetopsDatabase, SIGNAL(clicked()), SLOT(crownsTreetopsDatabaseClicked()));
    connect(btnCrownsCrownsRaster, SIGNAL(clicked()), SLOT(crownsCrownsRasterClicked()));
    connect(btnCrownsCrownsDatabase, SIGNAL(clicked()), SLOT(crownsCrownsDatabaseClicked()));
    connect(btnTopsTreetopsSRID, SIGNAL(clicked()), SLOT(topsTreetopsSRIDClicked()));

    connect(btnExit, SIGNAL(clicked()), SLOT(exitClicked()));
    connect(btnRun, SIGNAL(clicked()), SLOT(runClicked()));
    connect(btnCancel, SIGNAL(clicked()), SLOT(cancelClicked()));
    connect(btnHelp, SIGNAL(clicked()), SLOT(helpClicked()));
    
    if (m_callbacks) {
        connect((TreetopsCallbacks *) m_callbacks, SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
        connect((TreetopsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
    }
    connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

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
    m_config.smoothOriginalCHM = _getInputFile(this, "CHM for Smoothing", m_last, "GeoTiff (*.tif *.tiff)");
    txtSmoothOriginalCHM->setText(_str(m_config.smoothOriginalCHM));
    checkRun();
    if (m_config.topsOriginalCHM.empty()) {
        m_config.topsOriginalCHM = m_config.smoothOriginalCHM;
        txtTopsOriginalCHM->setText(_str(m_config.topsOriginalCHM));
    }
    updateView();
}

void TreetopsForm::smoothSmoothedCHMClicked() {
    m_config.smoothSmoothedCHM = _getOutputFile(this, "Smoothed CHM", m_last, "GeoTiff (*.tif *.tiff)");
    txtSmoothSmoothedCHM->setText(_str(m_config.smoothSmoothedCHM));
    checkRun();
    if (m_config.topsSmoothedCHM.empty()) {
        m_config.topsSmoothedCHM = m_config.smoothSmoothedCHM;
        txtTopsSmoothedCHM->setText(_str(m_config.topsSmoothedCHM));
    }
    if (m_config.crownsSmoothedCHM.empty()) {
        m_config.crownsSmoothedCHM = m_config.smoothSmoothedCHM;
        txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    }
}

void TreetopsForm::topsSmoothedCHMClicked() {
    m_config.topsSmoothedCHM = _getInputFile(this, "Smoothed CHM for Treetops", m_last, "GeoTiff (*.tif *.tiff)");
    txtTopsSmoothedCHM->setText(_str(m_config.topsSmoothedCHM));
    checkRun();
    if (m_config.crownsSmoothedCHM.empty()) {
        m_config.crownsSmoothedCHM = m_config.topsSmoothedCHM;
        txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    }
}

void TreetopsForm::topsOriginalCHMClicked() {
    m_config.topsOriginalCHM = _getInputFile(this, "Original CHM for Treetops", m_last, "GeoTiff (*.tif *.tiff)");
    txtTopsOriginalCHM->setText(_str(m_config.topsOriginalCHM));
    checkRun();
}

void TreetopsForm::topsTreetopsDatabaseClicked() {
    m_config.topsTreetopsDatabase = _getOutputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
    txtTopsTreetopsDatabase->setText(_str(m_config.topsTreetopsDatabase));
    checkRun();
    if (m_config.crownsTreetopsDatabase.empty()) {
        m_config.crownsTreetopsDatabase = m_config.topsTreetopsDatabase;
        txtCrownsTreetopsDatabase->setText(_str(m_config.crownsTreetopsDatabase));
    }
}

void TreetopsForm::crownsSmoothedCHMClicked() {
    m_config.crownsSmoothedCHM = _getInputFile(this, "Smoothed CHM for Crown Delineation", m_last, "GeoTiff (*.tif *.tiff)");
    txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseClicked() {
    m_config.crownsTreetopsDatabase = _getInputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
    txtCrownsTreetopsDatabase->setText(_str(m_config.crownsTreetopsDatabase));
    checkRun();
}

void TreetopsForm::crownsCrownsRasterClicked() {
    m_config.crownsCrownsRaster = _getOutputFile(this, "Crowns Raster", m_last, "GeoTiff (*.tif *.tiff)");
    txtCrownsCrownsRaster->setText(_str(m_config.crownsCrownsRaster));
    checkRun();
}

void TreetopsForm::crownsCrownsDatabaseClicked() {
    m_config.crownsCrownsDatabase = _getOutputFile(this, "Crowns Database", m_last, "SQLite (*.sqlite)");
    txtCrownsCrownsDatabase->setText(_str(m_config.crownsCrownsDatabase));
    checkRun();
}

void TreetopsForm::doSmoothChanged(bool v) {
    m_config.doSmoothing = v;
    checkRun();
}

void TreetopsForm::doTopsChanged(bool v) {
    m_config.doTops = v;
    checkRun();
}

void TreetopsForm::doCrownsChanged(bool v) {
    m_config.doCrowns = v;
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
    checkRun();
}

void TreetopsForm::done() {
    checkRun();
}

void TreetopsForm::exitClicked() {
    g_trace("quit");
    m_form->close();
}

void TreetopsForm::cancelClicked() {
    g_trace("cancel");
    m_cancel = true;
    checkRun();
}

void TreetopsForm::helpClicked() {
    g_trace("help");
    QDesktopServices::openUrl(QUrl("http://www.dijital.ca/geotools/help/treetops.html", QUrl::TolerantMode));
}

void TreetopsForm::checkRun() {
    btnRun->setEnabled(m_config.canRun() && !m_workerThread->isRunning());
    btnCancel->setEnabled(m_workerThread->isRunning());
    btnExit->setEnabled(!m_workerThread->isRunning());
}