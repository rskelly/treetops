/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <thread>

#include <QSettings>
#include <QFileDialog>
#include <QMessageBox>
#include <QDesktopServices>
#include <QUrl>

#include "pointnormalize.hpp"
#include "pointnormalize_ui.hpp"
#include "filelist.hpp"
#include "util.hpp"

QSettings _settings("PointNormalize", "Geotools");
QString _last_dir("last_dir");

using namespace geotools::ui;
using namespace geotools::util;

void _loadConfig(PointNormalizeConfig &config) {
    QSettings qs("PointNormalizeConfig", "GeoTools");
    config.dropNegative = qs.value(QString("dropNegative"), true).toBool();
    config.dropGround = qs.value(QString("dropGround"), true).toBool();
    config.threads = qs.value(QString("threads"), 1).toInt();
    config.overwrite = qs.value(QString("overwrite"), true).toBool();
    config.buffer = qs.value(QString("buffer"), 10.0).toDouble();
    config.outputDir = qs.value(QString("outputDir")).toString().toStdString();
    QStringList files = qs.value(QString("sourceFiles")).toStringList();
    for(const QString &file : files)
        config.sourceFiles.push_back(file.toStdString());
}

void _saveConfig(PointNormalizeConfig &config) {
    QSettings qs("PointNormalizeConfig", "GeoTools");
    qs.setValue(QString("dropNegative"), config.dropNegative);
    qs.setValue(QString("dropGround"), config.dropGround);
    qs.setValue(QString("threads"), config.threads);
    qs.setValue(QString("overwrite"), config.overwrite);
    qs.setValue(QString("buffer"), config.buffer);
    qs.setValue(QString("outputDir"), QString(config.outputDir.c_str()));
    QStringList files;
    for(const std::string &file : config.sourceFiles)
        files << QString(file.c_str());
    qs.setValue(QString("sourceFiles"), files);
}

void PointNormalizeCallbacks::stepCallback(float status) const {
    emit stepProgress((int) std::round(status * 100));
}

void PointNormalizeCallbacks::overallCallback(float status) const {
    emit overallProgress((int) std::round(status * 100));
}

PointNormalizeForm::PointNormalizeForm(QWidget *p) {
}

void PointNormalizeForm::setupUi(QWidget *form) {
    Ui::PointNormalizeForm::setupUi(form);
    m_form = form;
    m_filter = QString("LAS Files (*.las)");
    
    _loadConfig(m_config);
    
    if (_settings.contains(_last_dir)) {
        m_last.setPath(_settings.value(_last_dir).toString());
    } else {
        m_last = QDir::home();
    }

    m_bgrpNegative = new QButtonGroup();
    m_bgrpNegative->addButton(rdoKeepNegativePoints);
    m_bgrpNegative->addButton(rdoRemoveNegativePoints);
    m_bgrpGround = new QButtonGroup();
    m_bgrpGround->addButton(rdoKeepGroundPoints);
    m_bgrpGround->addButton(rdoRemoveGroundPoints);
    rdoKeepNegativePoints->setChecked(!m_config.dropNegative);
    rdoKeepGroundPoints->setChecked(!m_config.dropGround);
    
    m_fileList.init(this, btnAddFiles, btnRemoveAllFiles, btnRemoveSelectedFiles, 
            lstFiles, m_last, m_filter);
   
    m_callbacks = new PointNormalizeCallbacks();
    m_workerThread = new WorkerThread();
    m_workerThread->init(this);
    
    spnThreads->setMaximum(std::thread::hardware_concurrency());
    spnThreads->setValue(m_config.threads);
    spnBuffer->setValue(m_config.buffer);
    chkOverwrite->setChecked(m_config.overwrite);
    txtOutputFolder->setText(QString(m_config.outputDir.c_str()));
    m_fileList.setFiles(m_config.sourceFiles);
    
    connect(btnOutputFolder, SIGNAL(clicked()), this, SLOT(outputFolderClicked()));
    connect(btnExit, SIGNAL(clicked()), this, SLOT(exitClicked()));
    connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
    connect(btnHelp, SIGNAL(clicked()), this, SLOT(helpClicked()));
    connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
    connect(rdoKeepNegativePoints, SIGNAL(toggled(bool)), this, SLOT(keepNegativeToggled(bool)));
    connect(rdoKeepGroundPoints, SIGNAL(toggled(bool)), this, SLOT(keepGroundToggled(bool)));
    connect(&m_fileList, SIGNAL(fileListChanged()), this, SLOT(fileListChanged()));
    connect((PointNormalizeCallbacks *) m_callbacks, SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
    connect((PointNormalizeCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
    connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));
    connect(chkOverwrite, SIGNAL(toggled(bool)), this, SLOT(overwriteChanged(bool)));
    connect(spnThreads, SIGNAL(valueChanged(int)), this, SLOT(threadsChanged(int)));
    connect(spnBuffer, SIGNAL(valueChanged(double)), this, SLOT(bufferChanged(double)));
}

void PointNormalizeForm::bufferChanged(double b) {
    m_config.buffer = b;
    g_debug(" -- buffer " << b);
    checkRun();
}

void PointNormalizeForm::threadsChanged(int t) {
    m_config.threads = t;
    g_debug(" -- threads " << t);
    checkRun();
}

void PointNormalizeForm::overwriteChanged(bool o) {
    m_config.overwrite = o;
    g_debug(" -- overwrite " << o);
    checkRun();
}

void PointNormalizeForm::keepNegativeToggled(bool on) {
    m_config.dropNegative = !on;
    g_debug(" -- drop negative " << m_config.dropNegative);
    checkRun();
}

void PointNormalizeForm::keepGroundToggled(bool on) {
    m_config.dropGround = !on;
    g_debug(" -- drop ground " << m_config.dropGround);
    checkRun();
}

const PointNormalizeConfig& PointNormalizeForm::config() {
    return m_config;
}

Callbacks* PointNormalizeForm::callbacks() {
    return m_callbacks;
}

bool* PointNormalizeForm::cancel() {
    return &m_cancel;
}

void PointNormalizeForm::fileListChanged() {
    m_config.sourceFiles = m_fileList.files();
    checkRun();
}

void PointNormalizeForm::checkRun() {
    btnRun->setEnabled(!m_workerThread->isRunning() && 
        m_config.sourceFiles.size() > 0 && !m_config.outputDir.empty());
    btnCancel->setEnabled(m_workerThread->isRunning());
    btnExit->setEnabled(!m_workerThread->isRunning());
}

void PointNormalizeForm::outputFolderClicked() {
    QString res = QFileDialog::getExistingDirectory(this, "Output Directory", 
            m_last.path(),QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    m_config.outputDir = res.toStdString();
    m_last.setPath(res);
    _settings.setValue(_last_dir, m_last.path());
    txtOutputFolder->setText(res);
    checkRun();
}

void PointNormalizeForm::exitClicked() {
    m_form->close();
}

void PointNormalizeForm::cancelClicked() {
    m_cancel = true;
    prgStep->setValue(0);
    prgOverall->setValue(0);
    checkRun();
}

void PointNormalizeForm::helpClicked() {
  QDesktopServices::openUrl(QUrl("http://www.dijital.ca/geotools/help/pointnormalize.html", QUrl::TolerantMode));
}

void PointNormalizeForm::runClicked() {
    if (m_workerThread->isRunning())
        return;
    m_cancel = false;
    btnRun->setEnabled(false);
    btnCancel->setEnabled(true);
    btnExit->setEnabled(false);
    m_workerThread->start();
    checkRun();
}

void PointNormalizeForm::done() {
    g_debug(" -- done");
    if (m_workerThread->hasError()) {
        QMessageBox err((QWidget *) this);
        err.setText("Error");
        err.setInformativeText(QString(m_workerThread->getError().c_str()));
        err.exec();
    }
    checkRun();
}

PointNormalizeForm::~PointNormalizeForm() {
   _saveConfig(m_config);
}


void WorkerThread::init(PointNormalizeForm *parent) {
    m_parent = parent;
}

void WorkerThread::run() {
    try {
        PointNormalize pn;
        pn.normalize(m_parent->config(), m_parent->callbacks(), m_parent->cancel());
    } catch(const std::exception &ex) {
        m_error = ex.what();
    }
}

bool WorkerThread::hasError() {
    return !m_error.empty();
}

std::string WorkerThread::getError() {
    return m_error;
}