#include <QWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QSettings>

#include "geotools.h"
#include "pointstats.hpp"
#include "pointstats_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::point;
using namespace geotools::point::pointstats_config;

QSettings _settings("PointStats", "GeoTools");

std::string _str(const QVariant &v) {
    return v.toString().toStdString();
}

QString _str(const std::string &s) {
    return s.c_str();
}

void _loadConfig(PointStatsConfig &config) {
    PointStatsConfig dummy;
    QSettings qs("PointStatsConfig", "GeoTools");
    config.angleLimit = qs.value("angleLimit", dummy.angleLimit).toInt();
    config.attribute = qs.value("attribute", dummy.attribute).toInt();
    config.bounds.fromString(_str(qs.value("bounds", _str(dummy.bounds.toString()))));
    config.fill = qs.value("fill", dummy.fill).toBool();
    config.gapFractionType = qs.value("gapFractionType", dummy.gapFractionType).toInt();
    config.gapThreshold = qs.value("gapThreshold", dummy.gapThreshold).toDouble();
    config.hsrid = qs.value("hsrid", dummy.hsrid).toInt();
    config.normalize = qs.value("normalize", dummy.normalize).toBool();
    config.quantile = qs.value("quantile", dummy.quantile).toInt();
    config.quantileFilter = qs.value("quantileFilter", dummy.quantileFilter).toInt();
    config.quantileFilterFrom = qs.value("quantileFilterFrom", dummy.quantileFilterFrom).toInt();
    config.quantileFilterTo = qs.value("quantileFilterTo", dummy.quantileFilterTo).toInt();
    config.quantiles = qs.value("quantiles", dummy.quantiles).toInt();
    config.rebuild = qs.value("rebuild", dummy.rebuild).toBool();
    config.resolution = qs.value("resolution", dummy.resolution).toDouble();
    config.snap = qs.value("snap", dummy.snap).toBool();
    config.threads = qs.value("threads", dummy.threads).toInt();
    config.vsrid = qs.value("vsid", dummy.vsrid).toInt();
    QStringList files = qs.value("sourceFiles").toStringList();
    for(const QString &file : files)
        config.sourceFiles.push_back(file.toStdString());
    QList<QVariant> types = qs.value("types").toList();
    for(const QVariant &type : types)
        config.types.push_back(type.toInt());
    QList<QVariant> classes = qs.value("classes").toList();
    for(const QVariant &cls : classes)
        config.classes.insert(cls.toInt());
    QStringList dstFiles = qs.value("dstFiles").toStringList();
    for(const QString &file : dstFiles)
        config.dstFiles.push_back(_str(file));
}

void _saveConfig(PointStatsConfig &config) {
    QSettings qs("PointStatsConfig", "GeoTools");
    qs.setValue("angleLimit", config.angleLimit);
    qs.setValue("attribute", config.attribute);
    qs.setValue("bounds", _str(config.bounds.toString()));
    qs.setValue("fill", config.fill);
    qs.setValue("gapFractionType", config.gapFractionType);
    qs.setValue("gapThreshold", config.gapThreshold);
    qs.setValue("hsrid", config.hsrid);
    qs.setValue("normalize", config.normalize);
    qs.setValue("quantile", config.quantile);
    qs.setValue("quantileFilter", config.quantileFilter);
    qs.setValue("quantileFilterFrom", config.quantileFilterFrom);
    qs.setValue("quantileFilterTo", config.quantileFilterTo);
    qs.setValue("quantiles", config.quantiles);
    qs.setValue("rebuild", config.rebuild);
    qs.setValue("resolution", config.resolution);
    qs.setValue("snap", config.snap);
    qs.setValue("threads", config.threads);
    qs.setValue("vsid", config.vsrid);
    QList<QVariant> classes;
    for(const uint8_t &cls : config.classes)
        classes << QVariant(cls);
    qs.setValue("classes", classes);
    QStringList dstFiles;
    for(const std::string &file : config.dstFiles)
        dstFiles << file.c_str();
    qs.setValue("dstFiles", dstFiles);
    QStringList sourceFiles;
    for(const std::string &file : config.sourceFiles)
        sourceFiles << file.c_str();
    qs.setValue("sourceFiles", sourceFiles);
    QList<QVariant> types;
    for(const uint8_t &type : config.types)
        types << QVariant(type);
    qs.setValue("types", types);
}

void PointStatsCallbacks::stepCallback(float status) const {
    emit stepProgress((int) std::round(status * 100));
}

void PointStatsCallbacks::overallCallback(float status) const {
    emit overallProgress((int) std::round(status * 100));
}

void WorkerThread::init(PointStatsForm *parent, const geotools::util::Bounds &bounds) {
    m_parent = parent;
    m_bounds = bounds;
}

void WorkerThread::run() {
    try {
        m_error = "";
        geotools::point::PointStats l;
        l.pointstats(m_parent->m_config, m_parent->m_callbacks);
    } catch (const std::exception &e) {
        m_error = e.what();
    }
}

bool WorkerThread::hasError() {
    return !m_error.empty();
}

std::string WorkerThread::getError() {
    return m_error;
}

PointStatsForm::PointStatsForm(QWidget *p) :
    QWidget(p),
    m_form(nullptr), 
    m_last(QDir::home()),
    m_workerThread(nullptr),
    m_callbacks(nullptr) {
}

PointStatsForm::~PointStatsForm() {
    delete m_callbacks;
    if (m_workerThread) {
        m_workerThread->exit(0);
        delete m_workerThread;
    }
    _saveConfig(m_config);
}

void PointStatsForm::setupUi(QWidget *form) {

    Ui::PointStatsForm::setupUi(form);

    m_form = form;
    if (_settings.contains("last_dir")) {
        m_last.setPath(_settings.value("last_dir").toString());
    } else {
        m_last = QDir::home();
    }

    m_workerThread = new WorkerThread();
    m_callbacks = new PointStatsCallbacks();

    _loadConfig(m_config);
    
    g_debug("x");
    
    spnResolution->setValue(m_config.resolution);
    spnMaxAngle->setValue(m_config.angleLimit);
    spnThreads->setValue(m_config.threads);
    spnThreads->setMaximum(g_max(1, omp_get_num_procs()));
    chkSnapToGrid->setChecked(m_config.snap);
    spnQuantileFilter->setValue(m_config.quantileFilter);
    spnQuantileFilterFrom->setValue(m_config.quantileFilterFrom);
    spnQuantileFilterTo->setValue(m_config.quantileFilterTo);
    spnQuantiles->setValue(m_config.quantiles);
    spnQuantile->setValue(m_config.quantile);
    spnOriginX->setValue(m_config.originX);
    spnOriginY->setValue(m_config.originY);
    
    int i = 0;
    int defaultIdx = -1;
    for (const auto &it : types) {
        cboType->addItem(it.first.c_str(), QVariant(it.second));
        if (std::find(m_config.types.begin(), m_config.types.end(), it.second) != m_config.types.end()) {
            defaultIdx = i;
            break;
        }
        ++i;
    }
    cboType->setCurrentIndex(defaultIdx);

    i = 0;
    defaultIdx = -1;
    for (const auto &it : gapFractionTypes) {
        cboGapFunction->addItem(it.first.c_str(), QVariant(it.second));
        if (it.second == m_config.gapFractionType)
            defaultIdx = i;
        ++i;
    }
    cboGapFunction->setCurrentIndex(defaultIdx);

    i = 0;
    defaultIdx = -1;
    for (const auto &it : attributes) {
        cboAttribute->addItem(it.first.c_str(), QVariant(it.second));
        if (it.second == m_config.attribute)
            defaultIdx = i;
        ++i;
    }
    cboAttribute->setCurrentIndex(defaultIdx);

    for (i = 0; i < 255; ++i) {
        QString str;
        str.setNum(i);
        QListWidgetItem *item = new QListWidgetItem(str, lstClasses);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(m_config.classes.find(i) != m_config.classes.end() ? Qt::Checked : Qt::Unchecked);
        lstClasses->addItem(item);
    }

    connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFilesClicked()));
    connect(btnRemoveSelected, SIGNAL(clicked()), this, SLOT(removeFilesClicked()));
    connect(btnClearFiles, SIGNAL(clicked()), this, SLOT(clearFilesClicked()));
    connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
    connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
    connect(btnDestFile, SIGNAL(clicked()), this, SLOT(destFileClicked()));
    connect(lstFiles, SIGNAL(itemSelectionChanged()), this, SLOT(fileListSelectionChanged()));
    connect(lstClasses, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(classItemClicked(QListWidgetItem*)));
    connect(btnCRSConfig, SIGNAL(clicked()), this, SLOT(crsConfigClicked()));
    connect(cboType, SIGNAL(currentIndexChanged(int)), this, SLOT(typeSelected(int)));
    connect(spnResolution, SIGNAL(valueChanged(double)), this, SLOT(resolutionChanged(double)));
    connect(chkSnapToGrid, SIGNAL(toggled(bool)), this, SLOT(snapToGridChanged(bool)));
    connect(cboAttribute, SIGNAL(currentIndexChanged(int)), this, SLOT(attributeSelected(int)));
    connect(cboGapFunction, SIGNAL(currentIndexChanged(int)), this, SLOT(gapFunctionSelected(int)));
    connect(spnThreads, SIGNAL(valueChanged(int)), this, SLOT(threadsChanged(int)));
    connect(spnQuantile, SIGNAL(valueChanged(int)), this, SLOT(quantileChanged(int)));
    connect(spnQuantiles, SIGNAL(valueChanged(int)), this, SLOT(quantilesChanged(int)));
    connect(spnMaxAngle, SIGNAL(valueChanged(int)),this,  SLOT(maxAngleChanged(int)));
    connect(spnOriginX, SIGNAL(valueChanged(double)), this, SLOT(originXChanged(double)));
    connect(spnOriginY, SIGNAL(valueChanged(double)), this, SLOT(originYChanged(double)));
    connect(spnQuantileFilter, SIGNAL(valueChanged(int)), this, SLOT(quantileFilterChanged(int)));
    connect(spnQuantileFilterFrom, SIGNAL(valueChanged(int)), this, SLOT(quantileFilterFromChanged(int)));
    connect(spnQuantileFilterTo, SIGNAL(valueChanged(int)), this, SLOT(quantileFilterToChanged(int)));
    
    connect((PointStatsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
    connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

    updateTypeUi();
}

void PointStatsform::originXChanged(double x) {
    m_config.originX = x;
    g_debug(" -- originX " << x);
    checkRun();
}

void PointStatsform::originYChanged(double y) {
    m_config.originY = y;
    g_debug(" -- originY " << y);
    checkRun();
}

void PointStatsForm::quantileFilterChanged(int quantiles) {
    m_config.quantileFilter = quantiles;
    g_debug(" -- quantileFilter " << quantiles);
    checkRun();
}

void PointStatsForm::quantileFilterFromChanged(int quantiles) {
    m_config.quantileFilterFrom = quantiles;
    g_debug(" -- quantileFilterFrom " << quantiles);
    checkRun();
}

void PointStatsForm::quantileFilterToChanged(int quantiles) {
    m_config.quantileFilterTo = quantiles;
    g_debug(" -- quantileFilterTo " << quantiles);
    checkRun();
}

void PointStatsForm::classItemClicked(QListWidgetItem *item) {
    unsigned char c = (unsigned char) item->text().toUShort();
    if (item->checkState() == Qt::Checked) {
        m_config.classes.insert(c);
    } else {
        m_config.classes.erase(c);
    }
    g_debug(" -- classes " << m_config.classes.size());
}

void PointStatsForm::threadsChanged(int threads) {
    m_config.threads = threads;
    g_debug(" -- threads " << m_config.threads);
    checkRun();
}

void PointStatsForm::maxAngleChanged(int q) {
    m_config.angleLimit = q;
    g_debug(" -- angle " << m_config.angleLimit);
    checkRun();
}

void PointStatsForm::quantileChanged(int q) {
    m_config.quantile = q;
    g_debug(" -- quantile " << m_config.quantile);
    checkRun();
}

void PointStatsForm::quantilesChanged(int q) {
    m_config.quantiles = q;
    g_debug(" -- quantiles " << m_config.quantiles);
    checkRun();
}

void PointStatsForm::attributeSelected(int index) {
    std::string att = cboAttribute->itemText(index).toStdString();
    m_config.attribute = attributes[att];
    g_debug(" -- att " << m_config.attribute);
    checkRun();
}

void PointStatsForm::gapFunctionSelected(int index) {
    std::string gap = cboGapFunction->itemText(index).toStdString();
    m_config.gapFractionType = gapFractionTypes[gap];
    g_debug(" -- gap type " << m_config.gapFractionType);
    checkRun();
}

void PointStatsForm::snapToGridChanged(bool state) {
    m_config.snap = state;
    g_debug(" -- snap " << m_config.snap);
    checkRun();
}

void PointStatsForm::resolutionChanged(double res) {
    m_config.resolution = res;
    g_debug(" -- resolution " << m_config.resolution);
    checkRun();
}

void PointStatsForm::updateTypeUi() {
    // TODO: See state machine framework
    if(!m_config.types.size()) return;
    uint8_t type = m_config.types[0];
    spnQuantile->setVisible(type == TYPE_QUANTILE);
    spnQuantiles->setVisible(type == TYPE_QUANTILE);
    lblQuantile->setVisible(type == TYPE_QUANTILE);
    lblQuantiles->setVisible(type == TYPE_QUANTILE);
    cboGapFunction->setVisible(type == TYPE_GAP_FRACTION);
    lblGapFunction->setVisible(type == TYPE_GAP_FRACTION);
    lblAttribute->setVisible(type != TYPE_GAP_FRACTION);
    cboAttribute->setVisible(type != TYPE_GAP_FRACTION);
}

void PointStatsForm::typeSelected(int index) {
    std::string type = cboType->itemText(index).toStdString();
    if(m_config.types.size()) {
        m_config.types[0] = types[type];
    } else {
        m_config.types.push_back(types[type]);
    }
    g_debug(" -- type " << m_config.types[0]);
    updateTypeUi();
    checkRun();
}

void PointStatsForm::crsConfigClicked() {
    CRSSelector cs(m_form);
    cs.setHorizontalSRID(m_config.hsrid);
    cs.setVerticalSRID(m_config.vsrid);
    if (cs.exec()) {
        m_config.vsrid = cs.getVerticalSRID();
        m_config.hsrid = cs.getHorizontalSRID();
        std::stringstream ss;
        ss << "epsg:" << m_config.hsrid;
        if (m_config.vsrid > 0)
            ss << "+" << m_config.vsrid;
        txtCRSConfig->setText(QString(ss.str().c_str()));
    }
    g_debug(" -- crs " << m_config.hsrid << "; " << m_config.vsrid);
    checkRun();
}

void PointStatsForm::fileListSelectionChanged() {
    updateFileButtons();
    checkRun();
}

void PointStatsForm::destFileClicked() {
    QString res = QFileDialog::getSaveFileName(this, "Save File", m_last.path(), "GeoTiff (*.tif *.tiff)");
    if(m_config.dstFiles.size()) {
        m_config.dstFiles[0] = res.toStdString();
    } else {
        m_config.dstFiles.push_back(res.toStdString());
    }
    m_last.setPath(res);
    _settings.setValue("last_dir", m_last.path());
    txtDestFile->setText(res);
    g_debug(" -- dest file " << m_config.dstFiles[0]);
    checkRun();
}

void PointStatsForm::runClicked() {

    if (m_workerThread->isRunning())
        return;

    btnRun->setEnabled(false);
    btnCancel->setEnabled(false);

    // TODO: Bounds
    Bounds bounds;
    m_workerThread->init(this, bounds);
    m_workerThread->start();

}

void PointStatsForm::done() {
    if (m_workerThread->hasError()) {
        QMessageBox err((QWidget *) this);
        err.setText("Error");
        err.setInformativeText(QString(m_workerThread->getError().c_str()));
        err.exec();
    }
    checkRun();
    btnCancel->setEnabled(true);
}

void PointStatsForm::cancelClicked() {
    g_trace("quit");
    m_form->close();
}

void PointStatsForm::updateFileList() {
    while (lstFiles->count())
        lstFiles->takeItem(0);
    for (const std::string &file : m_config.sourceFiles)
        lstFiles->addItem(QString(file.c_str()));
    updateFileButtons();
    checkRun();
}

void PointStatsForm::updateFileButtons() {
    btnClearFiles->setEnabled(lstFiles->count() > 0);
    btnRemoveSelected->setEnabled(lstFiles->selectedItems().size() > 0);
}

void PointStatsForm::removeFilesClicked() {
    std::vector<std::string> lst;
    unsigned int i = 0;
    for (const std::string &file : m_config.sourceFiles) {
        QListWidgetItem *item = lstFiles->item(i);
        if (!item->isSelected())
            lst.push_back(file);
        ++i;
    }
    m_config.sourceFiles.clear();
    m_config.sourceFiles.assign(lst.begin(), lst.end());
    g_debug(" -- source files " << m_config.sourceFiles.size());
    updateFileList();
    checkRun();
}

void PointStatsForm::clearFilesClicked() {
    m_config.sourceFiles.clear();
    g_debug(" -- source files " << m_config.sourceFiles.size());
    updateFileList();
    checkRun();
}

void PointStatsForm::selectFilesClicked() {
    QFileDialog d(m_form);
    d.setDirectory(m_last);
    d.setFileMode(QFileDialog::ExistingFiles);
    d.setNameFilter(QString("LAS Files (*.las)"));
    if (d.exec()) {
        QStringList files = d.selectedFiles();
        m_last = d.directory();
        _settings.setValue("last_dir", m_last.path());
        std::set<std::string> tmp(m_config.sourceFiles.begin(), m_config.sourceFiles.end());
        for (int i = 0; i < files.size(); ++i)
            tmp.insert(files[i].toStdString());
        m_config.sourceFiles.clear();
        m_config.sourceFiles.assign(tmp.begin(), tmp.end());
    }
    g_debug(" -- source files " << m_config.sourceFiles.size());
    updateFileList();
    checkRun();
}

void PointStatsForm::checkRun() {
    // TODO: Check runnable.
    btnRun->setEnabled(true);
}
