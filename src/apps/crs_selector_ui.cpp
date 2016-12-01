#include <iterator>
#include <cstdlib>

#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QCompleter>
#include <QKeyEvent>

#include "cpl_conv.h"

#include "geotools.h"
#include "util.hpp"
#include "csv.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::util;
using namespace geotools::csv;

CRSSelector::CRSSelector(QWidget *parent, Qt::WindowFlags f) :
    QDialog(parent, f),
    m_vsrid(0), m_hsrid(0) {
    initUi();
}

void CRSSelector::loadCrs(std::map<int, std::string> &target, const std::string &filename) {
    g_debug(" -- loadCRS");
    std::string path(CPLFindFile("gdal", filename.c_str()));
    CSVReader csv(path);
    std::unordered_map<std::string, std::string> row;
    while(csv.next(row)) {
        if(row.find("COORD_REF_SYS_CODE") == row.end() || row.find("COORD_REF_SYS_NAME") == row.end())
            g_runerr("Missing fields from CRS database.");
        int srid = atoi(row["COORD_REF_SYS_CODE"].c_str());
        target[srid] = row["COORD_REF_SYS_NAME"];
    }
}

void CRSSelector::enableVertical(bool e) {
    txtVerticalDatum->setEnabled(e);
}

void CRSSelector::enableHorizontal(bool e) {
    txtHorizontalCRS->setEnabled(e);
}

void CRSSelector::initUi() {
    Ui::CRSSelector::setupUi((QDialog *) this);

    loadCrs(m_vcrs, "vertcs.csv");
    loadCrs(m_vcrs, "vertcs.override.csv");
    m_vlst << QString("0 - None");
    for (const std::pair<int, std::string> &item : m_vcrs)
        m_vlst << QString((std::to_string(item.first) + " - " + item.second).c_str());

    m_vcomp = new QCompleter(m_vlst, this);
    m_vcomp->setCaseSensitivity(Qt::CaseInsensitive);
    m_vcomp->setCompletionMode(QCompleter::InlineCompletion);
    m_vcomp->setFilterMode(Qt::MatchContains);
    txtVerticalDatum->setCompleter(m_vcomp);

    loadCrs(m_hcrs, "gcs.csv");
    loadCrs(m_hcrs, "pcs.csv");
    m_hlst << QString("0 - None");
    for (const std::pair<int, std::string> &item : m_hcrs)
        m_hlst << QString((std::to_string(item.first) + " - " + item.second).c_str());

    m_hcomp = new QCompleter(m_hlst, this);
    m_hcomp->setCaseSensitivity(Qt::CaseInsensitive);
    m_hcomp->setCompletionMode(QCompleter::InlineCompletion);
    m_hcomp->setFilterMode(Qt::MatchContains);
    txtHorizontalCRS->setCompleter(m_hcomp);

    connect(btnSelect, SIGNAL(clicked()), this, SLOT(selectClicked()));
    connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
    connect(txtHorizontalCRS, SIGNAL(textChanged(const QString &)), this, SLOT(textUpdate()));
    connect(txtVerticalDatum, SIGNAL(textChanged(const QString &)), this, SLOT(textUpdate()));

    updateFields();
}

void CRSSelector::updateFields() {
    //if (m_hcomp && m_hsrid > 0) {
        QString h;
        h.setNum(m_hsrid);
        txtHorizontalCRS->setText(h);
        m_hcomp->complete();
    //}
    //if (m_vcomp && m_vsrid > 0) {
        QString v;
        v.setNum(m_vsrid);
        txtVerticalDatum->setText(v);
        m_vcomp->complete();
    //}
}

void CRSSelector::textUpdate() {
    bool activate = m_hcomp->currentIndex().isValid();
    btnSelect->setEnabled(activate);
    /*
    // TODO: Takes focus away from text field.
    if(activate) {
        btnSelect->setFocus(Qt::FocusReason::OtherFocusReason);
    } else {
        btnCancel->setFocus(Qt::FocusReason::OtherFocusReason);
    }
     */
}

void CRSSelector::keyPressEvent(QKeyEvent *e) {
    if(e->key() != Qt::Key_Escape) {
        QDialog::keyPressEvent(e);
    } else {
        reject();
    }
}

void CRSSelector::cancelClicked() {
    reject();
}

void CRSSelector::selectClicked() {
    m_hsrid = 0;
    m_vsrid = 0;
    int hidx = m_hlst.indexOf(m_hcomp->currentCompletion());
    if (hidx >= 0)
        m_hsrid = std::next(m_hcrs.begin(), hidx)->first;
    int vidx = m_vlst.indexOf(m_vcomp->currentCompletion());
    if (vidx >= 0)
        m_vsrid = std::next(m_vcrs.begin(), vidx)->first;
    g_debug("vsrid: " << hidx << ", " << m_vsrid << "; hsrid: " << vidx << ", " << m_hsrid);
    accept();
}

void CRSSelector::setVerticalSRID(int srid) {
    m_vsrid = srid;
    updateFields();
}

void CRSSelector::setHorizontalSRID(int srid) {
    m_hsrid = srid;
    updateFields();
}

int CRSSelector::getVerticalSRID() {
    return m_vsrid;
}

int CRSSelector::getHorizontalSRID() {
    return m_hsrid;
}



