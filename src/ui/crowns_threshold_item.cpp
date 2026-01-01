/*
 * crowns_thresholds_ui.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#include <iostream>

#include "treetops.hpp"
#include "settings.hpp"

#include "ui/crowns_threshold_item.hpp"

using namespace tt::ui;
using namespace tt::config;

CrownsThresholdItem::CrownsThresholdItem(QWidget* parent) :
	QWidget(parent),
	m_index(0) {
	setFixedHeight(30);
	QGridLayout* layout = new QGridLayout();
	spnHeight = new QDoubleSpinBox();
	spnHeight->setDecimals(2);
	spnHeight->setMinimum(0);
	spnHeight->setMaximum(999);
	spnFraction = new QDoubleSpinBox();
	spnFraction->setDecimals(2);
	spnFraction->setMinimum(0.0);
	spnFraction->setMaximum(1.0);
	spnFraction->setSingleStep(0.1);
	spnRadius = new QDoubleSpinBox();
	spnRadius->setDecimals(2);
	spnRadius->setMinimum(0.0);
	spnRadius->setMaximum(999.0);
	btnDelete = new QToolButton();
	btnDelete->setText("X");
	layout->setContentsMargins(0, 0, 0, 0);
	layout->addWidget(spnHeight, 0, 0);
	layout->addWidget(spnFraction, 0, 1);
	layout->addWidget(spnRadius, 0, 2);
	layout->addWidget(btnDelete, 0, 3);
	setLayout(layout);
	connect(btnDelete, SIGNAL(clicked()), this, SLOT(itemDeleteClicked()));
	connect(spnHeight, SIGNAL(valueChanged(double)), this, SLOT(itemHeightChanged(double)));
	connect(spnRadius, SIGNAL(valueChanged(double)), this, SLOT(itemRadiusChanged(double)));
	connect(spnFraction, SIGNAL(valueChanged(double)), this, SLOT(itemFractionChanged(double)));
}

void CrownsThresholdItem::itemDeleteClicked() {
	emit itemDelete(this);
}

void CrownsThresholdItem::itemHeightChanged(double) {
	emit itemUpdate(this);
}

void CrownsThresholdItem::itemFractionChanged(double) {
	emit itemUpdate(this);
}

void CrownsThresholdItem::itemRadiusChanged(double) {
	emit itemUpdate(this);
}

void CrownsThresholdItem::set(int index, double height, double fraction, double radius) {
	m_index = index;
	spnHeight->blockSignals(true);
	spnFraction->blockSignals(true);
	spnRadius->blockSignals(true);
	spnHeight->setValue(height);
	spnFraction->setValue(fraction);
	spnRadius->setValue(radius);
	spnHeight->blockSignals(false);
	spnFraction->blockSignals(false);
	spnRadius->blockSignals(false);
}

int CrownsThresholdItem::index() const {
	return m_index;
}

double CrownsThresholdItem::height() const {
	return spnHeight->value();
}

double CrownsThresholdItem::fraction() const {
	return spnFraction->value();
}

double CrownsThresholdItem::radius() const {
	return spnRadius->value();
}

bool CrownsThresholdItem::operator<(const CrownsThresholdItem& other) const {
	return height() < other.height();
}

/*
void CrownsThresholdItem::itemUpdate(CrownsThresholdItem* item) {}
void CrownsThresholdItem::itemDelete(CrownsThresholdItem* item) {}
*/
