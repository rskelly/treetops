/*
 * tops_thresholds_ui.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#include "ui/tops_threshold_item.hpp"

using namespace tt::ui;

TopsThresholdItem::TopsThresholdItem(QWidget* parent) : 
		QWidget(parent),
		m_index(0) {
	setFixedHeight(30);
	QGridLayout* layout = new QGridLayout();
	spnHeight = new QDoubleSpinBox();
	spnHeight->setDecimals(2);
	spnHeight->setMinimum(0);
	spnHeight->setMaximum(999);
	spnWindow = new QSpinBox();
	spnWindow->setMinimum(3);
	spnWindow->setMaximum(99);
	spnWindow->setSingleStep(2);
	btnDelete = new QToolButton();
	btnDelete->setText("X");
	layout->setContentsMargins(0, 0, 0, 0);
	layout->addWidget(spnHeight, 0, 0);
	layout->addWidget(spnWindow, 0, 1);
	layout->addWidget(btnDelete, 0, 2);
	setLayout(layout);
	connect(btnDelete, SIGNAL(clicked()), this, SLOT(itemDeleteClicked()));
	connect(spnHeight, SIGNAL(valueChanged(double)), this, SLOT(itemHeightChanged(double)));
	connect(spnWindow, SIGNAL(valueChanged(int)), this, SLOT(itemWindowChanged(int)));
}

void TopsThresholdItem::itemDeleteClicked() {
	emit itemDelete(this);
}

void TopsThresholdItem::itemHeightChanged(double) {
	emit itemUpdate(this);
}

void TopsThresholdItem::itemWindowChanged(int) {
	emit itemUpdate(this);
}

void TopsThresholdItem::set(int index, double height, int window) {
	m_index = index;
	spnHeight->blockSignals(true);
	spnWindow->blockSignals(true);
	spnHeight->setValue(height);
	spnWindow->setValue(window);
	spnHeight->blockSignals(false);
	spnWindow->blockSignals(false);
}

int TopsThresholdItem::index() const {
	return m_index;
}

double TopsThresholdItem::height() const {
	return spnHeight->value();
}

int TopsThresholdItem::window() const {
	return spnWindow->value();
}

bool TopsThresholdItem::operator<(const TopsThresholdItem& other) const {
	return height() < other.height();
}

/*
void TopsThresholdItem::itemDelete(TopsThresholdItem* item) {

}

void TopsThresholdItem::itemUpdate(TopsThresholdItem* item) {

}
*/