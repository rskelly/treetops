/*
 * tops_thresholds_ui.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#include "tops_thresholds_ui.hpp"

TopsThresholdItem::TopsThresholdItem(QWidget *parent) :
	QWidget(parent) {
	QGridLayout *layout = new QGridLayout();
	spnHeight = new QDoubleSpinBox();
	spnHeight->setDecimals(2);
	spnHeight->setMinimum(0);
	spnHeight->setMaximum(999);
	spnWindow = new QSpinBox();
	spnWindow->setMinimum(3);
	spnWindow->setMaximum(99);
	btnDelete = new QToolButton();
	btnDelete->setText("X");
	layout->addWidget(spnHeight);
	layout->addWidget(spnWindow);
	layout->addWidget(btnDelete);
	setLayout(layout);
}

void TopsThresholdItem::set(float height, uint8_t window) {
	spnHeight->setValue(height);
	spnWindow->setValue(window);
}

float TopsThresholdItem::height() const {
	return spnHeight->value();
}

uint8_t TopsThresholdItem::window() const {
	return spnWindow->value();
}

bool TopsThresholdItem::operator<(const TopsThresholdItem &other) const {
	return window() < other.window();
}

TopsThresholdsForm::TopsThresholdsForm() :
	Ui::TopsThresholdsForm(),
	m_form(nullptr),
	scrollLayout(nullptr) {
}

void TopsThresholdsForm::setupUi(QWidget *form) {
	Ui::TopsThresholdsForm::setupUi(form);
	m_form = form;
	scrollLayout = new QVBoxLayout();
	scrollArea->setLayout(scrollLayout);
	connect(btnAddItem, SIGNAL(clicked()), this, SLOT(btnAddItemClicked()));
}

void TopsThresholdsForm::setThresholds(const std::map<float, uint8_t> &thresholds) {
	m_thresholds = thresholds;
	for(TopsThresholdItem *item : m_items)
		scrollLayout->removeWidget(item);
	m_items.clear();
	for(const auto &it : thresholds) {
		TopsThresholdItem *item = new TopsThresholdItem(scrollArea);
		item->set(it.first, it.second);
		m_items.push_back(item);
		scrollLayout->addWidget(item);
		connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
	}
	sortItems();
}

void TopsThresholdsForm::sortItems() {

}

std::map<float, uint8_t> TopsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void TopsThresholdsForm::btnAddItemClicked() {
	TopsThresholdItem *item = new TopsThresholdItem(scrollArea);
	m_items.push_back(item);
	scrollLayout->addWidget(item);
	connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
	connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
}

void TopsThresholdsForm::itemDelete(TopsThresholdItem *item) {
	m_items.remove(item);
	scrollLayout->removeWidget(item);
}

void TopsThresholdsForm::itemUpdate(TopsThresholdItem *item) {
	m_thresholds.clear();
	for(const TopsThresholdItem &it : m_items)
		m_thresholds[it.height()] = it.window();
	sortItems();
}

