/*
 * tops_thresholds_ui.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#include "tops_thresholds_ui.hpp"

TopsThresholdItem::TopsThresholdItem(QWidget *parent) :
	QWidget(parent) {
	setFixedHeight(30);
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
	layout->setMargin(0);
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
	scrollLayout(nullptr),
	m_confirm(false) {
	this->setupUi(new QWidget());
}

TopsThresholdsForm::~TopsThresholdsForm() {
	delete m_form;
}

void TopsThresholdsForm::setupUi(QWidget *form) {
	Ui::TopsThresholdsForm::setupUi(form);
	m_form = form;
	scrollLayout = new QVBoxLayout();
	scrollLayout->addItem(new QSpacerItem(1,1, QSizePolicy::Expanding, QSizePolicy::Expanding));
	contents->setLayout(scrollLayout);
	connect(btnExit, SIGNAL(clicked()), this, SLOT(btnExitClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(btnHelpClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(btnCancelClicked()));
	connect(btnAddItem, SIGNAL(clicked()), this, SLOT(btnAddItemClicked()));
}

void TopsThresholdsForm::setThresholds(const std::map<float, uint8_t> &thresholds) {
	m_thresholds = thresholds;
	sortItems();
}

void TopsThresholdsForm::sortItems() {
	for(TopsThresholdItem *item : m_items) {
		scrollLayout->removeWidget(item);
		disconnect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		disconnect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
		delete item;
	}
	m_items.clear();
	int idx = 0;
	for(const auto &it : m_thresholds) {
		TopsThresholdItem *item = new TopsThresholdItem(m_form);
		item->set(it.first, it.second);
		m_items.push_back(item);
		scrollLayout->insertWidget(idx++, item);
		connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
	}
}

std::map<float, uint8_t> TopsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void TopsThresholdsForm::btnAddItemClicked() {
	TopsThresholdItem *item = new TopsThresholdItem(m_form);
	m_items.push_back(item);
	scrollLayout->insertWidget(0, item);
	connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
	connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
}

void TopsThresholdsForm::itemDelete(TopsThresholdItem *item) {
	m_items.remove(item);
	disconnect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
	disconnect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
	scrollLayout->removeWidget(item);
	delete item;
}

void TopsThresholdsForm::itemUpdate(TopsThresholdItem *item) {
	m_thresholds.clear();
	for(const TopsThresholdItem &it : m_items)
		m_thresholds[it.height()] = it.window();
	sortItems();
	btnExit->setEnabled(true);
}

void TopsThresholdsForm::btnHelpClicked() {

}

void TopsThresholdsForm::btnCancelClicked() {
	m_confirm = false;
	m_form->close();
}

void TopsThresholdsForm::btnExitClicked() {
	m_confirm = true;
	m_form->close();
}

bool TopsThresholdsForm::isConfirm() {
	return m_confirm;
}

