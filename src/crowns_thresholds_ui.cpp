/*
 * crowns_thresholds_ui.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#include "crowns_thresholds_ui.hpp"

CrownsThresholdItem::CrownsThresholdItem(QWidget *parent) :
	QWidget(parent),
	m_index(0) {
	setFixedHeight(30);
	QGridLayout *layout = new QGridLayout();
	spnHeight = new QDoubleSpinBox();
	spnHeight->setDecimals(2);
	spnHeight->setMinimum(0);
	spnHeight->setMaximum(999);
	spnFraction = new QDoubleSpinBox();
	spnFraction->setDecimals(2);
	spnFraction->setMinimum(0.0);
	spnFraction->setMaximum(1.0);
	spnRadius = new QDoubleSpinBox();
	spnRadius->setDecimals(2);
	spnRadius->setMinimum(0.0);
	spnRadius->setMaximum(999.0);
	btnDelete = new QToolButton();
	btnDelete->setText("X");
	layout->setMargin(0);
	layout->addWidget(spnHeight, 0, 0);
	layout->addWidget(spnFraction, 0, 1);
	layout->addWidget(spnRadius, 0, 2);
	layout->addWidget(btnDelete, 0, 3);
	setLayout(layout);
	connect(btnDelete, SIGNAL(clicked()), this, SLOT(itemDeleteClicked()));
	connect(spnHeight, SIGNAL(valueChanged(double)), this, SLOT(itemHeightChanged(double)));
	connect(spnWindow, SIGNAL(valueChanged(int)), this, SLOT(itemWindowChanged(int)));
}

void CrownsThresholdItem::itemDeleteClicked() {
	emit itemDelete(this);
}

void CrownsThresholdItem::itemHeightChanged(double) {
	emit itemUpdate(this);
}

void CrownsThresholdItem::itemWindowChanged(int) {
	emit itemUpdate(this);
}

void CrownsThresholdItem::set(int index, double height, uint8_t window) {
	m_index = index;
	spnHeight->blockSignals(true);
	spnWindow->blockSignals(true);
	spnHeight->setValue(height);
	spnWindow->setValue(window);
	spnHeight->blockSignals(false);
	spnWindow->blockSignals(false);
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

bool CrownsThresholdItem::operator<(const CrownsThresholdItem &other) const {
	return height() < other.height();
}

CrownsThresholdsForm::CrownsThresholdsForm() :
	Ui::CrownsThresholdsForm(),
	m_form(nullptr),
	scrollLayout(nullptr),
	m_confirm(false) {
	this->setupUi(new QWidget());
}

CrownsThresholdsForm::~CrownsThresholdsForm() {
}

void CrownsThresholdsForm::setupUi(QWidget *form) {
	Ui::CrownsThresholdsForm::setupUi(form);
	m_form = form;
	scrollLayout = new QVBoxLayout();
	scrollLayout->addItem(new QSpacerItem(1,1, QSizePolicy::Expanding, QSizePolicy::Expanding));
	contents->setLayout(scrollLayout);
	connect(btnExit, SIGNAL(clicked()), this, SLOT(btnExitClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(btnHelpClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(btnCancelClicked()));
	connect(btnAddItem, SIGNAL(clicked()), this, SLOT(btnAddItemClicked()));
}

void CrownsThresholdsForm::setThresholds(const std::vector<std::tuple<double, double, double> > &thresholds) {
	m_thresholds = thresholds;
	sortItems();
}

void CrownsThresholdsForm::sortItems() {
	while(m_items.size() > m_thresholds.size()) {
		CrownsThresholdItem *item = m_items.back(); m_items.pop_back();
		scrollLayout->removeWidget(item);
		disconnect(item, SIGNAL(itemDelete(CrownsThresholdItem*)), this, SLOT(itemDelete(CrownsThresholdItem*)));
		disconnect(item, SIGNAL(itemUpdate(CrownsThresholdItem*)), this, SLOT(itemUpdate(CrownsThresholdItem*)));
		delete item;
	}
	while(m_items.size() < m_thresholds.size()) {
		CrownsThresholdItem *item = new CrownsThresholdItem(m_form);
		m_items.push_back(item);
		scrollLayout->insertWidget((int) m_items.size() - 1, item);
		connect(item, SIGNAL(itemDelete(CrownsThresholdItem*)), this, SLOT(itemDelete(CrownsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(CrownsThresholdItem*)), this, SLOT(itemUpdate(CrownsThresholdItem*)));
	}
	auto item = m_items.begin();
	int i = 0;
	for(auto &it : m_thresholds)
		(*item++)->set(i++, it.first, it.second);
}

std::vector<std::tuple<double, double, double> > CrownsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void CrownsThresholdsForm::btnAddItemClicked() {
	m_thresholds.push_back(std::make_tuple(0, 0, 0));
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::itemDelete(CrownsThresholdItem *item) {
	for(size_t i = item->index(); i < m_thresholds.size() - 1; ++i)
		m_thresholds[i] = m_thresholds[i + 1];
	m_thresholds.resize(m_thresholds.size() - 1);
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::itemUpdate(CrownsThresholdItem *item) {
	m_thresholds[item->index()] = std::make_tuple(item->height(), item->fraction(), item->radius());
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::updateButtons() {
	btnExit->setEnabled(valid());
}

bool CrownsThresholdsForm::valid() const {
	if(m_thresholds.empty())
		return false;
	for(const auto &it : m_thresholds) {
		if(std::get<0>(it) <= 0 || std::get<1>(it) <= 0 || std::get<2>(it) <= 0)
			return false;
	}
	return true;
}

void CrownsThresholdsForm::btnHelpClicked() {

}

void CrownsThresholdsForm::btnCancelClicked() {
	m_confirm = false;
	m_form->close();
}

void CrownsThresholdsForm::btnExitClicked() {
	m_confirm = true;
	m_form->close();
}

bool CrownsThresholdsForm::isConfirm() {
	return m_confirm;
}

