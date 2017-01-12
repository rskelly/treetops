/*
 * tops_thresholds_ui.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#include "tops_thresholds_ui.hpp"

TopsThresholdItem::TopsThresholdItem(QWidget *parent) :
	QWidget(parent),
	m_index(0) {
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

void TopsThresholdItem::set(int index, float height, uint8_t window) {
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

void TopsThresholdsForm::setThresholds(const std::vector<std::pair<float, uint8_t> > &thresholds) {
	m_thresholds = thresholds;
	sortItems();
}

void TopsThresholdsForm::sortItems() {
	while(m_items.size() > m_thresholds.size()) {
		TopsThresholdItem *item = m_items.back(); m_items.pop_back();
		scrollLayout->removeWidget(item);
		disconnect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		disconnect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
		delete item;
	}
	while(m_items.size() < m_thresholds.size()) {
		TopsThresholdItem *item = new TopsThresholdItem(m_form);
		m_items.push_back(item);
		scrollLayout->insertWidget(m_items.size() - 1, item);
		connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
	}
	auto item = m_items.begin();
	int i = 0;
	for(auto &it : m_thresholds)
		(*item++)->set(i++, it.first, it.second);
}

std::vector<std::pair<float, uint8_t> > TopsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void TopsThresholdsForm::btnAddItemClicked() {
	m_thresholds.push_back(std::make_pair(0, 0));
	sortItems();
	updateButtons();
}

void TopsThresholdsForm::itemDelete(TopsThresholdItem *item) {
	for(size_t i = item->index(); i < m_thresholds.size() - 1; ++i)
		m_thresholds[i] = m_thresholds[i + 1];
	m_thresholds.resize(m_thresholds.size() - 1);
	sortItems();
	updateButtons();
}

void TopsThresholdsForm::itemUpdate(TopsThresholdItem *item) {
	m_thresholds[item->index()] = std::make_pair(item->height(), item->window());
	sortItems();
	updateButtons();
}

void TopsThresholdsForm::updateButtons() {
	btnExit->setEnabled(valid());
}

bool TopsThresholdsForm::valid() const {
	if(m_thresholds.empty())
		return false;
	for(const auto &it : m_thresholds) {
		if(it.first <= 0 || it.second % 2 == 0 || it.second < 3)
			return false;
	}
	return true;
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

