/*
 * tops_thresholds_ui.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#include "treetops.hpp"
#include "tops_thresholds_ui.hpp"

using namespace geo::treetops::config;

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

bool TopsThresholdItem::operator<(const TopsThresholdItem &other) const {
	return height() < other.height();
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

void TopsThresholdsForm::setThresholds(const std::vector<TopThreshold> &thresholds) {
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
		scrollLayout->insertWidget((int) m_items.size() - 1, item);
		connect(item, SIGNAL(itemDelete(TopsThresholdItem*)), this, SLOT(itemDelete(TopsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(TopsThresholdItem*)), this, SLOT(itemUpdate(TopsThresholdItem*)));
	}
	auto item = m_items.begin();
	int i = 0;
	for(const TopThreshold& t : m_thresholds)
		(*item++)->set(i++, t.threshold, t.window);
}

std::vector<TopThreshold> TopsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void TopsThresholdsForm::btnAddItemClicked() {
	m_thresholds.emplace_back(0, 0);
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
	TopThreshold& t = m_thresholds[item->index()];
	t.threshold = item->height();
	t.window = item->window();
	sortItems();
	updateButtons();
}

void TopsThresholdsForm::updateButtons() {
	btnExit->setEnabled(valid());
}

bool TopsThresholdsForm::valid() const {
	if(m_thresholds.empty())
		return false;
	for(const TopThreshold& t : m_thresholds) {
		if(t.threshold <= 0 || t.window % 2 == 0 || t.window < 3)
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

