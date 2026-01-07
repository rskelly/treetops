/*
 * crowns_thresholds_ui.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#include <iostream>

#include "treetops.hpp"
#include "settings.hpp"

#include "ui/ui_crowns_thresholds.h"
#include "ui/crowns_thresholds_ui.hpp"


using namespace tt::ui;
using namespace tt::config;


/*
CrownsThresholdsForm::CrownsThresholdsForm() :
	Ui::CrownsThresholdsForm(),
	m_form(nullptr),
	scrollLayout(nullptr),
	m_confirm(false) {
	this->setupUi(new QWidget());
}
*/

void CrownsThresholdsForm::setupUi(QWidget* form) {
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

void CrownsThresholdsForm::setThresholds(const std::vector<CrownThreshold>& thresholds) {
	m_thresholds = thresholds;
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::sortItems() {
	while(m_items.size() > m_thresholds.size()) {
		CrownsThresholdItem* item = m_items.back(); m_items.pop_back();
		scrollLayout->removeWidget(item);
		disconnect(item, SIGNAL(itemDelete(CrownsThresholdItem*)), this, SLOT(itemDelete(CrownsThresholdItem*)));
		disconnect(item, SIGNAL(itemUpdate(CrownsThresholdItem*)), this, SLOT(itemUpdate(CrownsThresholdItem*)));
		delete item;
	}
	while(m_items.size() < m_thresholds.size()) {
		CrownsThresholdItem* item = new CrownsThresholdItem(m_form);
		m_items.push_back(item);
		scrollLayout->insertWidget((int) m_items.size() - 1, item);
		connect(item, SIGNAL(itemDelete(CrownsThresholdItem*)), this, SLOT(itemDelete(CrownsThresholdItem*)));
		connect(item, SIGNAL(itemUpdate(CrownsThresholdItem*)), this, SLOT(itemUpdate(CrownsThresholdItem*)));
	}
	auto item = m_items.begin();
	int i = 0;
	for(const CrownThreshold& c : m_thresholds)
		(*item++)->set(i++, c.threshold, c.fraction, c.radius);
}

std::vector<CrownThreshold> CrownsThresholdsForm::thresholds() const {
	return m_thresholds;
}

void CrownsThresholdsForm::btnAddItemClicked() {
	m_thresholds.emplace_back(0, 0, 0);
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::itemDelete(CrownsThresholdItem* item) {
	for(size_t i = item->index(); i < m_thresholds.size() - 1; ++i)
		m_thresholds[i] = m_thresholds[i + 1];
	m_thresholds.resize(m_thresholds.size() - 1);
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::itemUpdate(CrownsThresholdItem *item) {
	CrownThreshold& c = m_thresholds[item->index()];
	c.threshold = item->height();
	c.fraction = item->fraction();
	c.radius = item->radius();
	sortItems();
	updateButtons();
}

void CrownsThresholdsForm::updateButtons() {
	btnExit->setEnabled(valid());
}

bool CrownsThresholdsForm::valid() const {
	if(m_thresholds.empty())
		return false;
	for(const CrownThreshold& c : m_thresholds) {
		if(c.fraction <= 0 || c.radius <= 0 || c.threshold <= 0)
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

CrownsThresholdsForm::~CrownsThresholdsForm() = default;