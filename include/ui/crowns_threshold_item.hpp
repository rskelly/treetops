/*
 * crowns_thresholds_ui.hpp
 * 
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#ifndef __CROWNS_THRESHOLD_ITEM_HPP__
#define __CROWNS_THRESHOLD_ITEM_HPP__

#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>

#include "treetops.hpp"

#include "ui/ui_crowns_thresholds.h"

namespace tt {
namespace ui {

// Represents a single line in the crowns thresholds dialog.
class CrownsThresholdItem : public QWidget, Ui::CrownsThresholdsForm {
	Q_OBJECT
private:
	int m_index;

public:
	QDoubleSpinBox* spnHeight;
	QDoubleSpinBox* spnFraction;
	QDoubleSpinBox* spnRadius;
	QToolButton* btnDelete;

	CrownsThresholdItem(QWidget* parent = 0);
	void set(int index, double height, double fraction, double radius);
	double height() const;
	double fraction() const;
	double radius() const;
	bool operator<(const CrownsThresholdItem& other) const;
	int index() const;

public slots:
	void itemDeleteClicked();
	void itemHeightChanged(double);
	void itemFractionChanged(double);
	void itemRadiusChanged(double);

signals:
	void itemUpdate(CrownsThresholdItem* item);
	void itemDelete(CrownsThresholdItem* item);
};


} // ui
} // tt

#endif /* __CROWNS_THRESHOLD_ITEM_HPP__ */
