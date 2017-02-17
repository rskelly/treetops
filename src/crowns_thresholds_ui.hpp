/*
 * crowns_thresholds_ui.hpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#ifndef __CROWNS_THRESHOLDS_UI_HPP__
#define __CROWNS_THRESHOLDS_UI_HPP__

#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QToolButton>
#include <QWidget>

#include "ui_crowns_thresholds.h"

// Represents a single line in the crowns thresholds dialog.
class CrownsThresholdItem : public QWidget {
	Q_OBJECT
private:
	int m_index;

public:
	QDoubleSpinBox *spnHeight;
	QDoubleSpinBox *spnFraction;
	QDoubleSpinBox *spnRadius;
	QToolButton *btnDelete;
	CrownsThresholdItem(QWidget *parent = 0);
	void set(int index, double height, double fraction, double radius);
	double height() const;
	double fraction() const;
	double radius() const;
	bool operator<(const CrownsThresholdItem &other) const;
	int index() const;

public slots:
	void itemDeleteClicked();
	void itemHeightChanged(double);
	void itemFractionChanged(double);
	void itemRadiusChanged(double);

signals:
	void itemUpdate(CrownsThresholdItem *item);
	void itemDelete(CrownsThresholdItem *item);
};

// Represents the thresholds dialog.
class CrownsThresholdsForm : public QObject, public Ui::CrownsThresholdsForm {
	Q_OBJECT
private:
	void sortItems();
	bool valid() const;
	void updateButtons();

public:
	QWidget *m_form;
	QVBoxLayout *scrollLayout;
	std::list<CrownsThresholdItem*> m_items;
	std::vector<std::tuple<double, double, double> > m_thresholds;
	bool m_confirm;

	CrownsThresholdsForm();
	void setThresholds(const std::vector<std::tuple<double, double, double> > &thresholds);
	std::vector<std::tuple<double, double, double> > thresholds() const;
	void setupUi(QWidget *form);
	bool isConfirm();
	~CrownsThresholdsForm();

public slots:
	void btnAddItemClicked();
	void btnHelpClicked();
	void btnCancelClicked();
	void btnExitClicked();
	void itemDelete(CrownsThresholdItem *item);
	void itemUpdate(CrownsThresholdItem *item);
};

#endif /* __CROWNS_THRESHOLDS_UI_HPP__ */
