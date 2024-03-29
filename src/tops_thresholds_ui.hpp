/*
 * tops_thresholds_ui.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#ifndef __TOPS_THRESHOLDS_UI_HPP__
#define __TOPS_THRESHOLDS_UI_HPP__

#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QWidget>

#include "treetops.hpp"
#include "ui_tops_thresholds.h"

// Represents a single line in the tops thresholds dialog.
class TopsThresholdItem : public QWidget {
	Q_OBJECT
private:
	int m_index;

public:
	QDoubleSpinBox *spnHeight;
	QSpinBox *spnWindow;
	QToolButton *btnDelete;
	TopsThresholdItem(QWidget *parent = 0);
	void set(int index, double height, int window);
	double height() const;
	int window() const;
	bool operator<(const TopsThresholdItem &other) const;
	int index() const;

public slots:
	void itemDeleteClicked();
	void itemWindowChanged(int);
	void itemHeightChanged(double);

signals:
	void itemUpdate(TopsThresholdItem *item);
	void itemDelete(TopsThresholdItem *item);
};

// Represents the thresholds dialog.
class TopsThresholdsForm : public QObject, public Ui::TopsThresholdsForm {
	Q_OBJECT
private:
	void sortItems();
	bool valid() const;
	void updateButtons();

public:
	QWidget *m_form;
	QVBoxLayout *scrollLayout;
	std::list<TopsThresholdItem*> m_items;
	std::vector<geo::treetops::config::TopThreshold> m_thresholds;
	bool m_confirm;

	TopsThresholdsForm();
	void setThresholds(const std::vector<geo::treetops::config::TopThreshold> &thresholds);
	std::vector<geo::treetops::config::TopThreshold> thresholds() const;
	void setupUi(QWidget *form);
	bool isConfirm();
	~TopsThresholdsForm();

public slots:
	void btnAddItemClicked();
	void btnHelpClicked();
	void btnCancelClicked();
	void btnExitClicked();
	void itemDelete(TopsThresholdItem *item);
	void itemUpdate(TopsThresholdItem *item);
};

#endif /* __TOPS_THRESHOLDS_UI_HPP__ */
