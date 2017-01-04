/*
 * tops_thresholds_ui.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#ifndef __TOPS_THRESHOLDS_UI_HPP__
#define __TOPS_THRESHOLDS_UI_HPP__

#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QToolButton>
#include <QWidget>

#include "ui_tops_thresholds.h"

class TopsThresholdItem : public QWidget {
	Q_OBJECT
public:
	QDoubleSpinBox *spnHeight;
	QSpinBox *spnWindow;
	QToolButton *btnDelete;

	TopsThresholdItem(QWidget *parent = 0);
	void set(float height, uint8_t window);
	float height() const;
	uint8_t window() const;
	bool operator<(const TopsThresholdItem &other) const;

signals:
	void itemUpdate(TopsThresholdItem *item);
	void itemDelete(TopsThresholdItem *item);
};

class TopsThresholdsForm : public QObject, public Ui::TopsThresholdsForm {
	Q_OBJECT
private:
	void sortItems();

public:
	QWidget *m_form;
	QVBoxLayout *scrollLayout;
	std::list<TopsThresholdItem*> m_items;
	std::map<float, uint8_t> m_thresholds;

	TopsThresholdsForm();
	void setThresholds(const std::map<float, uint8_t> &thresholds);
	std::map<float, uint8_t> thresholds() const;
	void setupUi(QWidget *form);

public slots:
	void btnAddItemClicked();
	void itemDelete(TopsThresholdItem *item);
	void itemUpdate(TopsThresholdItem *item);
};

#endif /* __TOPS_THRESHOLDS_UI_HPP__ */
