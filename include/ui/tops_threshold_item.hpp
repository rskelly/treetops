/*
 * tops_thresholds_ui.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: rob
 */

#ifndef __TOPS_THRESHOLD_ITEM_HPP__
#define __TOPS_THRESHOLD_ITEM_HPP__

#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>


namespace tt {
namespace ui {

// Represents a single line in the tops thresholds dialog.
class TopsThresholdItem : public QWidget {
	Q_OBJECT
private:
	int m_index;
	QDoubleSpinBox* spnHeight;
	QSpinBox* spnWindow;
	QToolButton* btnDelete;

public:
	explicit TopsThresholdItem(QWidget* parent = nullptr);
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
	void itemUpdate(TopsThresholdItem* item);
	void itemDelete(TopsThresholdItem* item);
};

} // ui
} // tt

#endif /* __TOPS_THRESHOLD_ITEM_HPP__ */
