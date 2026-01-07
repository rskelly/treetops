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
#include <QtWidgets/QVBoxLayout>

#include "treetops.hpp"

#include "ui/ui_tops_thresholds.h"
#include "ui/tops_threshold_item.hpp"

namespace tt {
namespace ui {

// Represents the thresholds dialog.
class TopsThresholdsForm : public QWidget, Ui::TopsThresholdsForm {
	Q_OBJECT
private:
	void sortItems();
	bool valid() const;
	void updateButtons();
	
public:
	QWidget* m_form;
	QVBoxLayout* scrollLayout;
	std::list<TopsThresholdItem*> m_items;
	std::vector<tt::config::TopThreshold> m_thresholds;
	bool m_confirm;

	// TopsThresholdsForm(QWidget* parent = nullptr);
	void setThresholds(const std::vector<tt::config::TopThreshold>&);
	std::vector<tt::config::TopThreshold> thresholds() const;
	void setupUi(QWidget* form);
	bool isConfirm();

	~TopsThresholdsForm();

//signals:

public slots:
	void itemDelete(TopsThresholdItem*);
	void itemUpdate(TopsThresholdItem*);
	void btnAddItemClicked();
	void btnHelpClicked();
	void btnCancelClicked();
	void btnExitClicked();

};

} // ui
} // tt

#endif /* __TOPS_THRESHOLDS_UI_HPP__ */
