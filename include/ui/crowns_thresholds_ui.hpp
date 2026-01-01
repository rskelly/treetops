/*
 * crowns_thresholds_ui.hpp
 *
 *  Created on: Feb 20, 2017
 *      Author: rob
 */

#ifndef __CROWNS_THRESHOLDS_UI_HPP__
#define __CROWNS_THRESHOLDS_UI_HPP__

#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>

#include "treetops.hpp"

#include "ui/crowns_threshold_item.hpp"
#include "ui/ui_crowns_thresholds.h"

namespace tt {
namespace ui {

// Represents the thresholds dialog.
class CrownsThresholdsForm : public QObject, public Ui::CrownsThresholdsForm {
	Q_OBJECT
private:
	void sortItems();
	bool valid() const;
	void updateButtons();

public:
	QWidget* m_form;
	QVBoxLayout* scrollLayout;
	std::list<CrownsThresholdItem*> m_items;
	std::vector<tt::config::CrownThreshold> m_thresholds;
	bool m_confirm;

	CrownsThresholdsForm();
	void setThresholds(const std::vector<tt::config::CrownThreshold> &thresholds);
	std::vector<tt::config::CrownThreshold> thresholds() const;
	void setupUi(QWidget *form);
	bool isConfirm();
	~CrownsThresholdsForm();

public slots:
	void btnAddItemClicked();
	void btnHelpClicked();
	void btnCancelClicked();
	void btnExitClicked();

signals:
	void itemDelete(CrownsThresholdItem* item);
	void itemUpdate(CrownsThresholdItem* item);
};

} // ui
} // tt

#endif /* __CROWNS_THRESHOLDS_UI_HPP__ */
