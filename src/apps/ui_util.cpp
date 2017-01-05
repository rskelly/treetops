/*
 * ui_util.cpp
 *
 *  Created on: Dec 31, 2016
 *      Author: rob
 */

#include "ui_util.hpp"
#include "tops_thresholds_ui.hpp"
#include "geotools.hpp"

using namespace geotools::ui::util;

QString geotools::ui::util::qstr(const std::string &str) {
	return QString(str.c_str());
}

QString geotools::ui::util::qstr(int val) {
	QString s;
	s.setNum(val);
	return s;
}

std::string geotools::ui::util::stripBoost(const std::string &msg) {
	if (msg.substr(0, 7) == "boost::")
		return msg.substr(msg.find(" ", 0));
	return msg;
}

void geotools::ui::util::errorDialog(QWidget *parent, const std::string &title, const std::string &text,
	const std::string &detail) {
	QMessageBox err(parent);
	err.setWindowTitle(qstr(title));
	err.setText(qstr(text));
	if(!detail.empty())
		err.setDetailedText(qstr(detail));
	err.exec();
}

std::string geotools::ui::util::getInputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	QString res = QFileDialog::getOpenFileName(form, qstr(title), path.path(),
							qstr(filter));
	return res.toStdString();
}

std::string geotools::ui::util::getOutputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	QString res = QFileDialog::getSaveFileName(form, qstr(title), path.path(),
			qstr(filter));
	return res.toStdString();
}

void geotools::ui::util::getThresholds(QWidget *form, std::map<float, uint8_t> &thresholds) {
	TopsThresholdsForm *tf = new TopsThresholdsForm();
	tf->setThresholds(thresholds);
	QDialog *dlg = new QDialog();
	tf->setupUi(dlg);
	if(dlg->exec() && tf->isConfirm())
		thresholds = tf->thresholds();
}


