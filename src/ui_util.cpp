/*
 * ui_util.cpp
 *
 *  Created on: Dec 31, 2016
 *      Author: rob
 */

#include "ui_util.hpp"
#include "tops_thresholds_ui.hpp"
#include "crowns_thresholds_ui.hpp"
#include "geo.hpp"

using namespace geo::ui::util;

QString geo::ui::util::qstr(const std::string &str) {
	return QString(str.c_str());
}

QString geo::ui::util::qstr(int val) {
	QString s;
	s.setNum(val);
	return s;
}

std::string geo::ui::util::stripBoost(const std::string &msg) {
	if (msg.substr(0, 7) == "boost::")
		return msg.substr(msg.find(" ", 0));
	return msg;
}

void geo::ui::util::errorDialog(QWidget *parent, const std::string &title, const std::string &text,
	const std::string &detail) {
	QMessageBox err(parent);
	err.setWindowTitle(qstr(title));
	err.setText(qstr(text));
	if(!detail.empty())
		err.setDetailedText(qstr(detail));
	err.exec();
}

void geo::ui::util::getInputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter, std::string &filename) {
	QString res = QFileDialog::getOpenFileName(form, qstr(title), path.path(), qstr(filter));
	if(!res.isEmpty()) {
		path.setPath(res);
		filename = res.toStdString();
	}
}

void geo::ui::util::getOutputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter, std::string &filename) {
	QString res = QFileDialog::getSaveFileName(form, qstr(title), path.path(), qstr(filter));
	if(!res.isEmpty()) {
		path.setPath(res);
		filename = res.toStdString();
	}
}

void geo::ui::util::getTopsThresholds(QWidget *form, std::vector<std::tuple<double, uint8_t> > &thresholds) {
	TopsThresholdsForm tf;
	QDialog dlg;
	tf.setupUi(&dlg);
	tf.setThresholds(thresholds);
	dlg.exec(); // TODO: There is a correct way to manage accepted response from dialogs.
	if(tf.isConfirm())
		thresholds = tf.thresholds();
}

void geo::ui::util::getCrownsThresholds(QWidget *form, std::vector<std::tuple<double, double, double> > &thresholds) {
	CrownsThresholdsForm tf;
	QDialog dlg;
	tf.setupUi(&dlg);
	tf.setThresholds(thresholds);
	dlg.exec(); // TODO: There is a correct way to manage accepted response from dialogs.
	if(tf.isConfirm())
		thresholds = tf.thresholds();
}


