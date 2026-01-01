/*
 * ui_util.cpp
 *
 *  Created on: Dec 31, 2016
 *      Author: rob
 */

#ifndef _UI_UTIL_HPP_
#define _UI_UTIL_HPP_

#include "ui/ui_util.hpp"
#include "ui/crowns_thresholds_ui.hpp"
#include "ui/tops_thresholds_ui.hpp"

using namespace tt::ui;
using namespace tt::ui::util;
using namespace tt::config;

QString tt::ui::util::qstr(const std::string& str) {
	return QString(str.c_str());
}

QString tt::ui::util::qstr(int val) {
	QString s;
	s.setNum(val);
	return s;
}

std::string tt::ui::util::sstr(const QString& str) {
	std::vector<wchar_t> buf(str.size() + 1);
	str.toWCharArray(buf.data());
	std::wstring ws(buf.data());
	std::string s(ws.begin(), ws.end());
	return s;
}

std::string tt::ui::util::stripBoost(const std::string& msg) {
	if (msg.substr(0, 7) == "boost::")
		return msg.substr(msg.find(" ", 0));
	return msg;
}

void tt::ui::util::errorDialog(QWidget* parent, const std::string& title, const std::string& text,
	const std::string& detail) {
	QMessageBox err(parent);
	err.setWindowTitle(qstr(title));
	err.setText(qstr(text));
	if(!detail.empty())
		err.setDetailedText(qstr(detail));
	err.exec();
}

void tt::ui::util::infoDialog(QWidget* parent, const std::string& title, const std::string& text, const std::string& /*detail*/) {
	QMessageBox::information(parent, qstr(title), qstr(text), QMessageBox::Ok);
}

void tt::ui::util::getInputFile(QWidget* form, const std::string& title, std::string& path,
		const std::string& filter, std::string& filename) {
	QString res = QFileDialog::getOpenFileName(form, qstr(title), qstr(path), qstr(filter));
	if(!res.isEmpty()) {
		path = sstr(res);
		filename = sstr(res);
	}
}

void tt::ui::util::getOutputFile(QWidget* form, const std::string& title, std::string& path,
		const std::string& filter, std::string& filename, bool confirmOverwrite) {
	QFileDialog::Option opt = (QFileDialog::Option) 0;
	if(!confirmOverwrite)
		opt = QFileDialog::DontConfirmOverwrite;
	QString res = QFileDialog::getSaveFileName(form, qstr(title), QString(path.c_str()), qstr(filter), nullptr, opt);
	if(!res.isEmpty()) {
		path = res.toStdString();
		filename = res.toStdString();
	}
}

void tt::ui::util::getTopsThresholds(QWidget* /*form*/, std::vector<TopThreshold>& thresholds) {
	/*
	tt::ui::TopsThresholdsForm tf;
	QDialog dlg;
	tf.setupUi(&dlg);
	tf.setThresholds(thresholds);
	dlg.exec(); // TODO: There is a correct way to manage accepted response from dialogs.
	if(tf.isConfirm())
		thresholds = tf.thresholds();
	*/
}

void tt::ui::util::getCrownsThresholds(QWidget* /*form*/, std::vector<CrownThreshold>& thresholds) {
	/*
	tt::ui::CrownsThresholdsForm tf;
	QDialog dlg;
	tf.setupUi(&dlg);
	tf.setThresholds(thresholds);
	dlg.exec(); // TODO: There is a correct way to manage accepted response from dialogs.
	if(tf.isConfirm())
		thresholds = tf.thresholds();
	*/
}


#endif // _UI_UTIL_HPP_
