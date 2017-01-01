/*
 * ui_util.cpp
 *
 *  Created on: Dec 31, 2016
 *      Author: rob
 */

#include "ui_util.hpp"

using namespace geotools::ui::util;

// Convenience: returns a QString from a std string.
QString geotools::ui::util::qstr(const std::string &str) {
	return QString(str.c_str());
}

// Convenience: returns a QString from an int.
QString geotools::ui::util::qstr(int val) {
	QString s;
	s.setNum(val);
	return s;
}

// Hack to strip the boost::... part of an error message
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

// Convenience: open an input file dialog and return the string.
std::string geotools::ui::util::getInputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	return QFileDialog::getOpenFileName(form, qstr(title), path.path(),
							qstr("All Files (*)")).toStdString();
}

// Convenience: open an output file dialog and return the string.
std::string geotools::ui::util::getOutputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	return QFileDialog::getSaveFileName(form, qstr(title), path.path(),
			qstr("All Files (*)")).toStdString();
}


