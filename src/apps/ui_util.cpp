/*
 * ui_util.cpp
 *
 *  Created on: Dec 31, 2016
 *      Author: rob
 */

#include "ui_util.hpp"

using namespace geotools::ui::util;

// Convenience: returns a QString from a std string.
QString qstr(const std::string &str) {
	return QString(str.c_str());
}

// Convenience: returns a QString from an int.
QString qstr(int val) {
	QString s;
	s.setNum(val);
	return s;
}

// Hack to strip the boost::... part of an error message
std::string stripBoost(const std::string &msg) {
	if (msg.substr(0, 7) == "boost::")
		return msg.substr(msg.find(" ", 0));
	return msg;
}

void errorDialog(QWidget *parent, const std::string &title, const std::string &text,
	const std::string &detail = "") {
	QMessageBox err(parent);
	err.setWindowTitle(qstr(title));
	err.setText(qstr(text));
	if(!detail.empty())
		err.setDetailedText(qstr(detail));
	err.exec();
}

// Convenience: open an input file dialog and return the string.
std::string getInputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	/*
	QFileDialog *d = new QFileDialog(form, qstr(title), path.path(), qstr("All Files (*)"));
	d->setAcceptMode(QFileDialog::AcceptOpen);
	d->setFileMode(QFileDialog::ExistingFile);
	d->setViewMode(QFileDialog::Detail);
	d->setDirectory(path);
	//d->setOption( QFileDialog::DontUseNativeDialog, true);
	QString res;
	if(d->exec()) {
		res = d->selectedFiles().first();
		path.setPath(res);
	}
	return res.toStdString();
	*/
	return QFileDialog::getOpenFileName(form, qstr(title), path.path(),
							qstr("All Files (*)")).toStdString();
}

// Convenience: open an output file dialog and return the string.
std::string getOutputFile(QWidget *form, const std::string &title, QDir &path,
		const std::string &filter) {
	/*
	QFileDialog *d = new QFileDialog(form, qstr(title), path.path(), qstr("All Files (*)"));
	d->setAcceptMode(QFileDialog::AcceptSave);
	d->setFileMode(QFileDialog::AnyFile);
	d->setViewMode(QFileDialog::Detail);
	d->setDirectory(path);
	//d->setOption( QFileDialog::DontUseNativeDialog, true);
	QString res;
	if(d->exec()) {
		res = d->selectedFiles().first();
		path.setPath(res);
	}
	return res.toStdString();
	}
	*/
	return QFileDialog::getSaveFileName(form, qstr(title), path.path(),
			qstr("All Files (*)")).toStdString();
}


