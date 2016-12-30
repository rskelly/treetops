/*
 * ui_util.hpp
 *
 *  Created on: Dec 18, 2016
 *      Author: rob
 */

#ifndef __UI_UTIL_HPP__
#define __UI_UTIL_HPP__

#include <string>

#include <QString>
#include <QWidget>
#include <QFileDialog>
#include <QDir>

namespace geotools {
	namespace ui {
		namespace util {

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

			// Convenience: open an input file dialog and return the string.
			std::string getInputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter) {
				QFileDialog *d = new QFileDialog(form, qstr(title), path.path());//, qstr(filter));
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
			}

			// Convenience: open an output file dialog and return the string.
			std::string getOutputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter) {
				QFileDialog *d = new QFileDialog(form, qstr(title), path.path());//, qstr(filter));
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

		}
	}
}

#endif /* __UI_UTIL_HPP__ */
