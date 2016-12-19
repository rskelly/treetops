/*
 * ui_util.hpp
 *
 *  Created on: Dec 18, 2016
 *      Author: rob
 */

#ifndef SRC_APPS_UI_UTIL_CPP_
#define SRC_APPS_UI_UTIL_CPP_

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

			// Convenience: open an input file dialog and return the string.
			std::string getInputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter) {
				QString res = QFileDialog::getOpenFileName(form, qstr(title), path.path(), qstr(filter));
				path.setPath(res);
				return res.toStdString();
			}

			// Convenience: open an output file dialog and return the string.
			std::string getOutputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter) {
				QString res = QFileDialog::getSaveFileName(form, qstr(title), path.path(), qstr(filter));
				path.setPath(res);
				return res.toStdString();
			}

		}
	}
}

#endif /* SRC_APPS_UI_UTIL_CPP_ */
