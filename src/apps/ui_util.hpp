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
#include <QMessageBox>
#include <QDir>

#define VECTOR_PATTERN "Vector Files (*.sqlite *.shp)"
#define RASTER_PATTERN "Raster Files (*.tif *.tiff *.dat)"

namespace geotools {

	namespace ui {

		namespace util {

			// Convenience: returns a QString from a std string.
			QString qstr(const std::string &str);

			// Convenience: returns a QString from an int.
			QString qstr(int val);

			// Hack to strip the boost::... part of an error message
			std::string stripBoost(const std::string &msg);

			void errorDialog(QWidget *parent, const std::string &title, const std::string &text, 
				const std::string &detail = "");

			// Convenience: open an input file dialog and return the string.
			std::string getInputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter);

			// Convenience: open an output file dialog and return the string.
			std::string getOutputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter);

			void getThresholds(QWidget *form, std::map<float, uint8_t> &thresholds);
		}

	}

}

#endif /* __UI_UTIL_HPP__ */
