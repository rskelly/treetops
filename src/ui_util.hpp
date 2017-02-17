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
#define ALL_PATTERN "All Files (*)"

namespace geo {

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
			void getInputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter, std::string &filename);

			// Convenience: open an output file dialog and return the string.
			void getOutputFile(QWidget *form, const std::string &title, QDir &path,
					const std::string &filter, std::string &filename);

			// Get the list of thresholds from the thresholds dialog.
			void getTopsThresholds(QWidget *form, std::vector<std::tuple<double, uint8_t> > &thresholds);

			void getCrownsThresholds(QWidget *form, std::vector<std::tuple<double, double, double> > &thresholds);
		}

	}

}

#endif /* __UI_UTIL_HPP__ */
