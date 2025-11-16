/*
 * ui_util.hpp
 *
 *  Created on: Dec 18, 2016
 *      Author: rob
 */

#ifndef __UI_UTIL_HPP__
#define __UI_UTIL_HPP__

#include <string>

#include <QtCore/QString>
#include <QtCore/QDir>
#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>

#include "treetops.hpp"

#define VECTOR_PATTERN "Vector Files (*.sqlite *.shp)"
#define RASTER_PATTERN "Raster Files (*.tif *.tiff *.dat)"
#define ALL_PATTERN "All Files (*)"

namespace geo {

	namespace ui {

		namespace util {

			// Convenience: returns a QString from a std string.
			QString G_DLL_EXPORT qstr(const std::string& str);

			// Convenience: returns a QString from an int.
			QString G_DLL_EXPORT qstr(int val);

			// Convenience: returns a std::string from a QString -- avoids ABI incompatibility issues.
			std::string G_DLL_EXPORT sstr(const QString& str);

			// Hack to strip the boost::... part of an error message
			std::string G_DLL_EXPORT stripBoost(const std::string& msg);

			/**
			 * Displays a modal dialog indicating an error.
			 *
			 * @param parent The parent widget.
			 * @param title The dialog title.
			 * @param text The dialog text.
			 * @param detail A detail string, optional.
			 */
			void G_DLL_EXPORT errorDialog(QWidget* parent, const std::string& title,
					const std::string& text, const std::string& detail = "");

			/**
			 * Displays an informative modal dialog.
			 *
			 * @param parent The parent widget.
			 * @param title The dialog title.
			 * @param text The dialog text.
			 * @param detail A detail string, optional.
			 */
			void G_DLL_EXPORT infoDialog(QWidget* parent, const std::string& title,
					const std::string& text, const std::string& detail = "");

			// Convenience: open an input file dialog and return the string.
			void G_DLL_EXPORT getInputFile(QWidget* form, const std::string& title, std::string& path,
					const std::string& filter, std::string& filename);

			// Convenience: open an output file dialog and return the string.
			void G_DLL_EXPORT getOutputFile(QWidget* form, const std::string& title, std::string& path,
					const std::string& filter, std::string& filename, bool confirmOverwrite = true);

			// Get the list of thresholds from the thresholds dialog.
			void G_DLL_EXPORT getTopsThresholds(QWidget* form, std::vector<geo::treetops::config::TopThreshold>& thresholds);

			void G_DLL_EXPORT getCrownsThresholds(QWidget* form, std::vector<geo::treetops::config::CrownThreshold>& thresholds);
		}

	}

}

#endif /* __UI_UTIL_HPP__ */
