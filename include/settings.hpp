/*
 * settings.hpp
 *
 *  Created on: Jun 25, 2018
 *      Author: rob
 */

#ifndef LIBTREETOPS_INCLUDE_SETTINGS_HPP_
#define LIBTREETOPS_INCLUDE_SETTINGS_HPP_

#include <string>

#include <QtCore/QSettings>

#include "treetops.hpp"

namespace geo {
namespace treetops {
namespace config {

typedef std::unordered_map<std::string, std::string> smap;

/**
 * A class for loading and saving settings.
 */
class Settings {
private:
	QSettings* m_settings;	///<! The settings storage object.
	smap m_localSettings;	///<! A map of local settings.
	smap m_userSettings;	///<! A map of transient user settings.
	std::string m_lastDir;	///<! The last-used directory.

public:

	Settings();

	/**
	 * Return the last-used directory.
	 *
	 * \return The last-used directory.
	 */
	std::string& lastDir();

	/**
	 * Load the settings contained in filename into the TreetopsConfig object.
	 * If false is returned, there is no settings file available
	 *
	 * \param config A TreetopsConfig instance.
	 * \param filename A path to a settings file.
	 * \return False if no file is available or there's a failure. True otherwise.
	 */
	bool load(geo::treetops::config::TreetopsConfig& config, const std::string& filename);

	/**
	 * Save the settings contained in the TreetopsConfig object to the
	 * settings file contained in the config object.
	 *
	 * \param config A TreetopsConfig instance.
	 */
	void save(geo::treetops::config::TreetopsConfig& config);

	~Settings();
};

} // config
} // treetops
} // geo

#endif /* LIBTREETOPS_INCLUDE_SETTINGS_HPP_ */
