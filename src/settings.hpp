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
	QSettings m_settings;
	std::string m_lastFile;
	smap m_localSettings;
	smap m_userSettings;

public:
	std::string topsDatabaseLastDir;
	std::string originalCHMLastDir;
	std::string smoothedCHMLastDir;
	std::string crownsRasterLastDir;
	std::string crownsDatabaseLastDir;
	std::string outputLastDir;

	Settings();

	/**
	 * Return the most recently used settings file.
	 *
	 * @return The most recently used settings file.
	 */
	const std::string& lastFile() const;

	/**
	 * Load the settings from the previously-accessed settings file.
	 * If there isn't one, do nothing and return false.
	 */
	bool load(geo::treetops::config::TreetopsConfig& config);

	/**
	 * Load the settings contained in filename into the TreetopsConfig object.
	 * If false is returned, there is no settings file available
	 *
	 * @param config A TreetopsConfig instance.
	 * @param filename A path to a settings file.
	 * @return False if no file is available or there's a failure. True otherwise.
	 */
	bool load(geo::treetops::config::TreetopsConfig& config, const std::string& filename);

	/**
	 * Save the settings contained in the TreetopsConfig object to the filename.
	 *
	 * @param config A TreetopsConfig instance.
	 */
	void save(geo::treetops::config::TreetopsConfig& config);

	~Settings();
};

} // config
} // treetops
} // geos

#endif /* LIBTREETOPS_INCLUDE_SETTINGS_HPP_ */
