/*
 * settings.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: rob
 */

#include <fstream>
#include <unordered_map>

#include <QtCore/QSettings>

#include "settings.hpp"
#include "ui_util.hpp"

using namespace geo::treetops::config;

namespace {

	/**
	 * Load the key-value file contents into a map.
	 *
	 * \param filename The target file. Will be overwritten.
	 * \param map The map.
	 */
	bool loadMap(const std::string& filename, smap& map) {
		if (filename.empty())
			return false;
		std::ifstream ins(filename, std::ios::in);
		if (ins.bad())
			return false;
		constexpr int len = 1024;
		char key[len];
		char val[len];
		while (ins.good()) {
			ins.getline(key, len, ':');
			ins.getline(val, len, '\n');
			map[key] = val;
		}
		return true;
	}

	/**
	 * Save the map as a key-value file.
	 *
	 * \param filename The target file.
	 * \param map The map.
	 */
	void saveMap(const std::string& filename, smap& map) {
		if (!filename.empty()) {
			std::ofstream ofs(filename, std::ios::out);
			if (ofs.good()) {
				for (auto& it : map)
					ofs << it.first << ":" << it.second << "\n";
			}
		}
	}

	/**
	 * Return the map value corresponding to the given key as an integer.
	 * Return the alternate if the key doesn't exist.
	 *
	 * \param map The map.
	 * \param key The key.
	 * \param alt The alternate value.
	 */
	int geti(const smap& map, const std::string& key, int alt) {
		if (map.find(key) == map.end()) {
			return alt;
		}
		else {
			return atoi(map.at(key).c_str());
		}
	}

	/**
	 * Return the map value corresponding to the given key as a double.
	 * Return the alternate if the key doesn't exist.
	 *
	 * \param map The map.
	 * \param key The key.
	 * \param alt The alternate value.
	 */
	double getf(const smap& map, const std::string& key, double alt) {
		if (map.find(key) == map.end()) {
			return alt;
		}
		else {
			return atof(map.at(key).c_str());
		}
	}

	/**
	 * Return the map value corresponding to the given key as a boolean.
	 * Return the alternate if the key doesn't exist.
	 * Allowed values are anything that starts with 't' or 'T' or '1'.
	 *
	 * \param map The map.
	 * \param key The key.
	 * \param alt The alternate value.
	 */
	bool getb(const smap& map, const std::string& key, bool alt) {
		if (map.find(key) == map.end()) {
			return alt;
		}
		else {
			const std::string& val = map.at(key);
			return val[0] == 't' || val[0] == 'T' || val[0] == '1';
		}
	}

	/**
	 * Return the map value corresponding to the given key as a string.
	 * Return the alternate if the key doesn't exist.
	 *
	 * \param map The map.
	 * \param key The key.
	 * \param alt The alternate value.
	 */
	std::string gets(const smap& map, const std::string& key, const std::string& alt) {
		if (map.find(key) == map.end()) {
			return alt;
		}
		else {
			return map.at(key);
		}
	}

} // anon

Settings::Settings() :
	m_settings(new QSettings("ui_settings.txt", QSettings::Format::IniFormat)),
	m_lastDir("") {
	m_lastDir = geo::ui::util::sstr(m_settings->value("local/lastDir", "").toString());
}

std::string& Settings::lastDir() {
	return m_lastDir;
}

bool Settings::load(TreetopsConfig& config, const std::string& filename) {
	m_settings->setValue("local/settings", QString(filename.c_str()));

	smap map;
	if(!loadMap(filename, map))
		return false;

	config.setBuildIndex(getb(map, "buildIndex", config.buildIndex()));
	config.setTableCacheSize(geti(map, "tableCacheSize", config.tableCacheSize()));
	config.setRowCacheSize(geti(map, "rowCacheSize", config.rowCacheSize()));

	config.setOriginalCHM(gets(map, "originalCHM", config.originalCHM()));
	config.setOriginalCHMBand(geti(map, "originalCHMBand", config.originalCHMBand()));
	config.setSmoothedCHM(gets(map, "smoothedCHM", config.smoothedCHM()));
	config.setSmoothedCHMDriver(gets(map, "smoothedCHMDriver", config.smoothedCHMDriver()));
	config.setTreetopsDatabase(gets(map, "treetopsDatabase", config.treetopsDatabase()));
	config.setTreetopsDatabaseDriver(gets(map, "treetopsDatabaseDriver", config.treetopsDatabaseDriver()));
	config.setCrownsRaster(gets(map, "crownsRaster", config.crownsRaster()));
	config.setCrownsRasterDriver(gets(map, "crownsRasterDriver", config.crownsRasterDriver()));
	config.setCrownsDatabase(gets(map, "crownsDatabase", config.crownsDatabase()));
	config.setCrownsDatabaseDriver(gets(map, "crownsDatabaseDriver", config.crownsDatabaseDriver()));

	config.setDoSmoothing(getb(map, "doSmoothing", config.doSmoothing()));
	config.setSmoothWindowSize(geti(map, "smoothWindowSize", config.smoothWindowSize()));
	config.setSmoothSigma(getf(map, "smoothSigma", config.smoothSigma()));

	config.setDoTops(getb(map, "doTops", config.doTops()));
	config.parseTopsThresholds(gets(map, "topsThresholds", ""));
	config.setTopsMaxNulls(getf(map, "topsMaxNulls", config.topsMaxNulls()));

	config.setDoCrowns(getb(map, "doCrowns", config.doCrowns()));
	config.parseCrownsThresholds(gets(map, "crownsThresholds", ""));
	config.setCrownsUpdateHeights(getb(map, "crownsUpdateHeights", config.crownsUpdateHeights()));
	config.setCrownsDoDatabase(getb(map, "crownsDoDatabase", config.crownsDoDatabase()));
	config.setCrownsRemoveHoles(getb(map, "crownsRemoveHoles", config.crownsRemoveHoles()));
	config.setCrownsRemoveDangles(getb(map, "crownsRemoveDangles", config.crownsRemoveDangles()));
	config.setCrownsKeepSmoothed(getb(map, "crownsKeepSmoothed", config.crownsKeepSmoothed()));

	return true;
}

void Settings::save(TreetopsConfig& config) {
	smap map;
	map["settingsFile"] = config.settings();
	map["buildIndex"] = std::to_string(config.buildIndex());
	map["tableCacheSize"] = std::to_string(config.tableCacheSize());
	map["rowCacheSize"] = std::to_string(config.rowCacheSize()); //(24 * 1024 * 1024),

	map["originalCHM"] = config.originalCHM();
	map["originalCHMBand"] = std::to_string(config.originalCHMBand());
	map["smoothedCHM"] = config.smoothedCHM();
	map["smoothedCHMDriver"] = config.smoothedCHMDriver();
	map["treetopsDatabase"] = config.treetopsDatabase();
	map["treetopsDatabaseDriver"] = config.treetopsDatabaseDriver();
	map["crownsRaster"] = config.crownsRaster();
	map["crownsRasterDriver"] = config.crownsRasterDriver();
	map["crownsDatabase"] = config.crownsDatabase();
	map["crownsDatabaseDriver"] = config.crownsDatabaseDriver();

	map["doSmoothing"] = std::to_string(config.doSmoothing());
	map["smoothWindowSize"] = std::to_string(config.smoothWindowSize());
	map["smoothSigma"] = std::to_string(config.smoothSigma());

	map["doTops"] = std::to_string(config.doTops());
	map["topsThresholds"] = config.topsThresholdsList();
	map["topsMaxNulls"] = std::to_string(config.topsMaxNulls());

	map["doCrowns"] = std::to_string(config.doCrowns());
	map["crownsThresholds"] = config.crownsThresholdsList();
	map["crownsUpdateHeights"] = std::to_string(config.crownsUpdateHeights());
	map["crownsDoDatabase"] = std::to_string(config.crownsDoDatabase());
	map["crownsRemoveHoles"] = std::to_string(config.crownsRemoveHoles());
	map["crownsRemoveDangles"] = std::to_string(config.crownsRemoveDangles());
	map["crownsKeepSmoothed"] = std::to_string(config.crownsKeepSmoothed());

	saveMap(config.settings(), map);
}

Settings::~Settings() {
	m_settings->setValue("local/lastDir", QString(lastDir().c_str()));
	delete m_settings;
}

