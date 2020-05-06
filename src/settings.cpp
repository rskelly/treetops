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

using namespace geo::treetops::config;

/**
 * Load the key-value file contents into a map.
 *
 * \param filename The target file. Will be overwritten.
 * \param map The map.
 */
bool _load(const std::string& filename, smap& map) {
	if(filename.empty())
		return false;
	std::ifstream ins(filename, std::ios::in);
	if(ins.bad())
		return false;
	constexpr int len = 1024;
	char key[len];
	char val[len];
	while(ins.good()) {
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
void _save(const std::string& filename, smap& map) {
	if(!filename.empty()) {
		std::ofstream ofs(filename, std::ios::out);
		if(ofs.good()) {
			for(auto& it : map)
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
int _geti(const smap& map, const std::string& key, int alt) {
	if(map.find(key) == map.end()) {
		return alt;
	} else {
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
double _getf(const smap& map, const std::string& key, double alt) {
	if(map.find(key) == map.end()) {
		return alt;
	} else {
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
bool _getb(const smap& map, const std::string& key, bool alt) {
	if(map.find(key) == map.end()) {
		return alt;
	} else {
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
std::string _gets(const smap& map, const std::string& key, const std::string& alt) {
	if(map.find(key) == map.end()) {
		return alt;
	} else {
		return map.at(key);
	}
}

Settings::Settings() {
	m_lastDir = m_settings.value("local/lastDir", "").toString().toStdString();
}

std::string& Settings::lastDir() {
	return m_lastDir;
}

bool Settings::load(TreetopsConfig& config, const std::string& filename) {
	m_settings.setValue("local/settings", QString(filename.c_str()));

	smap map;
	if(!_load(filename, map))
		return false;

	config.setBuildIndex(_getb(map, "buildIndex", config.buildIndex()));
	config.setTableCacheSize(_geti(map, "tableCacheSize", config.tableCacheSize()));
	config.setRowCacheSize(_geti(map, "rowCacheSize", config.rowCacheSize()));

	config.setOriginalCHM(_gets(map, "originalCHM", config.originalCHM()));
	config.setOriginalCHMBand(_geti(map, "originalCHMBand", config.originalCHMBand()));
	config.setSmoothedCHM(_gets(map, "smoothedCHM", config.smoothedCHM()));
	config.setSmoothedCHMDriver(_gets(map, "smoothedCHMDriver", config.smoothedCHMDriver()));
	config.setTreetopsDatabase(_gets(map, "treetopsDatabase", config.treetopsDatabase()));
	config.setTreetopsDatabaseDriver(_gets(map, "treetopsDatabaseDriver", config.treetopsDatabaseDriver()));
	config.setCrownsRaster(_gets(map, "crownsRaster", config.crownsRaster()));
	config.setCrownsRasterDriver(_gets(map, "crownsRasterDriver", config.crownsRasterDriver()));
	config.setCrownsDatabase(_gets(map, "crownsDatabase", config.crownsDatabase()));
	config.setCrownsDatabaseDriver(_gets(map, "crownsDatabaseDriver", config.crownsDatabaseDriver()));

	config.setDoSmoothing(_getb(map, "doSmoothing", config.doSmoothing()));
	config.setSmoothWindowSize(_geti(map, "smoothWindowSize", config.smoothWindowSize()));
	config.setSmoothSigma(_getf(map, "smoothSigma", config.smoothSigma()));

	config.setDoTops(_getb(map, "doTops", config.doTops()));
	config.parseTopsThresholds(_gets(map, "topsThresholds", ""));
	config.setTopsMaxNulls(_getf(map, "topsMaxNulls", config.topsMaxNulls()));

	config.setDoCrowns(_getb(map, "doCrowns", config.doCrowns()));
	config.parseCrownsThresholds(_gets(map, "crownsThresholds", ""));
	config.setCrownsUpdateHeights(_getb(map, "crownsUpdateHeights", config.crownsUpdateHeights()));
	config.setCrownsDoDatabase(_getb(map, "crownsDoDatabase", config.crownsDoDatabase()));
	config.setCrownsRemoveHoles(_getb(map, "crownsRemoveHoles", config.crownsRemoveHoles()));
	config.setCrownsRemoveDangles(_getb(map, "crownsRemoveDangles", config.crownsRemoveDangles()));

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

	_save(config.settings(), map);
}

Settings::~Settings() {
	m_settings.setValue("local/lastDir", QString(lastDir().c_str()));
}

