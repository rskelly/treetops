/*
 * settings.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: rob
 */

#include <fstream>
#include <unordered_map>

#include <QSettings>
#include <QDir>

#include <nlohmann/json.hpp>

#include "settings.hpp"
#include "util.hpp"

using namespace tt::config;
using namespace tt::util;

using json = nlohmann::json;

namespace nlohmann {

template <>
struct adl_serializer<TopThreshold> {

    static TopThreshold from_json(const json& j) {
        return {j.at("threshold"), j.at("window")};
    }

    static void to_json(json& j, TopThreshold t) {
        j["threshold"] = t.threshold;
        j["window"] = t.window;
    }
};

template <>
struct adl_serializer<CrownThreshold> {

    static CrownThreshold from_json(const json& j) {
        return {j.at("threshold"), j.at("fraction"), j.at("radius")};
    }

    static void to_json(json& j, CrownThreshold t) {
        j["threshold"] = t.threshold;
        j["fraction"] = t.fraction;
		j["radius"] = t.radius;
    }

};

} // namespace nlohmann

Settings::Settings() :
	QObject(),
	m_lastDir("") {
	load();
}

void Settings::parseTopThresholds(const std::string& t) {

}

void Settings::parseCrownThresholds(const std::string& t) {

}

void Settings::crownThresholds(const std::vector<CrownThreshold>& ct) {
	m_crownThresholds.assign(ct.begin(), ct.end());
	emit settingsUpdate("crownThresholds");
}

const std::vector<CrownThreshold>& Settings::crownThresholds() {
	return m_crownThresholds;
}

void Settings::topThresholds(const std::vector<TopThreshold>& tt) {
	m_topThresholds.assign(tt.begin(), tt.end());
	emit settingsUpdate("topThresholds");
}

const std::vector<TopThreshold>& Settings::topThresholds() {
	return m_topThresholds;
}

const std::string& Settings::get(const std::string& k, const char* d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = d;
		emit settingsUpdate(k.c_str());
	}
	return m_settings[k];
}

const std::string& Settings::get(const std::string& k, const std::string& d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = d;
		emit settingsUpdate(k.c_str());
	}
	return m_settings[k];
}

bool Settings::get(const std::string& k, bool d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = std::to_string(d);
		emit settingsUpdate(k.c_str());
	}
	return m_settings[k] == "true";
}

int Settings::get(const std::string& k, int d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = std::to_string(d);
		emit settingsUpdate(k.c_str());
	}
	return std::stoi(m_settings[k]);
}

float Settings::get(const std::string& k, float d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = std::to_string(d);
		emit settingsUpdate(k.c_str());
	}
	return std::stof(m_settings[k]);
}

double Settings::get(const std::string& k, double d) {
	if(m_settings.count(k) == 0) {
		m_settings[k] = std::to_string(d);
		emit settingsUpdate(k.c_str());
	}
	return std::stod(m_settings[k]);
}

void Settings::set(const std::string& k, const std::string& v) {
	m_settings[k] = v;
	emit settingsUpdate(k.c_str());
}

void Settings::set(const std::string& k, const char* v) {
	m_settings[k] = v;
	emit settingsUpdate(k.c_str());
}

void Settings::set(const std::string& k, bool v) {
	m_settings[k] = v == true ? "true" : "false";
	emit settingsUpdate(k.c_str());
}

void Settings::set(const std::string& k, int v) {
	m_settings[k] = std::to_string(v);
	emit settingsUpdate(k.c_str());
}

void Settings::set(const std::string& k, float v) {
	m_settings[k] = std::to_string(v);
	emit settingsUpdate(k.c_str());
}

void Settings::set(const std::string& k, double v) {
	m_settings[k] = std::to_string(v);
	emit settingsUpdate(k.c_str());
}

void Settings::settingsFile(const std::string& file) {
	// Save to the existing configuration, set the file, 
	// then save again to create the new file.
	save();
	QSettings settings("treetops.ini", QSettings::NativeFormat);
	settings.setValue("settingsFile", QString(file.c_str()));
	load();
	emit settingsUpdate("settingsFile");
}

std::string Settings::settingsFile() {
	QSettings settings("treetops.ini", QSettings::NativeFormat);
	return settings.value("settingsFile", QDir::home().filePath("tt_settings.json")).toString().toStdString();
}

void Settings::load() {
	std::string path = settingsFile();
	try {
		std::ifstream str(path);
		json data = json::parse(str);
		for (json::iterator it = data.begin(); it != data.end(); ++it) {
			std::string key = it.key();
			if(key == "topsThresholds") {

				m_topThresholds.clear();
				json thresh = it.value();
				for(json::iterator it_ = thresh.begin(); it_ != thresh.end(); ++it_)
					m_topThresholds.push_back(it_.value().get<TopThreshold>());

			} else if(it.key() == "crownsThresholds") {

				m_crownThresholds.clear();
				json thresh = it.value();
				for(json::iterator it_ = thresh.begin(); it_ != thresh.end(); ++it_)
					m_crownThresholds.push_back(it_.value().get<CrownThreshold>());

			} else {

				m_settings[it.key()] = it.value();

			}

			emit settingsUpdate(it.key());

		}
	} catch(const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
}

void Settings::save() {
	std::ofstream str(settingsFile());
	json data;

	data["topsThresholds"] = json::array();
	for(std::vector<TopThreshold>::iterator it = m_topThresholds.begin(); it < m_topThresholds.end(); ++it)
		data["topsThresholds"].push_back(json(*it));
	
		data["crownsThresholds"] = json::array();
	for(std::vector<CrownThreshold>::iterator it = m_crownThresholds.begin(); it < m_crownThresholds.end(); ++it)
		data["crownsThresholds"].push_back(json(*it));
	
	for(smap::iterator it = m_settings.begin(); it != m_settings.end(); ++it)
		data[it->first] = it->second;
	
	str << data.dump();
}


void Settings::lastDir(const std::string& path) {
	std::string dir = path;
	if (isfile(dir))
		dir = tt::util::parent(dir);
	m_settings["lastDir"] = dir;
	emit settingsUpdate("lastDir");
}


const std::string& Settings::lastDir() {
	return m_settings["lastDir"];
}

/*
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
*/

Settings::~Settings() {
	save();
}

