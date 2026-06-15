#include "config.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include "util.hpp"

using namespace tt::config;
using namespace tt::util;

using json = nlohmann::json;

namespace nlohmann {

template <>
struct adl_serializer<TopThreshold> {
	static TopThreshold from_json(const json& j) {
		return {j.at("threshold").get<double>(), j.at("window").get<int>()};
	}
	static void to_json(json& j, TopThreshold t) {
		j = json{{"threshold", t.threshold}, {"window", t.window}};
	}
};

template <>
struct adl_serializer<CrownThreshold> {
	static CrownThreshold from_json(const json& j) {
		return {
			j.at("threshold").get<double>(),
			j.at("fraction").get<double>(),
			j.at("radius").get<double>()
		};
	}
	static void to_json(json& j, CrownThreshold t) {
		j = json{
			{"threshold", t.threshold},
			{"fraction", t.fraction},
			{"radius", t.radius}
		};
	}
};

} // namespace nlohmann

Config::Config() {
	setDefaults();
}

bool Config::canRun() const {
	return exists(get("originalCHM", ""));
}

void Config::setDefaults() {
	set("originalCHMBand", 1);
	set("doSmoothing", true);
	set("smoothWindowSize", 3);
	set("smoothSigma", 1.0f);
	set("doTops", true);
	set("doCrowns", true);
	set("topsMaxNulls", 0.5);
	set("crownsUpdateHeights", true);
	set("crownsDoDatabase", true);
	set("crownsRemoveHoles", true);
	set("crownsRemoveDangles", true);
	set("crownsKeepSmoothed", true);
	set("smoothedCHMDriver", "GTiff");
	set("crownsRasterDriver", "GTiff");
	set("treetopsDatabaseDriver", "Spatialite");
	set("crownsDatabaseDriver", "Spatialite");

	if(m_topThresholds.empty()) {
		m_topThresholds = {
			{3.0, 3},
			{5.0, 7}
		};
	}
	if(m_crownThresholds.empty()) {
		m_crownThresholds = {
			{5.0, 0.4, 5.0},
			{8.0, 0.4, 7.0}
		};
	}
}

void Config::deriveOutputPaths(const std::string& chmPath) {
	std::string dir = parent(chmPath);
	std::string stem = basename(chmPath);

	set("originalCHM", chmPath);
	set("smoothedCHM", join(dir, stem + "_smooth.tif"));
	set("treetopsDatabase", join(dir, "tops.sqlite"));
	set("crownsRaster", join(dir, stem + "_crowns.tif"));
	set("crownsDatabase", join(dir, stem + "_crowns.sqlite"));
	set("topsWindowsRaster", join(dir, "tops_windows.tif"));
	set("topsIdsRaster", join(dir, "tops_ids.tif"));
}

void Config::parseTopThresholds(const std::string& text) {
	m_topThresholds.clear();
	if(text.empty())
		return;

	std::stringstream ss(text);
	std::string item;
	while(std::getline(ss, item, ',')) {
		std::stringstream pair(item);
		std::string height;
		std::string window;
		if(!std::getline(pair, height, ':') || !std::getline(pair, window, ':'))
			continue;
		m_topThresholds.emplace_back(std::stod(height), std::stoi(window));
	}
}

void Config::parseCrownThresholds(const std::string& text) {
	m_crownThresholds.clear();
	if(text.empty())
		return;

	std::stringstream ss(text);
	std::string item;
	while(std::getline(ss, item, ',')) {
		std::stringstream triple(item);
		std::string height;
		std::string fraction;
		std::string radius;
		if(!std::getline(triple, height, ':')
				|| !std::getline(triple, fraction, ':')
				|| !std::getline(triple, radius, ':'))
			continue;
		m_crownThresholds.emplace_back(
			std::stod(height),
			std::stod(fraction),
			std::stod(radius));
	}
}

void Config::crownThresholds(const std::vector<CrownThreshold>& ct) {
	m_crownThresholds.assign(ct.begin(), ct.end());
}

const std::vector<CrownThreshold>& Config::crownThresholds() const {
	return m_crownThresholds;
}

void Config::topThresholds(const std::vector<TopThreshold>& tt) {
	m_topThresholds.assign(tt.begin(), tt.end());
}

const std::vector<TopThreshold>& Config::topThresholds() const {
	return m_topThresholds;
}

std::string Config::get(const std::string& k, const char* d) const {
	return get(k, std::string(d));
}

std::string Config::get(const std::string& k, const std::string& d) const {
	auto it = m_settings.find(k);
	if(it == m_settings.end())
		return d;
	return it->second;
}

bool Config::get(const std::string& k, bool d) const {
	auto it = m_settings.find(k);
	if(it == m_settings.end())
		return d;
	return it->second == "true";
}

int Config::get(const std::string& k, int d) const {
	auto it = m_settings.find(k);
	if(it == m_settings.end())
		return d;
	return std::stoi(it->second);
}

float Config::get(const std::string& k, float d) const {
	auto it = m_settings.find(k);
	if(it == m_settings.end())
		return d;
	return std::stof(it->second);
}

double Config::get(const std::string& k, double d) const {
	auto it = m_settings.find(k);
	if(it == m_settings.end())
		return d;
	return std::stod(it->second);
}

void Config::set(const std::string& k, const char* v) {
	m_settings[k] = v;
}

void Config::set(const std::string& k, const std::string& v) {
	m_settings[k] = v;
}

void Config::set(const std::string& k, bool v) {
	m_settings[k] = v ? "true" : "false";
}

void Config::set(const std::string& k, int v) {
	m_settings[k] = std::to_string(v);
}

void Config::set(const std::string& k, float v) {
	m_settings[k] = std::to_string(v);
}

void Config::set(const std::string& k, double v) {
	m_settings[k] = std::to_string(v);
}

void Config::load(const std::string& path) {
	std::ifstream str(path);
	if(!str)
		throw std::runtime_error("Failed to open settings file: " + path);

	json data = json::parse(str);

	for(json::iterator it = data.begin(); it != data.end(); ++it) {
		const std::string key = it.key();
		if(key == "topsThresholds") {
			m_topThresholds.clear();
			for(const json& entry : it.value())
				m_topThresholds.push_back(entry.get<TopThreshold>());
		} else if(key == "crownsThresholds") {
			m_crownThresholds.clear();
			for(const json& entry : it.value())
				m_crownThresholds.push_back(entry.get<CrownThreshold>());
		} else if(it.value().is_string()) {
			m_settings[key] = it.value().get<std::string>();
		}
	}
}

void Config::save(const std::string& path) const {
	std::ofstream str(path);
	if(!str)
		throw std::runtime_error("Failed to write settings file: " + path);

	json data;
	data["topsThresholds"] = json::array();
	for(const TopThreshold& t : m_topThresholds)
		data["topsThresholds"].push_back(json(t));

	data["crownsThresholds"] = json::array();
	for(const CrownThreshold& c : m_crownThresholds)
		data["crownsThresholds"].push_back(json(c));

	for(smap::const_iterator it = m_settings.begin(); it != m_settings.end(); ++it)
		data[it->first] = it->second;

	str << data.dump(2);
}

std::string Config::vectorDriver(const std::string& key) const {
	std::string driver = lowercase(get(key, "SQLite"));
	if(driver == "spatialite" || driver == "sqlite")
		return "SQLite";
	if(driver == "esri shapefile")
		return "ESRI Shapefile";
	return get(key, "SQLite");
}
