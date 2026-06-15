#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace tt {
namespace config {

class TopThreshold {
public:
	double threshold;
	int window;

	TopThreshold(double threshold, int window) :
		threshold(threshold), window(window) {}

	TopThreshold() : TopThreshold(0, 0) {}
};

class CrownThreshold {
public:
	double threshold;
	double fraction;
	double radius;

	CrownThreshold(double threshold, double fraction, double radius) :
		threshold(threshold), fraction(fraction), radius(radius) {}

	CrownThreshold() : CrownThreshold(0, 0, 0) {}
};

typedef std::unordered_map<std::string, std::string> smap;

/**
 * Runtime configuration for the treetops processing pipeline.
 * Loads and saves JSON settings without any UI dependencies.
 */
class Config {
private:
	smap m_settings;
	std::vector<CrownThreshold> m_crownThresholds;
	std::vector<TopThreshold> m_topThresholds;

public:
	Config();

	bool canRun() const;

	void parseTopThresholds(const std::string&);
	void parseCrownThresholds(const std::string&);

	void crownThresholds(const std::vector<CrownThreshold>&);
	const std::vector<CrownThreshold>& crownThresholds() const;

	void topThresholds(const std::vector<TopThreshold>&);
	const std::vector<TopThreshold>& topThresholds() const;

	std::string get(const std::string&, const char*) const;
	std::string get(const std::string&, const std::string&) const;
	bool get(const std::string&, bool) const;
	int get(const std::string&, int) const;
	float get(const std::string&, float) const;
	double get(const std::string&, double) const;

	void set(const std::string&, const char*);
	void set(const std::string&, const std::string&);
	void set(const std::string&, bool);
	void set(const std::string&, int);
	void set(const std::string&, float);
	void set(const std::string&, double);

	void setDefaults();
	void deriveOutputPaths(const std::string& chmPath);

	void load(const std::string& path);
	void save(const std::string& path) const;

	std::string vectorDriver(const std::string& key) const;
};

} // config
} // tt

#endif
