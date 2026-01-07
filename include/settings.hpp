/*
 * settings.hpp
 *
 *  Created on: Jun 25, 2018
 *      Author: rob
 */

#ifndef _SETTINGS_HPP_
#define _SETTINGS_HPP_

#include <string>
#include <fstream>
#include <unordered_map>

#include <QObject>

namespace tt {
namespace config {

/**
 * Represents a configurable treetop threshold.
 */
class TopThreshold {
public:
	double threshold; 	///<! The minimum height for a pixel to be considered a top.
	int window;			///<! The size of the window used to search for maxima.

	/**
	 * Construct a TopThreshold instance with the minimum height threshold
	 * and kernel size.
	 *
	 * \param threshold The minimum top height threshold.
	 * \param window The size of the kernel used to search for maxima (an odd number from 1 to n).
	 */
	TopThreshold(double threshold, int window) :
		threshold(threshold), window(window) {}

	/**
	 * Construct a TopThreshold with defaults (0, 0).
	 */
	TopThreshold() : TopThreshold(0, 0) {}
};

/**
 * Represents a configurable tree crown threshold.
 */
class CrownThreshold {
public:
	double threshold;	///<! The minimum height for a pixel to be considered a member of a crown.
	double fraction;	///<! The minimum height as a proportion of the top height for a pixel to be considered a member of a crown.
	double radius;		///<! The maximum radius, in map units, of a crown.

	/**
	 * Construct a CrownThreshold with the given threshold, fraction and radius.
	 *
	 * \param threshold The minimum height for a pixel to be considered a member of a crown.
	 * \param fraction The minimum height as a proportion of the top height for a pixel to be considered a member of a crown.
	 * \param radius The maximum radius, in map units, of a crown.
	 */
	CrownThreshold(double threshold, double fraction, double radius) :
		threshold(threshold), fraction(fraction), radius(radius) {}

	/**
	 * Construct a default CrownThreshold with default values (0, 0, 0).
	 */
	CrownThreshold() : CrownThreshold(0, 0, 0) {}
};

typedef std::unordered_map<std::string, std::string> smap;

/**
 * A class for loading and saving settings.
 */
class Settings : public QObject {
	Q_OBJECT
private:
	smap m_settings;			///<! A map of local settings.
	std::string m_lastDir;		///<! The last-used directory.
	std::vector<CrownThreshold> m_crownThresholds;
	std::vector<TopThreshold> m_topThresholds;

signals:
	/**
	 * Emitted when the named setting has been updated.
	 */
	void settingsUpdate(const std::string&);

public:

	explicit Settings();

	/**
	 * Return true if the settings are such that the process can run.
	 */
	bool canRun();

	void parseTopThresholds(const std::string&);

	void parseCrownThresholds(const std::string&);

	void crownThresholds(const std::vector<CrownThreshold>&);

	const std::vector<CrownThreshold>& crownThresholds();

	void topThresholds(const std::vector<TopThreshold>&);

	const std::vector<TopThreshold>& topThresholds();

	const std::string& get(const std::string&, const char*);

	const std::string& get(const std::string&, const std::string&);

	bool get(const std::string&, bool);

	int get(const std::string&, int);

	float get(const std::string&, float);

	double get(const std::string&, double);

	void set(const std::string&, const char*);

	void set(const std::string&, const std::string&);

	void set(const std::string&, bool);

	void set(const std::string&, int);

	void set(const std::string&, float);

	void set(const std::string&, double);

	void settingsFile(const std::string&);

	std::string settingsFile();

	/**
	 * Load the settings from the file.
	 */
	void load();

	/**
	 * Save the settings to the file.
	 */
	void save();

	/**
	 * Return the last-used directory.
	 *
	 * \return The last-used directory.
	 */
	const std::string& lastDir();

	/**
	 * Set the last-used directory on the settings object. If the filename is a file,
	 * takes the parent directory. Otherwise uses the file itself. Returns the original
	 * filename as passed.
	 *
	 * \param path The directory or file path to use as the last-used directory.
	 * \return The directory path.
	 */
	void lastDir(const std::string&);

	~Settings();

};

} // config
} // tt

#endif /* _SETTINGS_HPP_ */
