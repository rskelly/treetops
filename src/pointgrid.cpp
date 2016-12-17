/*
 * This progam computes some basic statistics for differences between
 * pairs of rasters, aggregated by class.
 *
 * The classes come from a classification raster, which is the first argument.
 * Then, each of the given files is differenced with each other file
 * to produce pairwise mean, variance and standard deviation within
 * the class.
 * 
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <map>

#include "geotools.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "lasreader.hpp"
#include "laspoint.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::las;

void mean(const std::vector<std::string> &sourceFiles, const std::string &destFile, 
		const std::set<int> &classes, double resolution, double alignX, double alignY) {
	LASMultiReader lr(sourceFiles, resolution, -resolution);
	Bounds bounds = lr.bounds();
	bounds.align(alignX, alignY, resolution, -resolution);
	lr.setBounds(bounds);

	std::cerr << "bounds " << bounds.print() << "\n";
	Raster<float> out(destFile, 1, bounds, resolution, -resolution, 0);
	out.setNodata(-9999);
	MemRaster<float> sum(out.cols(), out.rows());
	sum.fill(0);
	MemRaster<uint32_t> count(out.cols(), out.rows());
	count.fill(0);

	LASPoint pt;
	bool file;
	int f = 1;
	while(lr.next(pt, nullptr, nullptr, &file)) {
		if(file)
			std::cerr << "file " << ++f << "\n";
		if(classes.find(pt.cls) == classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		sum.set(col, row, sum.get(col, row) + pt.z);
		count.set(col, row, count.get(col, row) + 1);
	}

	for(uint64_t i = 0; i < sum.size(); ++i) {
		int ct = count.get(i);
		sum.set(i, ct == 0 ? -9999 : sum.get(i) / ct);
	}

	out.writeBlock(sum);
}

void density(const std::vector<std::string> &sourceFiles, const std::string &destFile, 
		const std::set<int> &classes, double resolution, double alignX, double alignY) {
	LASMultiReader lr(sourceFiles, resolution, -resolution);
	Bounds bounds = lr.bounds();
	bounds.align(alignX, alignY, resolution, -resolution);
	lr.setBounds(bounds);

	Raster<float> out(destFile, 1, bounds, resolution, -resolution, 0);
	out.setNodata(-9999);
	MemRaster<float> sum(out.cols(), out.rows());
	sum.fill(0);
	MemRaster<uint32_t> count(out.cols(), out.rows());
	count.fill(0);

	LASPoint pt;
	bool file;
	int f = 1;
	while(lr.next(pt, nullptr, nullptr, &file)) {
		if(file)
			std::cerr << "file " << ++f << "\n";
		if(classes.find(pt.cls) == classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		count.set(col, row, count.get(col, row) + 1);
	}

	for(uint64_t i = 0; i < sum.size(); ++i) {
		int ct = count.get(i);
		sum.set(i, (float) ct / resolution * resolution);
	}

	out.writeBlock(sum);
}

void usage() {
	std::cerr << "usage\n";
}

int main(int argc, char ** argv) {

	if(argc == 1) {
		usage();
		return 1;
	}

	g_loglevel(0);

	using namespace geotools::raster;

	std::map<std::string, uint8_t> statTypes;
	statTypes["mean"] = 1;
	statTypes["median"] = 2;
	statTypes["stddev"] = 8;
	statTypes["min"] = 9;
	statTypes["max"] = 10;
	statTypes["variance"] = 11;
	statTypes["density"] = 12;

	std::set<int> classes;
	std::string destFile;
	std::vector<std::string> sourceFiles;
	int method;
	double resolution;
	double alignX = 0;
	double alignY = 0;

	try {

		for (int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-c") {
				Util::intSplit(classes, argv[++i]);
			} else if(arg == "-o") {
				destFile = argv[++i];
			} else if(arg == "-r") {
				resolution = atof(argv[++i]);
			} else if(arg == "-ax") {
				alignX = atof(argv[++i]);
			} else if(arg == "-ay") {
				alignY = atof(argv[++i]);
			} else if(arg == "-m") {
				std::string t = argv[++i];
				if(statTypes.find(t) == statTypes.end())
					g_argerr("Unknown method: " << t);
				method = statTypes[t];
			} else {
				sourceFiles.push_back(argv[i]);
			}
		}

		switch(method) {
		case 1:
			mean(sourceFiles, destFile, classes, resolution, alignX, alignY);
			break;
		case 12:
			density(sourceFiles, destFile, classes, resolution, alignX, alignY);
			break;
		default:
			g_argerr("Unknown method: " << method);
		}

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
