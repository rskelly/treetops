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
#include "pointgrid.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::las;
using namespace geotools::point;

bool __pg_cancel = false;

void min(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	cancel = cancel == nullptr ? &__pg_cancel : cancel;

	if(callbacks) {
		callbacks->stepCallback(0.01f);
		callbacks->statusCallback("Preparing...");
	}

	LASMultiReader lr(config.sourceFiles, config.resolutionX,
			config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY,
			config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds,
			config.resolutionX, config.resolutionY, config.srid);
	out.setNodata(-9999);

	MemRaster<float> min(out.cols(), out.rows());
	min.fill(-9999.0);

	std::vector<bool> visited((size_t) out.cols() * out.rows());
	std::fill(visited.begin(), visited.end(), false);

	if(callbacks) {
		callbacks->stepCallback(0.02f);
		callbacks->statusCallback("Processing...");
	}

	LASPoint pt;
	bool file;
	size_t f = 1;
	while(lr.next(pt, nullptr, nullptr, &file)) {
		if(file) {
			if(callbacks)
				callbacks->stepCallback(0.02f + ((float) f++ / config.sourceFiles.size()) * 0.96);
		}
		if(config.classes.find(pt.cls) == config.classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		size_t idx = (size_t) row * out.cols() + col;
		if(!visited[idx]) {
			min.set(idx, pt.z);
			visited[idx] = true;
		} else if(pt.z < min.get(idx)) {
			min.set(idx, pt.z);
		}
	}

	if(callbacks) {
		callbacks->stepCallback(0.99f);
		callbacks->statusCallback("Writing...");
	}

	out.writeBlock(min);

	if(callbacks) {
		callbacks->stepCallback(1.0f);
		callbacks->statusCallback("Done...");
	}

}

void max(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	cancel = cancel == nullptr ? &__pg_cancel : cancel;

	if(callbacks) {
		callbacks->stepCallback(0.01f);
		callbacks->statusCallback("Preparing...");
	}

	LASMultiReader lr(config.sourceFiles, config.resolutionX,
			config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY,
			config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds,
			config.resolutionX, config.resolutionY, config.srid);
	out.setNodata(-9999);

	MemRaster<float> max(out.cols(), out.rows());
	max.fill(-9999.0);

	std::vector<bool> visited((size_t) out.cols() * out.rows());
	std::fill(visited.begin(), visited.end(), false);

	if(callbacks) {
		callbacks->stepCallback(0.02f);
		callbacks->statusCallback("Processing...");
	}

	LASPoint pt;
	bool file;
	size_t f = 1;
	while(lr.next(pt, nullptr, nullptr, &file)) {
		if(file) {
			if(callbacks)
				callbacks->stepCallback(0.02f + ((float) f++ / config.sourceFiles.size()) * 0.96);
		}
		if(config.classes.find(pt.cls) == config.classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		size_t idx = (size_t) row * out.cols() + col;
		if(!visited[idx]) {
			max.set(idx, pt.z);
			visited[idx] = true;
		} else if(pt.z > max.get(idx)) {
			max.set(idx, pt.z);
		}
	}

	if(callbacks) {
		callbacks->stepCallback(0.99f);
		callbacks->statusCallback("Writing...");
	}

	out.writeBlock(max);

	if(callbacks) {
		callbacks->stepCallback(1.0f);
		callbacks->statusCallback("Done...");
	}

}

void mean(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	cancel = cancel == nullptr ? &__pg_cancel : cancel;

	if(callbacks) {
		callbacks->stepCallback(0.01f);
		callbacks->statusCallback("Preparing...");
	}

	LASMultiReader lr(config.sourceFiles, config.resolutionX,
			config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY,
			config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds,
			config.resolutionX, config.resolutionY, config.srid);
	out.setNodata(-9999);

	MemRaster<float> sum(out.cols(), out.rows());
	sum.fill(0);
	MemRaster<uint32_t> count(out.cols(), out.rows());
	count.fill(0);

	if(callbacks) {
		callbacks->stepCallback(0.02f);
		callbacks->statusCallback("Processing...");
	}

	LASPoint pt;
	bool file;
	size_t f = 1;
	while(lr.next(pt, nullptr, nullptr, &file)) {
		if(file) {
			if(callbacks)
				callbacks->stepCallback(0.02f + ((float) f++ / config.sourceFiles.size()) * 0.95);
		}
		if(config.classes.find(pt.cls) == config.classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		sum.set(col, row, sum.get(col, row) + pt.z);
		count.set(col, row, count.get(col, row) + 1);
	}

	if(callbacks) {
		callbacks->stepCallback(0.98f);
		callbacks->statusCallback("Computing...");
	}

	for(uint64_t i = 0; i < sum.size(); ++i) {
		int ct = count.get(i);
		sum.set(i, ct == 0 ? -9999 : sum.get(i) / ct);
	}

	if(callbacks) {
		callbacks->stepCallback(0.99f);
		callbacks->statusCallback("Writing...");
	}

	out.writeBlock(sum);

	if(callbacks) {
		callbacks->stepCallback(1.0f);
		callbacks->statusCallback("Done...");
	}

}

void variance(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	cancel = cancel == nullptr ? &__pg_cancel : cancel;

	if(callbacks) {
		callbacks->stepCallback(0.01f);
		callbacks->statusCallback("Preparing...");
	}

	LASMultiReader lr(config.sourceFiles, config.resolutionX,
			config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY,
			config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds,
			config.resolutionX, config.resolutionY, config.srid);
	out.setNodata(-9999);

	MemRaster<float> sum(out.cols(), out.rows());
	sum.fill(0);

	if(callbacks) {
		callbacks->stepCallback(0.02f);
		callbacks->statusCallback("Processing...");
	}

	LASPoint pt;
	bool file;
	size_t f = 1;
	{

		MemRaster<uint32_t> count(out.cols(), out.rows());
		count.fill(0);

		while(lr.next(pt, nullptr, nullptr, &file)) {
			if(file) {
				if(callbacks)
					callbacks->stepCallback(0.02f + ((float) f++ / config.sourceFiles.size()) * 0.45);
			}
			if(config.classes.find(pt.cls) == config.classes.end())
				continue;
			int col = out.toCol(pt.x);
			int row = out.toRow(pt.y);
			sum.set(col, row, sum.get(col, row) + pt.z);
			count.set(col, row, count.get(col, row) + 1);
		}

		if(callbacks) {
			callbacks->stepCallback(0.48f);
			callbacks->statusCallback("Computing...");
		}

		for(uint64_t i = 0; i < sum.size(); ++i) {
			int ct = count.get(i);
			sum.set(i, ct == 0 ? -9999.0f : sum.get(i) / ct);
		}

	}

	{
		MemRaster<float> var(out.cols(), out.rows());
		var.fill(-9999.0);

		f = 1;
		lr.reset();
		while(lr.next(pt, nullptr, nullptr, &file)) {
			if(file) {
				if(callbacks)
					callbacks->stepCallback(0.48f + ((float) f++ / config.sourceFiles.size()) * 0.50);
			}
			if(config.classes.find(pt.cls) == config.classes.end())
				continue;
			int col = out.toCol(pt.x);
			int row = out.toRow(pt.y);
			double v;
			if((v = sum.get(col, row)) != -9999.0)
				var.set(col, row, g_sq(pt.z - v));
		}

		if(callbacks) {
			callbacks->stepCallback(0.99f);
			callbacks->statusCallback("Writing...");
		}

		out.writeBlock(var);

	}


	if(callbacks) {
		callbacks->stepCallback(1.0f);
		callbacks->statusCallback("Done...");
	}

}

void stddev(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	cancel = cancel == nullptr ? &__pg_cancel : cancel;

	if(callbacks) {
		callbacks->stepCallback(0.01f);
		callbacks->statusCallback("Preparing...");
	}

	LASMultiReader lr(config.sourceFiles, config.resolutionX,
			config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY,
			config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds,
			config.resolutionX, config.resolutionY, config.srid);
	out.setNodata(-9999);

	MemRaster<float> sum(out.cols(), out.rows());
	sum.fill(0);

	if(callbacks) {
		callbacks->stepCallback(0.02f);
		callbacks->statusCallback("Processing...");
	}

	LASPoint pt;
	bool file;
	size_t f = 1;

	{

		MemRaster<uint32_t> count(out.cols(), out.rows());
		count.fill(0);

		while(lr.next(pt, nullptr, nullptr, &file)) {
			if(file) {
				if(callbacks)
					callbacks->stepCallback(0.02f + ((float) f++ / config.sourceFiles.size()) * 0.45);
			}
			if(config.classes.find(pt.cls) == config.classes.end())
				continue;
			int col = out.toCol(pt.x);
			int row = out.toRow(pt.y);
			sum.set(col, row, sum.get(col, row) + pt.z);
			count.set(col, row, count.get(col, row) + 1);
		}

		if(callbacks) {
			callbacks->stepCallback(0.48f);
			callbacks->statusCallback("Computing...");
		}

		for(uint64_t i = 0; i < sum.size(); ++i) {
			int ct = count.get(i);
			sum.set(i, ct == 0 ? -9999 : sum.get(i) / ct);
		}

	}

	{
		MemRaster<float> var(out.cols(), out.rows());
		var.fill(-9999);

		f = 1;

		lr.reset();
		while(lr.next(pt, nullptr, nullptr, &file)) {
			if(file) {
				if(callbacks)
					callbacks->stepCallback(0.48f + ((float) f++ / config.sourceFiles.size()) * 0.50);
			}
			if(config.classes.find(pt.cls) == config.classes.end())
				continue;
			int col = out.toCol(pt.x);
			int row = out.toRow(pt.y);
			double v;
			if((v = sum.get(col, row)) != -9999.0)
				var.set(col, row, std::sqrt(g_sq(pt.z - v)));
		}

		if(callbacks) {
			callbacks->stepCallback(0.99f);
			callbacks->statusCallback("Writing...");
		}

		out.writeBlock(var);

	}


	if(callbacks) {
		callbacks->stepCallback(1.0f);
		callbacks->statusCallback("Done...");
	}

}

void density(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	LASMultiReader lr(config.sourceFiles, config.resolutionX, config.resolutionY);
	Bounds bounds = lr.bounds();
	bounds.align(config.alignX, config.alignY, config.resolutionX, config.resolutionY);
	lr.setBounds(bounds);

	Raster<float> out(config.destFile, 1, bounds, config.resolutionX, config.resolutionY, config.srid);
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
		if(config.classes.find(pt.cls) == config.classes.end())
			continue;
		int col = out.toCol(pt.x);
		int row = out.toRow(pt.y);
		count.set(col, row, count.get(col, row) + 1);
	}

	float res = g_sq(config.resolutionX);
	for(uint64_t i = 0; i < sum.size(); ++i) {
		int ct = count.get(i);
		sum.set(i, (float) ct / res);
	}

	out.writeBlock(sum);
}

std::map<std::string, uint8_t> statTypes;

void PointGrid::process(const PointgridConfig &config, Callbacks *callbacks, bool *cancel) {

	switch(config.method) {
	case 1:
		mean(config, callbacks, cancel);
		break;
	case 2:
		stddev(config, callbacks, cancel);
		break;
	case 3:
		variance(config, callbacks, cancel);
		break;
	case 4:
		min(config, callbacks, cancel);
		break;
	case 5:
		max(config, callbacks, cancel);
		break;
	case 6:
		density(config, callbacks, cancel);
		break;
	default:
		g_argerr("Unknown method: " << config.method);
	}
}

void usage() {
	std::cerr << "This program computes simple cell stats for point clouds.\n"
			<< "Usage: pointgrid <options> <source file [source file [...]]>\n"
			<< " -m       The method. One of:\n";
			for(const auto &it : statTypes)
				std::cerr << "          " << it.first << "\n";
	std::cerr << " -o       The output file.\n"
			<< " -ax, -ay The vertical and horizontal alignment.\n"
			<< " -rx, -ry The vertical and horizontal resolution.\n"
			<< " -c       A comma-separated list of classes.\n"
			<< " -s       The SRID\n";
}

int main(int argc, char ** argv) {

	statTypes["mean"] = 1;
	statTypes["stddev"] = 2;
	statTypes["variance"] = 3;
	statTypes["min"] = 4;
	statTypes["max"] = 5;
	statTypes["density"] = 6;

	if(argc == 1) {
		usage();
		return 1;
	}

	g_loglevel(0);

	using namespace geotools::raster;

	PointgridConfig config;

	try {

		for (int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-c") {
				Util::intSplit(config.classes, argv[++i]);
			} else if(arg == "-o") {
				config.destFile = argv[++i];
			} else if(arg == "-rx") {
				config.resolutionX = atof(argv[++i]);
			} else if(arg == "-ry") {
				config.resolutionY = atof(argv[++i]);
			} else if(arg == "-ax") {
				config.alignX = atof(argv[++i]);
			} else if(arg == "-ay") {
				config.alignY = atof(argv[++i]);
			} else if(arg == "-m") {
				std::string t = argv[++i];
				if(statTypes.find(t) == statTypes.end())
					g_argerr("Unknown method: " << t);
				config.method = statTypes[t];
			} else if(arg == "-s") {
				config.srid = atoi(argv[++i]);
			} else {
				config.sourceFiles.push_back(argv[i]);
			}
		}

		PointGrid pg;
		pg.process(config);

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
