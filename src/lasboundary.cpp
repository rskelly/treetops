/*
 * lasboundary.cpp
 *
 * This tool creates a shapefile representing the boundary of a LiDAR
 * point cloud as represented by an alpha shape. See CGAL
 * alpha shape docs for more info.
 *
 *  Created on: Mar 13, 2015
 *      Author: rob
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <math.h>

#include <omp.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "geotools.hpp"
#include "util.hpp"
#include "lasutil.hpp"
#include "raster.hpp"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

#include "raster.hpp"
#include "lasreader.hpp"

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

	namespace las {

		void lasboundary(const std::vector<std::string> &srcFiles, const std::string &dstFile,
				const std::string &polyFile, int srid, double res, std::set<int> &classes) {

			if (srcFiles.empty())
				g_argerr("No source dir given.");
			if (dstFile.empty())
				g_argerr("No dest file given.");

			g_loglevel(G_LOG_DEBUG);
			g_debug("reader");
			LASMultiReader lr(srcFiles, res, res);
			Bounds bounds = lr.bounds();

			g_debug("raster");
			MemRaster<uint16_t> rast(bounds.maxCol(res) + 1, bounds.maxRow(res) + 1);
			rast.fill(0);
			rast.setNodata(0);

			g_debug("reading");
			LASPoint pt;
			while(lr.next(pt)) {
				if(classes.find(pt.cls) == classes.end()) continue;
				int col = bounds.toCol(pt.x, res);
				int row = bounds.toRow(pt.y, res);
				rast.set(col, row, rast.get(col, row) + 1);
			}

			for(uint64_t i = 0; i < rast.size(); ++i)
				rast.set(i, rast.get(i) <= 2 ? 0 : 1);

			g_debug("filling");
			for(uint16_t r = 0; r < rast.rows(); ++r) {
				if(rast.get(0, r) == 0)
					rast.floodFill(0, r, 0, 2);
				if(rast.get(rast.cols() - 1, r) == 0)
					rast.floodFill(rast.cols() - 1, r, 0, 2);
			}

			g_debug("cleaning");
			for(uint64_t i = 0; i < rast.size(); ++i)
				rast.set(i, rast.get(i) == 2 ? 0 : 1);

			g_debug("polygonizing");
			Raster<uint16_t> out(dstFile, 1, bounds, res, res, 0);
			out.writeBlock(rast);
			if(!polyFile.empty())
				out.polygonize(polyFile, 1);
		}

	} // las

} // geotools

void usage() {
	std::cerr << "Usage: lasboundary [options] -i <src dir> -o <dst file>\n"
			<< "	This program creates a Shapefile containing the boundary \n"
			<< "  of a point cloud represented by a set of LAS files. The \n"
			<< "  boundary is an alpha shape based on the Delaunay triangulation  \n"
			<< "  with alpha as the square of the given radius.\n"
			<< "  src dir  - The source directory or a single LAS file.\n"
			<< "  dst file - The output Shapefile.\n"
			<< "  -c       - A comma-delimited list of integers indicating \n"
			<< "             which classes to preserve (e.g. 2 = ground). Defaults to all.\n"
			<< "  -r       - The resolution of the grid. Default 1.0m.\n"
			<< "  -s       - The integer EPSG ID of the coordinate reference system.\n";
}

int main(int argc, char **argv) {

	std::string polyFile;
	std::string dstFile;
	int srid = 0;
	double res = 20.0;
	std::set<int> classes;
	std::vector<std::string> srcFiles;

	for (int i = 1; i < argc; ++i) {
		std::string s(argv[i]);
		if (s == "-c") {
			// Gets the set of classes to keep
			Util::intSplit(classes, argv[++i]);
		} else if (s == "-r") {
			res = atof(argv[++i]);
		} else if (s == "-s") {
			srid = atoi(argv[++i]);
		} else if (s == "-o") {
			dstFile = argv[++i];
		} else if (s == "-p") {
			polyFile = argv[++i];
		} else {
			srcFiles.push_back(argv[i]);
		}
	}

	try {

		geotools::las::lasboundary(srcFiles, dstFile, polyFile, srid, res, classes);

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;

}
