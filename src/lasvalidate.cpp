/*
 * This tool reads survey points from a CSV file and LiDAR returns
 * from a set of LAS files, in order to compute the difference between
 * the survey elevation and the barycentric elevation produced by a TIN
 * of the LiDAR points.
 *
 * The user can configure the classification of points to use, and the radius 
 * to search for points. The radius should include enough points to produce a triangulation
 * that encompasses the survey point. (If this cannot be acheived, the interpolated
 * elevation is NaN.)
 *
 * The output is a CSV file containing the 3D survey position and the inteprolated elevation.
 * Optionally, the output can include all of the LiDAR returns which were used to calculate
 * the TIN surface.
 *
 * Authored by: Rob Skelly rob@dijital.ca
 */

 #include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include <geos/triangulate/DelaunayTriangulationBuilder.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/Coordinate.h>


#include <liblas/liblas.hpp>

#include "geotools.h"
#include "Util.hpp"

namespace las = liblas;
namespace geom = geos::geom;

/**
 * A simple class for storing LiDAR returns.
 */
class Pnt {
public:
	double x;
	double y;
	double z;
	int cls;
	Pnt(double x, double y, double z, int cls) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->cls = cls;
	}
};

/**
 *  A survey point. Contains the position, the interpolated
 * elevation and a list of related LiDAR returns.
 */
class Sample {
public:
	double x;
	double y;
	double z;
	double interpZ;
	std::vector<Pnt> returns;
	Sample(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->interpZ = nan("");
	}
};

void usage() {
	std::cerr << "Usage: lasvalidate [options] <lasfiles>" << std::endl
			<< "    This program takes a spreadsheed of 3D survey locations and a LiDAR point cloud, in the form " << std::endl
			<< "    of a list of LAS files, and interpolates the elevation of each survey point from the resulting " << std::endl
			<< "    TIN. This allows the survey elevation to be compared to the LiDAR surface elevaton." << std::endl
			<< "    <lasfiles>  Is a list of las files." << std::endl
			<< "    -p  --  A csv file with x, y and z columns. These are the survey locations." << std::endl
			<< "    -o  --  The output file." << std::endl
			<< "    -r  --  The radius to search for lidar returns." << std::endl
			<< "    -c  --  The classes to include; comma-separated list." << std::endl
			<< "    -i  --  Include the points used to calculate the interpolated elevation in the output." << std::endl
			<< "    -v  --  Verbose mode." << std::endl;
}

/**
 * Find the distance between two coordinates.
 */
double dist(double x1, double y1, double x2, double y2) {
	return sqrt(_sq(x2-x1) + _sq(y2-y1));
}

/**
 * Compute a boundary that contains the sample points. Extend it by the
 * search radius. This is used as a naive LAS file filter.
 */
void computeBounds(std::vector<double> &bounds, std::vector<Sample> &points, double distance) {
	for(int i = 0; i < points.size(); ++i) {
		bounds[0] = _min(points[i].x, bounds[0]);
		bounds[1] = _min(points[i].y, bounds[1]);
		bounds[2] = _max(points[i].x, bounds[2]);
		bounds[3] = _max(points[i].y, bounds[3]);
	}	
	bounds[0] -= distance;
	bounds[1] -= distance;
	bounds[2] += distance;
	bounds[3] += distance;
}

/**
 * Load samples (surveys) from the CSV file.
 */
void loadSamples(std::string &datafile, std::vector<Sample> &points) {
	std::vector<std::tuple<double, double, double> > pts;
	Util::loadXYZSamples(datafile, pts);
	for(auto it = pts.begin(); it != pts.end(); ++it) 
		points.push_back(Sample(std::get<0>(*it), std::get<1>(*it), std::get<2>(*it)));
}

/**
 * Write the output file. Optionally include the LiDAR returns.
 */
void writeOutput(std::string &outfile, std::vector<Sample> &samples, bool keepPoints) {
	std::ofstream out;
	out.open(outfile.c_str());
	out << "station_index,station_x,station_y,station_z,station_interp_z";
	if(keepPoints)
		out << ",lidar_x,lidar_y,lidar_z,lidar_class";
	out << std::endl << std::setprecision(9);
	for(int i = 0; i < samples.size(); ++i) {
		Sample samp = samples[i];
		out << i << "," << samp.x << "," << samp.y << "," << samp.z << "," << samp.interpZ;
		if(keepPoints) {
			for(int j = 0; j < samples[i].returns.size(); ++j) {
				Pnt ret = samp.returns[j];
				out << ret.x << "," << ret.y << "," << ret.z << "," << ret.cls;
			}
		}
		out << std::endl;
	}
	out.close();
}

/**
 * Distance between two coordinates (triangle side length.)
 */
double sideLength(const geom::Coordinate &a, const geom::Coordinate &b) {
	return sqrt(_sq(a.x - b.x) + _sq(a.y - b.y));	
}

/**
 * Area of the triangle described by 3 coordinates.
 */
double triArea(const geom::Coordinate &c0, const geom::Coordinate &c1, const geom::Coordinate &c2) {
	double a = sideLength(c0, c1);
	double b = sideLength(c1, c2);
	double c = sideLength(c2, c0);
	double s = (a + b + c) / 2.0;
	return sqrt(s * (s - a) * (s - b) * (s - c));
}

/**
 * Interpolate the value of the triangle at point cs.
 */
double interpolateTriangle(const geom::Coordinate *cs, const geom::Geometry *tri) {
	const geom::Coordinate c0 = tri->getCoordinates()->getAt(0);
	const geom::Coordinate c1 = tri->getCoordinates()->getAt(1);
	const geom::Coordinate c2 = tri->getCoordinates()->getAt(2);
	double tat = triArea(c0, c1, c2);
	double ta2 = triArea(c0, c1, *cs);
	double ta1 = triArea(c0, c2, *cs);
	double ta0 = triArea(c1, c2, *cs);
	return (ta0 / tat) * c0.z + (ta1 / tat) * c1.z + (ta2 / tat) * c2.z;
}

/**
 * Find the interpolated elevation of the sample using a TIN of
 * its nearby LiDAR returns.
 */
void interpolateSampleZ(Sample &sample) {
	using namespace geos::geom;
	using namespace geos::triangulate;
	GeometryFactory gf;
	Coordinate sc(sample.x, sample.y, sample.z);
	// Convert returns to Points.
	std::vector<Geometry *> points;
	for(auto it = sample.returns.begin(); it != sample.returns.end(); ++it)
		points.push_back(gf.createPoint(Coordinate(it->x, it->y, it->z)));
	GeometryCollection *mp = gf.createGeometryCollection(points);
	// Build Delaunay.
	DelaunayTriangulationBuilder dtb;
	dtb.setSites(*mp);
	// Interpolate triangles.
	std::auto_ptr<GeometryCollection> tris = dtb.getTriangles(gf);
	for(size_t i = 0; i < tris->getNumGeometries(); ++i) {
		const Geometry *tri = tris->getGeometryN(i);
		if(tri->contains(gf.createPoint(sc))) {
			const Coordinate *scc = &sc;
			sample.interpZ = interpolateTriangle(scc, tri);
			break;
		}
	}
}

/**
 * Validate the LiDAR point cloud using surveys stored in outfile.
 */
void validate(std::string &outfile, std::string &datafile, std::vector<std::string> &lasfiles, 
	std::set<int> &classes, double distance, bool keepPoints, bool quiet) {

	if(outfile.empty())
		throw "Outfile not given.";
	if(datafile.empty())
		throw "Data file not given.";
	if(lasfiles.size() == 0)
		throw "No las files.";
	if(distance < 0.0)
		throw "Distance must be greater than zero.";

	if(!quiet && classes.size() == 0)
		std::cerr << "No classes given; matching all classes." << std::endl;

	std::vector<Sample> samples;
	loadSamples(datafile, samples);

	if(!quiet)
		std::cerr << samples.size() << " samples found." << std::endl;

	// Set up search bounds for filtering files.
	std::vector<double> bounds = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX };
	computeBounds(bounds, samples, distance);

	if(!quiet)
		std::cerr << "Bounds: " << bounds[0] << "," << bounds[1] << "," << bounds[2] << "," << bounds[3] << std::endl;

	las::ReaderFactory rf;
	for(auto it = lasfiles.begin(); it != lasfiles.end(); ++it) {

		std::ifstream in(it->c_str());
		las::Reader r = rf.CreateWithStream(in);
		las::Header h = r.GetHeader();

		// Check that this file is relevant.
		std::vector<double> bounds0 = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX };
		if(!Util::computeLasBounds(h, bounds0, 2))
			Util::computeLasBounds(r, bounds0, 2); // If the header bounds are bogus.

		if(!quiet)
			std::cerr << "LAS bounds: " << bounds0[0] << "," << bounds0[1] << "," << bounds0[2] << "," << bounds0[3] << std::endl;

		if(!Util::intersects(bounds, bounds0, 2)) {
			if(!quiet)
				std::cerr << "Skipping..." << std::endl;
			continue;
		}

		if(!quiet)
			std::cerr << *it << std::endl;

		while(r.ReadNextPoint()) {
			las::Point pt = r.GetPoint();
			// If this point is not in the class list, skip it.
			int cls = pt.GetClassification().GetClass();
			if(classes.size() > 0 && !Util::inList(classes, cls))
				continue;

			double lx = pt.GetX();
			double ly = pt.GetY();

			for(auto it = samples.begin(); it != samples.end(); ++it) {
				double px = it->x;
				double py = it->y;
				// If outside of the radius, skip it.
				if(dist(px, py, lx, ly) > distance)
					continue;

				double lz = pt.GetZ();
				it->returns.push_back(Pnt(lx, ly, lz, cls));

			}
		}
	}

	std::cerr << "Interpolating..." << std::endl;
	for(auto it = samples.begin(); it != samples.end(); ++it) 
		interpolateSampleZ(*it);
	
	std::cerr << "Writing..." << std::endl;
	writeOutput(outfile, samples, keepPoints);

}

int main(int argc, char **argv) {

 	try {

 		std::vector<std::string> lasfiles;
  		std::string outfile;
  		std::string datafile;
  		std::set<int> classes;
 		double distance = 0.0;
 		bool quiet = true;
 		bool keepPoints = false;

 		for(int i = 1; i < argc; ++i) {
 			std::string p(argv[i]);
 			if(p == "-o") {
 				outfile = argv[++i];
 			} else if(p == "-p") {
 				datafile = argv[++i];
 			} else if(p == "-i") {
 				keepPoints = true;
 			} else if(p == "-r") {
 				distance = atof(argv[++i]);
 			} else if(p == "-v") {
 				quiet = false;
 			} else if(p == "-c") {
				Util::intSplit(classes, argv[++i]);
 			} else {
 				lasfiles.push_back(argv[i]);
 			}
 		}

		validate(outfile, datafile, lasfiles, classes, distance, keepPoints, quiet);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



