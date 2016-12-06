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

#include "lasvalidate.hpp"
#include "util.hpp"
#include "lasutil.hpp"

namespace geom = geos::geom;

using namespace geotools::util;

namespace geotools {

namespace las {

namespace util {

/**
 * Find the distance between two coordinates.
 */
double dist(double x1, double y1, double x2, double y2) {
	return sqrt(g_sq(x2 - x1) + g_sq(y2 - y1));
}

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
	std::string id;
	double x;
	double y;
	double z;
	double interpZ;
	std::vector<Pnt> returns;

	Sample(std::string &id, double x, double y, double z) {
		this->id = id;
		this->x = x;
		this->y = y;
		this->z = z;
		this->interpZ = nan("");
	}

	Pnt* nearest() {
		Pnt *near = nullptr;
		double dst = G_DBL_MAX_POS;
		for (Pnt &p : returns) {
			double dst0 = dist(p.x, p.y, x, y);
			if (dst0 < dst) {
				dst = dst0;
				near = &p;
			}
		}
		return near;
	}
};

/**
 * Compute a boundary that contains the sample points. Extend it by the
 * search radius. This is used as a naive LAS file filter.
 */
void computeBounds(Bounds &bounds, Sample &sample, double distance) {
	bounds.minx(sample.x - distance);
	bounds.miny(sample.y - distance);
	bounds.maxx(sample.x + distance);
	bounds.maxy(sample.y + distance);
}

/**
 * Load samples (surveys) from the CSV file.
 */
void loadSamples(std::string &datafile, std::vector<Sample> &points) {
	std::vector<std::tuple<std::string, double, double, double> > pts;
	Util::loadIDXYZSamples(datafile, pts);
	for (auto it = pts.begin(); it != pts.end(); ++it)
		points.push_back(
				Sample(std::get < 0 > (*it), std::get < 1 > (*it),
						std::get < 2 > (*it), std::get < 3 > (*it)));
}

/**
 * Write the output file. Optionally include the LiDAR returns.
 */
void writeOutput(std::string &outfile, std::string &pointfile,
		std::vector<Sample> &samples) {
	// Open survey output file.
	std::ofstream sout;
	sout.open(outfile.c_str());
	sout
			<< "station_index,station_id,station_x,station_y,station_z,station_interp_z,nearest_x,nearest_y,nearest_z,nearest_dist"
			<< std::endl << std::setprecision(9);
	// Open point output file.
	bool writePoints = !pointfile.empty();
	std::ofstream pout;
	if (writePoints) {
		pout.open(pointfile.c_str());
		pout << "station_index,lidar_x,lidar_y,lidar_z,lidar_class";
		pout << std::endl << std::setprecision(9);
	}
	for (unsigned int i = 0; i < samples.size(); ++i) {
		Sample samp = samples[i];
		Pnt *near = samp.nearest();
		sout << i << "," << samp.id << "," << samp.x << "," << samp.y << ","
				<< samp.z << "," << samp.interpZ;
		// If there is a nearest point, output its coords.
		if (near != nullptr) {
			sout << "," << near->x << "," << near->y << "," << near->z << ","
					<< dist(samp.x, samp.y, near->x, near->y);
		} else {
			sout << ",,,,";
		}
		sout << std::endl;
		if (writePoints) {
			for (unsigned int j = 0; j < samp.returns.size(); ++j) {
				Pnt ret = samp.returns[j];
				pout << i << "," << ret.x << "," << ret.y << "," << ret.z << ","
						<< ret.cls << std::endl;
			}
		}
	}
	// Files close on destruction.
}

/**
 * Distance between two coordinates (triangle side length.)
 */
double sideLength(const geom::Coordinate &a, const geom::Coordinate &b) {
	return sqrt(g_sq(a.x - b.x) + g_sq(a.y - b.y));
}

/**
 * Area of the triangle described by 3 coordinates.
 */
double triArea(const geom::Coordinate &c0, const geom::Coordinate &c1,
		const geom::Coordinate &c2) {

	std::vector<double> sides = { sideLength(c0, c1), sideLength(c1, c2),
			sideLength(c2, c0) };
	std::sort(sides.begin(), sides.end());
	double a = sides[0];
	double b = sides[1];
	double c = sides[2];
	return 0.25
			* sqrt(
					(a + (b + c)) * (c - (a - b)) * (c + (a - b))
							* (a + (b - c)));
}

/**
 * Interpolate the value of the triangle at point cs.
 */
double interpolateTriangle(const geom::Coordinate *cs,
		const geom::Geometry *tri) {
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
	const geom::GeometryFactory *gf =
			geom::GeometryFactory::getDefaultInstance();
	Coordinate sc(sample.x, sample.y, sample.z);
	// Convert returns to Points.
	std::vector<Geometry *> points;
	for (Pnt &pt : sample.returns)
		points.push_back(gf->createPoint(Coordinate(pt.x, pt.y, pt.z)));
	GeometryCollection *mp = gf->createGeometryCollection(points);
	// Build Delaunay.
	DelaunayTriangulationBuilder dtb;
	dtb.setSites(*mp);
	// Interpolate triangles.
	std::unique_ptr<GeometryCollection> tris = dtb.getTriangles(*gf);
	for (size_t i = 0; i < tris->getNumGeometries(); ++i) {
		const Geometry *tri = tris->getGeometryN(i);
		if (tri->contains(gf->createPoint(sc))) {
			const Coordinate *scc = &sc;
			sample.interpZ = interpolateTriangle(scc, tri);
			break;
		}
	}
}

} // util

/**
 * Validate the LiDAR point cloud using surveys stored in outfile.
 */
void validate(std::string &outfile, std::string &pointfile,
		std::string &datafile, std::vector<std::string> &lasfiles,
		std::set<int> &classes, double distance) {

	if (outfile.empty())
		g_argerr("Outfile not given.");
	if (datafile.empty())
		g_argerr("Data file not given.");
	if (lasfiles.size() == 0)
		g_argerr("No LAS files.");
	if (distance < 0.0)
		g_argerr("Distance must be greater than zero.");
	if (outfile == pointfile)
		g_argerr("Output file and LiDAR point file must be different.");
	if (classes.size() == 0)
		g_warn("No classes given; matching all classes.");

	using namespace geotools::las::util;

	std::vector<Sample> samples;
	loadSamples(datafile, samples);

	g_trace(samples.size() << " samples found.");

	liblas::ReaderFactory rf;
	//for(auto it = lasfiles.begin(); it != lasfiles.end(); ++it) {
	for (std::string &lasfile : lasfiles) {

		g_trace(lasfile);

		std::ifstream in(lasfile.c_str());
		liblas::Reader r = rf.CreateWithStream(in);
		liblas::Header h = r.GetHeader();

		// Check that this file is relevant.
		Bounds bounds0;
		bounds0.collapse(2);
		if (!LasUtil::computeLasBounds(h, bounds0, 2))
			LasUtil::computeLasBounds(r, bounds0, 2); // If the header bounds are bogus.

		g_trace(
				"LAS bounds: " << bounds0[0] << "," << bounds0[1] << ","
						<< bounds0[2] << "," << bounds0[3]);

		Bounds bounds;
		bounds.collapse();

		bool inBounds = false;
		for (Sample &sample : samples) {
			computeBounds(bounds, sample, distance);
			if (bounds.intersects(bounds0, 2)) {
				inBounds = true;
				break;
			}
		}

		if (!inBounds) {
			g_trace("No samples in bounds. Skipping...");
			continue;
		}

		while (r.ReadNextPoint()) {
			liblas::Point pt = r.GetPoint();
			// If this point is not in the class list, skip it.
			int cls = pt.GetClassification().GetClass();
			if (classes.size() > 0 && !Util::inList(classes, cls))
				continue;

			double lx = pt.GetX();
			double ly = pt.GetY();

			for (Sample &sample : samples) {
				double px = sample.x;
				double py = sample.y;
				// If outside of the radius, skip it.
				if (dist(px, py, lx, ly) > distance)
					continue;
				// Add to sample returns.
				double lz = pt.GetZ();
				sample.returns.push_back(Pnt(lx, ly, lz, cls));
			}
		}
	}

	g_trace("Interpolating...");
	for (Sample &sample : samples)
		interpolateSampleZ(sample);

	g_trace("Writing...");
	writeOutput(outfile, pointfile, samples);

}

} // las

} // geotools

