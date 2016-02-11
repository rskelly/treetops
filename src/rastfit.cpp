/*
 * rastfit
 *
 * This program adjusts a raster's elevations to match another's using natural
 * neighbour interpolation. It does this by collecting a configurable number of
 * random samples (constrained by an optional mask) and using these as sites in the 
 * triangulation. 
 *
 *  Created on: Jan 16, 2016
 *      Author: rob
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/ch_jarvis.h>

#include "Util.hpp"
#include "Raster.hpp"
#include "ShapeWriter.hpp"

typedef CGAL::Simple_cartesian<float>                                                 K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>                  Vb; // Vertex can store its area
typedef CGAL::Triangulation_data_structure_2<Vb>                                      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                                        Delaunay;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Delaunay>                    AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<Delaunay>    AP;
typedef CGAL::Voronoi_diagram_2<Delaunay,AT,AP>                                       Voronoi;
typedef CGAL::Polygon_2<K>                                                            Polygon_2;
typedef CGAL::Iso_rectangle_2<K>                                                      Iso_rectangle_2;
typedef CGAL::Direction_2<K>                                                          Direction_2;

typedef Delaunay::Vertex_handle             DVertex_handle;
typedef Delaunay::Edge                      DEdge;
typedef Voronoi::Face_handle                VFace_handle;
typedef Voronoi::Face                       VFace;
typedef Voronoi::Ccb_halfedge_circulator    VCcb_halfedge_circulator;
typedef Voronoi::Locate_result              VLocate_result;
typedef Voronoi::Halfedge                   VHalfedge;
typedef Voronoi::Halfedge_handle            VHalfedge_handle;
typedef Voronoi::Vertex_handle              VVertex_handle;
typedef K::Point_2                          Point_2;
typedef K::Ray_2                            Ray_2;
typedef K::Segment_2                        Segment_2;

float _abs(float v) {
	return v < 0.0 ? -v : v;
}

class Point {
public:
	float x, y, diff;
	Point() : x(nan("")), y(nan("")), diff(nan("")) {}
	Point(float x, float y) :
		x(x),
		y(y),
		diff(nan("")) {
	}
	Point(float x, float y, float z) :
		x(x),
		y(y),
		diff(z) {
	}
	Point(Point &p) :
		Point(p.x, p.y, p.diff) {
	}
	bool equals(const Point &p) const {
		return x == p.x && y == p.y;
	}
	~Point() {
	}
};

/**
 * Sorts points so that they're in row/col order, 
 * which improves the read efficiency on the raster.
 */
class samplesort {
public:
	samplesort(Raster<float> &r) : r(r) {
	}
	bool operator()(std::unique_ptr<Point> &a, std::unique_ptr<Point> &b) const {
		int ar = r.toBlockRow(a->y);
		int br = r.toBlockRow(b->y);
		int ac = r.toBlockCol(a->x);
		int bc = r.toBlockCol(b->x);
		if(ar == br) {
			return ac < bc;
		} else {
			return ar < br;
		}
	}
private:
	Raster<float> &r;
};

/**
 * Shuffle and truncate the list of points, adding them to samples.
 */
void shuffle(std::vector<std::unique_ptr<Point> > &samp, std::list<std::unique_ptr<Point> > &samples, int numSamples) {
	std::random_shuffle(samp.begin(), samp.end());
	for(int i = 0; i < numSamples; ++i)
		samples.push_back(std::move(samp[i])); 
}

/**
 * Computes the raster differences for each sample point.
 */
void computeSampleDifferences(Raster<float> &base, Raster<float> &adj, std::list<std::unique_ptr<Point> > &samples) {
	const samplesort ssort(base);
	samples.sort(ssort); // Sort because we want to read the raster efficiently
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		float a = base.get((*pt)->x, (*pt)->y);
		float b = adj.get((*pt)->x, (*pt)->y);
		if(a != base.nodata() && b != adj.nodata()) {
			(*pt)->diff = a - b;
		} else {
			(*pt)->diff = adj.nodata();
		}
	}
}

/**
 * Generates a set of sample points, limited to valid pixels in the mask.
 */
void generateMaskSamples(std::list<std::unique_ptr<Point> > &samples, Raster<char> &mask, int numSamples) {
	std::vector<std::unique_ptr<Point> > pts;
	for(int r = 0; r < mask.rows(); ++r) {
		for(int c = 0; c < mask.cols(); ++c) {
			if((int) mask.get(c, r) != 0) 
				pts.push_back(std::unique_ptr<Point>(new Point(mask.toX(c), mask.toY(r))));
		}
	}
	shuffle(pts, samples, numSamples);
}

float _random() {
	float r = ((float) std::rand()) / RAND_MAX;
	return r;
}

/**
 * Generates a list of random samples within the bounds of a raster.
 */
void generateRandomSamples(std::list<std::unique_ptr<Point> > &samples, Raster<float> &a, Raster<float> &b, 
	float *bounds, int numSamples) {
	do {
		float x = bounds[0] + (bounds[2] - bounds[0]) * _random();
		float y = bounds[1] + (bounds[3] - bounds[1]) * _random();
		if(a.getOrNodata(x, y) != a.nodata() && b.getOrNodata(x, y) != b.nodata()) {
			samples.push_back(std::unique_ptr<Point>(new Point(x, y)));
			--numSamples;
		}
	} while(numSamples);
}

/**
 * Get a segment that represents the halfedge. If the halfedge is
 * finite, return an equivalent segment. If it is not, return a finite segment.
 * Either way, crop to the bounding rectangle.
 */
Segment_2 *getSegment(const VHalfedge &he, const Iso_rectangle_2 &boundary) {

	if(!he.is_ray()) {

		Segment_2 seg(he.source()->point(), he.target()->point());
		return new Segment_2(seg);
		/*
		CGAL::Object inter;
		CGAL::Point_2<K> pi;
		CGAL::Segment_2<K> si;

		// If the segment is inside the rectangle, return it.
		if((boundary.bounded_side(seg.source()) == CGAL::ON_BOUNDED_SIDE && 
			boundary.bounded_side(seg.target()) == CGAL::ON_BOUNDED_SIDE)) {
			return new Segment_2(seg);
		}

		bool bothOutside = (boundary.bounded_side(seg.source()) == CGAL::ON_UNBOUNDED_SIDE && 
			boundary.bounded_side(seg.target()) == CGAL::ON_UNBOUNDED_SIDE);

		// The segment still might intersect the bounding rectangle. If it does,
		// return the clipped segment. It is possible for both ends to be outside
		// the rectangle and for an intersection to occur.
		for(int i = 0; i < 4; ++i) {
			// Rectangle side.
			Segment_2 s2(boundary.vertex(i), boundary.vertex((i + 1) % 4));
			if((inter = CGAL::intersection(seg, s2))) {
				if(CGAL::assign(pi, inter)) {
					// If the intersection is a point, make a segment starting with the
					// point that is inside the boundary. If they're both outside, return null.
					if(boundary.bounded_side(seg.source()) == CGAL::ON_BOUNDED_SIDE) {
						return new Segment_2(seg.source(), pi);
					} else if(boundary.bounded_side(seg.target()) == CGAL::ON_BOUNDED_SIDE) {
						return new Segment_2(pi, seg.target());
					}
				} else if(CGAL::assign(si, inter)) {
					// If the intersection is a segment, just return it.
					return new Segment_2(si);
				}
			}
		}

		//throw "Failed to build segment from bounded edge.";
		*/

	} else {

		// Get the dual of the halfedge and find its midpoint. This
		// helps determine the orientation of the ray.
		DEdge e = he.dual();
		DVertex_handle v1 = e.first->vertex((e.second + 1) % 3);
		DVertex_handle v2 = e.first->vertex((e.second + 2) % 3);

		// Direction perpendicular and to the right of v1, v2.
		Direction_2 dir(v1->point().y() - v2->point().y(), v2->point().x() - v1->point().x());
		// If the ray has a target, reverse the dir.
		Ray_2 ray(he.has_source() ? he.source()->point() : he.target()->point(), he.has_source() ? dir : -dir);

		// Check and return the intersection.
		CGAL::Object inter = CGAL::intersection(boundary, ray);
		CGAL::Point_2<K> p;
		CGAL::Segment_2<K> s;
		if(CGAL::assign(p, inter)) {
			return new Segment_2(ray.source(), p);
		} else if(CGAL::assign(s, inter)) {
			return new Segment_2(s);
		}

		//throw "Failed to build segment from unbounded edge.";
	}

	return nullptr;
}

/**
 * Returns the area of the clipped face.
 */
float faceArea(VFace f, Iso_rectangle_2 &boundary, Voronoi &vor) {
	std::list<Point_2> pts;
	VCcb_halfedge_circulator hc = f.ccb(), done(hc);
	do {
		Segment_2 *seg = getSegment(*hc, boundary);
		if(seg != nullptr) {
			pts.push_back(seg->source());
			pts.push_back(seg->target());
			delete seg;
		}
	} while(++hc != done);

	// Add the boundary vertices that are inside the current face.
	for(int i = 0; i < 4; ++i) {
		VLocate_result lr = vor.locate(boundary.vertex(i));
		VFace_handle *fh = boost::get<VFace_handle>(&lr);
		if(fh && **fh == f)
			pts.push_back(boundary.vertex(i));
	}

	// Build a convex hull of the points, and use the polygon to check the area.
	std::vector<Point_2> hull;
	CGAL::ch_jarvis(pts.begin(), pts.end(), std::back_inserter(hull));
	Polygon_2 poly(hull.begin(), hull.end());
	
	float area = _abs(poly.area());
	return area;

}

void printSamples(std::list<std::unique_ptr<Point> > &samples) {
	std::cout << "x,y,diff" << std::endl;
	for(std::list<std::unique_ptr<Point> >::iterator it = samples.begin(); it != samples.end(); ++it) {
		std::cout << (*it)->x << "," << (*it)->y << "," << (*it)->diff << std::endl;
	}
}

float _min(float a, float b) {
	return a > b ? b : a;
}

float _max(float a, float b) {
	return a < b ? b : a;
}

void computeBounds(Raster<float> &a, Raster<float> &b, float *bounds) {
	bounds[0] = _max(a.toX(0), b.toX(0));
	bounds[1] = _max(a.toY(a.rows()), b.toY(b.rows()));
	bounds[2] = _min(a.toX(a.cols()), b.toX(b.cols()));
	bounds[3] = _min(a.toY(0), b.toY(0));
}

/**
 * Compute an adjustment for adjfile, based on basefile, maskfile and a number of samples. Write the 
 * result to outfile.
 */
void adjust(std::string &basefile, std::string &adjfile, std::string &maskfile, std::string &outfile, int numSamples, float resolution) {

	if(basefile == outfile || adjfile == outfile || maskfile == outfile)
		throw "The output file must not be the same as any input file.";

	if(numSamples < 2)
		throw "Too few samples.";

	if(resolution <= 0.0)
		throw "Invalid resolution.";

	// Initializes source, destination rasters.
	Raster<float> base(basefile, 1, false);
	Raster<float> adj(adjfile, 1, false);
	std::string proj = adj.projection();
	Raster<float> out(outfile, adj.minx(), adj.miny(), adj.maxx(), adj.maxy(), resolution, proj);

	std::list<std::unique_ptr<Point> > samples;

	// Generate a list of samples, shuffle it and clip it to the
	// desired count.
	if(maskfile.empty() || maskfile == "-") {
		float bounds[4];
		computeBounds(base, adj, bounds);
		generateRandomSamples(samples, adj, base, bounds, numSamples);
	} else {
		Raster<char> mask(maskfile, 1);
		generateMaskSamples(samples, mask, numSamples);
	}

	// Compute sample differences.
	computeSampleDifferences(base, adj, samples);

	printSamples(samples);

	// Build a boundary for clipping the voronoi.
	Iso_rectangle_2 boundary(
		Point_2(adj.toX(0) - 1000.0, adj.toY(0) + 1000.0),
		Point_2(adj.toX(adj.cols()) + 1000.0, adj.toY(adj.rows()) - 1000.0)
	);

	// Start a delaunay triangulation.
	Delaunay dt;

	// Maps for site differences and face areas (by vertex ID).
	std::map<unsigned int, float> diffs;
	std::map<unsigned int, float> areas;

	// Build the delaunay triangulation on the samples, and associate the differences
	// with the vertices.
	unsigned int id = 0;
	for(auto pt = samples.begin(); pt != samples.end(); ++pt) {
		Point_2 p((*pt)->x, (*pt)->y);
		DVertex_handle h = dt.insert(p);
		h->info() = id;
		diffs[id] = (*pt)->diff;
		++id;
	}

	// Pre-compute the areas of the original faces.
	Voronoi vt(dt);
	Voronoi::Face_iterator bf = vt.faces_begin();
	do {
		DVertex_handle h = bf->dual();
		float area = faceArea(*bf, boundary, vt);
		areas[h->info()] = area;
	} while(++bf != vt.faces_end());

	// Add the point for the initial cell to the triangulation.
	Point_2 vc(out.toX(0), out.toY(0));
	DVertex_handle vh = dt.insert(vc);

	for(int r = 0; r < out.rows(); ++r) {
		std::cerr << "Row " << r << " of " << out.rows() << std::endl;
		for(int c = 0; c < out.cols(); ++c) {

			Point_2 vc(out.toX(c), out.toY(r));
			dt.move_if_no_collision(vh, vc);

			Voronoi vt0(dt);
			VLocate_result lr = vt0.locate(vh->point());
			VFace_handle fh = boost::get<VFace_handle>(lr);

			float area = faceArea(*fh, boundary, vt0);
			std::cerr << "Area: " << area << std::endl;
			// float z = adj.get(out.toX(c), out.toY(r)), z0 = z;
			// TODO: Not adjusting now; just making a difference raster.
			float z0 = 0.0;
			if(area > 0.0) {
				VCcb_halfedge_circulator hc = vt0.ccb_halfedges(fh), done(hc);
				do {
					VFace_handle f = hc->twin()->face();
					DVertex_handle v = f->dual();
					float a = faceArea(*f, boundary, vt0);
					unsigned int id = v->info();
					std::cerr << "Area " << id << ": " << (areas[id] - a) << ", " << areas[id] << ", " << a << ", " << ((areas[id] - a) / area) << std::endl;
					z0 += ((areas[id] - a) / area) * diffs[id];
				} while(++hc != done);
			} else {
				ShapeWriter sw;
				sw.put(boundary);
				sw.put(dt);
				sw.put(vt, boundary);
				sw.put(vt0, boundary);
				sw.put(vc);
				sw.write();
				std::cerr << "Area is zero at " << c << ", " << r << std::endl;
				throw "Zero area.";
			}
			out.set(c, r, z0);
		}
	}
}

void usage() {
	std::cout << "Usage: rastfit <base file> <file to adjust> <mask file | \"-\"> <output file> <samples> <resolution>" << std::endl;
	std::cout << "	This program will adjust one raster's elevations to match another's " << std::endl;
	std::cout << "	using a natural neighbours interpoled adjustment to all pixels. The mask is " << std::endl;
	std::cout << "  optional and used to limit the placement of samples to where elevations are " << std::endl;
	std::cout << "	known to be good. The output is an adjustment raster with the same extent " << std::endl;
	std::cout << "	as the input raster, with the new resolution. " << std::endl;
}

int main(int argc, char **argv) {

 	try {

 		if(argc < 7)
 			throw "Too few arguments.";

 		std::string basefile = argv[1];
 		std::string adjfile = argv[2];
 		std::string maskfile = argv[3];
  		std::string outfile = argv[4];
 		int samples = atoi(argv[5]);
 		float resolution = atof(argv[6]);

 		adjust(basefile, adjfile, maskfile, outfile, samples, resolution);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }



