#include <atomic>
#include <thread>
#include <vector>

#include <boost/filesystem.hpp>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <liblas/liblas.hpp>

#include "geotools.hpp"
#include "pointnormalize.hpp"
#include "raster.hpp"
#include "util.hpp"
#include "lasreader.hpp"

using namespace geotools::point;
using namespace geotools::raster;
using namespace geotools::las;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::All_faces_iterator All_faces_iterator;
typedef Delaunay::Face Face;
typedef Delaunay::Face_handle Face_handle;

PointNormalizeConfig::PointNormalizeConfig() :
		dropNegative(true), dropGround(true), threads(1), overwrite(false), buffer(
				10.0) {
}

class FileSorter {
private:
	double m_colSize;
	double m_rowSize;
public:

	FileSorter(double colSize, double rowSize) :
			m_colSize(colSize), m_rowSize(rowSize) {
	}

	bool operator()(const std::string &a, const std::string &b) {
		LASReader ar(a);
		LASReader br(b);
		Bounds ab = ar.bounds();
		Bounds bb = br.bounds();
		int idxa = ((int) (ab.miny() / m_rowSize))
				* ((int) (ab.width() / m_colSize))
				+ ((int) (ab.minx() / m_colSize));
		int idxb = ((int) (bb.miny() / m_rowSize))
				* ((int) (bb.width() / m_colSize))
				+ ((int) (bb.minx() / m_colSize));
		return idxa < idxb;
	}
};

class FileSearcher {
private:
	std::vector<std::string> m_files;
	std::vector<std::pair<std::string, Bounds> > m_bounds;
public:
	FileSearcher(const std::vector<std::string> &files) :
			m_files(files) {
		init();
	}

	void init() {
		for (const std::string &file : m_files) {
			LASReader r(file);
			Bounds b = r.bounds();
			m_bounds.push_back(std::make_pair(file, b));
		}
	}

	std::vector<std::string> getFiles(const Bounds &bounds) {
		std::vector<std::string> files;
		for (const auto &pr : m_bounds) {
			if (pr.second.intersects(bounds))
				files.push_back(pr.first);
		}
		return files;
	}
};

bool _cancel = false;

void _setup(const PointNormalizeConfig &config) {
	// Sanity checks.
	if (config.outputDir.empty())
		g_argerr("A point output directory must be provided.");
	if (!Util::mkdir(config.outputDir))
		g_argerr("Couldn't create output directory: " << config.outputDir);
	uint32_t threads = config.threads;
	uint32_t tmax = std::thread::hardware_concurrency();
	if (threads < 1) {
		g_warn("Too few threads specified. Using 1.");
		threads = 1;
	}
	if (threads > tmax && tmax > 0) {
		g_warn("Too many threads specified. Using 1.");
		threads = 1;
	}
	// Set the thread count.
	omp_set_dynamic(1);
	omp_set_num_threads(threads);
}

void PointNormalize::normalize(const PointNormalizeConfig &config,
		const Callbacks *callbacks, bool *cancel) {

	using namespace boost::filesystem;

	// If the cancel flag isn't given, take a reference to a dummy
	// that is always false.
	if (!cancel)
		cancel = &_cancel;

	// Set up
	_setup(config);

	// Prepare some necessary variables
	path outDir(config.outputDir);
	liblas::ReaderFactory rf;
	FileSearcher fileSearch(config.sourceFiles);
	uint32_t numFiles = config.sourceFiles.size();
	uint32_t curFile = 0;

	for (const std::string &pointFile : config.sourceFiles) {
		if (*cancel)
			return;
		if (callbacks)
			callbacks->overallCallback(((float) ++curFile - 0.5f) / numFiles);

		// Original and new file paths
		path oldPath(pointFile);
		path newPath(outDir / oldPath.filename());
		if (!config.overwrite && exists(newPath)) {
			g_warn("The output file " << newPath << " already exists.");
			continue;
		}
		g_debug("Normalize: processing " << pointFile << " to " << newPath);

		std::vector<LASPoint> objPts;
		std::list<Point_3> meshPts;

		// Build a mesh from the buffered ground points.
		LASReader lr(pointFile);
		Bounds bounds = lr.bounds();
		bounds.extend(bounds.minx() - config.buffer,
				bounds.miny() - config.buffer);
		bounds.extend(bounds.maxx() + config.buffer,
				bounds.maxy() + config.buffer);
		std::vector<std::string> bufFiles = fileSearch.getFiles(bounds);

		for (const std::string &file : bufFiles) {
			if (*cancel)
				return;
			LASReader r(file);
			LASPoint pt;
			bool isPtFile = file == pointFile;
			while (r.next(pt)) {
				if (*cancel)
					return;
				if (pt.cls == 2 && bounds.contains(pt.x, pt.y)) {
					Point_3 p(pt.x, pt.y, pt.z);
					meshPts.push_back(std::move(p));
				}
				if (isPtFile && (!config.dropGround || pt.cls != 2))
					objPts.push_back(pt);
			}
		}

		// Build the Delaunay triangulation.
		Delaunay dt(meshPts.begin(), meshPts.end());
		meshPts.clear();
		if (*cancel)
			return;

		// Prepare point counts.
		uint64_t finalPtCount = 0;
		uint64_t finalPtCountByReturn[255];
		for (int i = 0; i < 5; ++i)
			finalPtCountByReturn[i] = 0;

		// Final bounds for the new file
		double lbounds[6];
		for (int i = 0; i < 3; ++i) {
			lbounds[i] = G_DBL_MAX_POS;
			lbounds[i + 3] = G_DBL_MAX_NEG;
		}

		// Hit for the locate function
		Face_handle hint;

		std::ofstream outstr(newPath.c_str(), std::ios::binary);
		std::ifstream instr(pointFile.c_str(), std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header outHeader(lasReader.GetHeader());
		liblas::Writer lasWriter(outstr, outHeader);
		instr.close();

		uint64_t ptCount = objPts.size();
		std::atomic<uint64_t> curPt(0);

		#pragma omp parallel for
		for (uint64_t i = 0; i < ptCount; ++i) {
			if (*cancel)
				continue;

			if (callbacks)
				callbacks->stepCallback((float) ++curPt / ptCount);

			LASPoint &pt = objPts[i];
			if (pt.cls == 2) {
				pt.z = 0.0;
			} else {
				Point_3 p(pt.x, pt.y, pt.z);
				hint = dt.locate(p, hint);
				if (hint == NULL || dt.is_infinite(hint)) {
					g_warn("Not found in mesh: " << pt.x << ", " << pt.y);
					continue;
				}
				double area = 0.0;
				double total = 0.0;
				for (int i = 0; i < 3; ++i) {
					Point_3 p1 = hint->vertex(i)->point();
					Point_3 p2 = hint->vertex((i + 1) % 3)->point();
					Point_3 p3 = hint->vertex((i + 2) % 3)->point();
					double h = Util::computeArea(p1.x(), p1.y(), p1.z(), p2.x(),
							p2.y(), p2.z(), p.x(), p.y(), p.z());
					area += h;
					total += h * p3.z();
				}
				if (std::isnan(area)) {
					//g_warn("Area is nan: " << area << "; " << total);
					continue;
				}
				double z = p.z() - total / area;
				if (z < 0.0 && config.dropNegative)
					continue;
				pt.z = z;
			}
			#pragma omp critical(__las_write)
			{
				lbounds[0] = g_min(lbounds[0], pt.x);
				lbounds[1] = g_min(lbounds[1], pt.y);
				lbounds[2] = g_min(lbounds[2], pt.z);
				lbounds[3] = g_max(lbounds[3], pt.x);
				lbounds[4] = g_max(lbounds[4], pt.y);
				lbounds[5] = g_max(lbounds[5], pt.z);
				liblas::Point opt(&outHeader);
				pt.toLibLAS(opt);
				lasWriter.WritePoint(opt);
				finalPtCountByReturn[pt.returnNum - 1]++;
				++finalPtCount;
			}
		}

		outHeader.SetPointRecordsCount(finalPtCount);
		for (int i = 0; i < 5; ++i)
			outHeader.SetPointRecordsByReturnCount(i + 1,
					finalPtCountByReturn[i]);
		outHeader.SetMin(lbounds[0], lbounds[1], lbounds[2]);
		outHeader.SetMax(lbounds[3], lbounds[4], lbounds[5]);
		lasWriter.WriteHeader();

		outstr.close();

		if (callbacks)
			callbacks->overallCallback((float) curFile / numFiles);
	}

}
