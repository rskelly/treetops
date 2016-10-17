/*
 * Grids a point cloud represented by one or more LAS files.
 * Can produce grids of from intensity and elevation, using
 * minimum, maximum, mean, std dev, density, variance and count.
 *
 * Authored by: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <cstdio>
#include <math.h>
#include <exception>
#include <unordered_set>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "pointstats.hpp"
#include "lasutil.hpp"
#include "pointstream.hpp"
#include "cellstats.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::point;
using namespace geotools::las;
using namespace geotools::point::stats;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Projection_traits_xy_3<K>  					Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> 					Delaunay;
typedef K::Point_3 											Point_3;
typedef K::Plane_3 											Plane_3;
typedef Delaunay::Finite_faces_iterator						Finite_faces_iterator;
typedef Delaunay::Face										Face;

namespace geotools {

	namespace point {

		namespace pointstats_config {
			
			double defaultResolution = 10.0; 
			bool defaultSnapToGrid = true;
			unsigned char defaultType = TYPE_MEAN;
			unsigned char defaultAttribute = ATT_HEIGHT;
			unsigned char defaultAngleLimit = 180;
			unsigned char defaultGapFraction = GAP_BLB;
			unsigned int defaultQuantile = 49;
			unsigned int defaultQuantiles = 100;
			unsigned int defaultFilterQuantiles = 100;
			unsigned int defaultFilterQuantileFrom = 0;
			unsigned int defaultFilterQuantileTo = 99;
			unsigned int defaultThreads = 1;

			std::set<unsigned char> defaultClasses = {2}; 

			std::map<std::string, unsigned char> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, 
				{"Median", TYPE_MEDIAN}, {"Rugosity", TYPE_RUGOSITY}, {"Kurtosis", TYPE_KURTOSIS}, {"Skewness", TYPE_SKEW},
				{"Gap Fraction", TYPE_GAP_FRACTION}
			};
			std::map<std::string, unsigned char> attributes = {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};
			std::map<std::string, unsigned char> gapFractionTypes = {
				{"IR", GAP_IR}, {"BLa", GAP_BLA}, {"BLb", GAP_BLB}, {"RR", GAP_RR}, {"FR", GAP_FR}
			};

		} // config
		
		unsigned char PointStatsConfig::parseAtt(const std::string &attStr) {
			if("intensity" == attStr) {
				return ATT_INTENSITY;
			} else if("height" == attStr) {
				return ATT_HEIGHT;
			} 
			return 0;
		}

		unsigned char PointStatsConfig::parseType(const std::string &typeStr) {
			if("min" == typeStr) {
				return TYPE_MIN;
			} else if("max" == typeStr) {
				return TYPE_MAX;
			} else if("mean" == typeStr) {
				return TYPE_MEAN;
			} else if("density" == typeStr) {
				return TYPE_DENSITY;
			} else if("variance" == typeStr) {
				return TYPE_VARIANCE;
			} else if("stddev" == typeStr) {
				return TYPE_STDDEV;
			} else if("pvariance" == typeStr) {
				return TYPE_PVARIANCE;
			} else if("pstddev" == typeStr) {
				return TYPE_PSTDDEV;
			} else if("count" == typeStr) {
				return TYPE_COUNT;
			} else if("median" == typeStr) {
				return TYPE_MEDIAN;
			} else if("skew" == typeStr) {
				return TYPE_SKEW;
			} else if("rugosity" == typeStr) {
				return TYPE_RUGOSITY;
			} else if("kurtosis" == typeStr) {
				return TYPE_KURTOSIS;
			}
			return 0;
		}

		void PointStats::checkConfig(const PointStatsConfig &config) {
			if(config.resolution <= 0.0)
				g_argerr("Resolution must be > 0: " << config.resolution);
			if(config.sourceFiles.size() == 0)
				g_argerr("At least one input file is required.");
			if(config.dstFile.empty()) 
				g_argerr("An output file is required.");
			if(config.attribute == 0)
				g_argerr("An attribute is required.");
			if(config.type == 0)
				g_argerr("A valid type is required.");
			if(config.classes.size() == 0)
				g_warn("No classes given. Matching all classes.");
			if(config.angleLimit <= 0)
				g_argerr("Angle limit must be greater than zero.");

			g_debug("Resolution: " << config.resolution);
			g_debug("Files: " << config.sourceFiles.size());
			g_debug("Destination: " << config.dstFile);
			g_debug("Attribute: " << config.attribute);
			g_debug("Type: " << config.type);
			g_debug("Classes: " << config.classes.size());
			g_debug("Angle Limit: " << config.angleLimit);
		}

		void PointStats::computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds,
			std::set<std::string> &selectedFiles, Bounds &workBounds, unsigned long *pointCount) {
			g_debug(" -- computeWorkBounds - work bounds initial: " << workBounds.print());
			liblas::ReaderFactory rf;
			unsigned long count = 0;
			for(const std::string &file : files) {
				g_debug(" -- computeWorkBounds - checking file " << file);
				geotools::las::PointStream ps(file);
				count += ps.pointCount();
				if(bounds.intersects(ps.fileBounds(), 2)) {
					selectedFiles.insert(file);
					workBounds.extend(ps.fileBounds());
				}
			}
			*pointCount = count;
			g_debug(" -- computeWorkBounds - work bounds final: " << workBounds.print() << "; point count " << *pointCount);
		}

		CellStats* PointStats::getComputer(const PointStatsConfig &config) {
			using namespace geotools::point::stats;
			switch(config.type) {
			case TYPE_MEAN: 		return new CellMean();
			case TYPE_MEDIAN: 		return new CellMedian();
			case TYPE_COUNT: 		return new CellCount();
			case TYPE_STDDEV: 		return new CellSampleStdDev();
			case TYPE_VARIANCE: 	return new CellSampleVariance();
			case TYPE_PSTDDEV: 		return new CellPopulationStdDev();
			case TYPE_PVARIANCE:	return new CellPopulationVariance();
			case TYPE_DENSITY: 		return new CellDensity(g_sq(config.resolution));
			case TYPE_RUGOSITY: 	return new CellRugosity();
			case TYPE_MAX: 			return new CellMax();
			case TYPE_MIN: 			return new CellMin();
			case TYPE_KURTOSIS: 	return new CellKurtosis();
			case TYPE_SKEW: 		return new CellSkewness();
			case TYPE_QUANTILE: 	return new CellQuantile(config.quantile, config.quantiles);
			case TYPE_GAP_FRACTION:	return new CellGapFraction(config.gapFractionType);
			default:
				g_argerr("Invalid statistic type: " << config.type);
			}
		}

		void PointStats::pointstats(const PointStatsConfig &config, const Callbacks *callbacks) {

			checkConfig(config);
			
		 	if(config.threads > 0) {
		 		g_debug(" -- pointstats running with " << config.threads << " threads");
		 		omp_set_dynamic(1);
		 		omp_set_num_threads(config.threads);
		 	} else {
		 		g_argerr("Run with >=1 threads.");
		 	}

			// Compute the work bounds, and store the list
			// of relevant files. Snap the bounds if necessary.
			std::set<std::string> files;
			unsigned long pointCount;
			Bounds workBounds;
			workBounds.collapse();
			computeWorkBounds(config.sourceFiles, config.bounds, files, workBounds, &pointCount);
			if(config.snap) {
				workBounds.snap(config.resolution);
				g_debug(" -- pointstats - snapped work bounds: " << workBounds.print());
			}

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			grid.fill(-9999.0);
			g_debug(" -- pointstats - raster size: " << grid.cols() << ", " << grid.rows());

			std::unique_ptr<CellStats> computer(getComputer(config));
			CellStatsFilter *filter = nullptr, *tmp;
			if(config.hasClasses()) 
				tmp = filter = new ClassFilter(config.classes);
			if(config.hasQuantileFilter())
				tmp = tmp->chain(new QuantileFilter(config.quantileFilter, config.quantileFilterFrom, config.quantileFilterTo));
			if(filter)
				computer->setFilter(filter);

			std::vector<std::string> filesv(files.begin(), files.end());
			unsigned int fileIdx = 0;

			std::map<unsigned long, std::list<std::shared_ptr<LASPoint> > > edges;

			#pragma omp parallel for
			for(unsigned int f = 0; f < filesv.size(); ++f) {

				g_debug(" -- pointstats - file " << filesv[f]);

				if(callbacks)
					callbacks->overallCallback((fileIdx + 0.5f) / filesv.size());

				const std::string &file = filesv[f];
				std::map<unsigned long, std::list<std::shared_ptr<LASPoint> > > values;

				LASPoint pt;
				PointStream ps(file);

				while(ps.next(pt)) {

					unsigned long idx = grid.toRow(pt.y) * grid.cols() + grid.toCol(pt.x);
					std::shared_ptr<LASPoint> pv(new LASPoint(pt));
					values[idx].push_back(std::move(pv));

				}

				g_debug(" -- pointstats - computing results");
				std::list<unsigned long> edgeIds;
				for(const auto &it : values) {
					int col = it.first % grid.cols();
					int row = it.first / grid.cols();
					if(ps.contains(grid.toX(col), grid.toY(row)) 
							&& ps.contains(grid.toX(col + 1), grid.toY(row + 1))) {
						if(filter)
							filter->setPoints(it.second);
						grid.set(it.first, computer->compute(it.second));
					} else {
						edgeIds.push_back(it.first);
					}
				}

				g_debug(" -- pointstats - " << edgeIds.size() << " edge cells deferred.");
				// Write to the edge pixel map to solve later.
				#pragma omp critical
				{
					for(const unsigned long &id : edgeIds) {
						for(std::shared_ptr<LASPoint> &pt : values[id])
							edges[id].push_back(std::move(pt));
					}
	
					++fileIdx; // Do it here to save another critical/atomic
				}

				if(callbacks)
					callbacks->overallCallback((float) fileIdx / filesv.size()); // N.B. fileIdx already incremented

			}

			unsigned int i = 0;
			std::vector<unsigned long> ids(edges.size());
			for(const auto &it : edges)
				ids[i++] = it.first;
			std::sort(ids.begin(), ids.end());

			#pragma omp parallel for
			for(unsigned int i = 0; i < ids.size(); ++i) {
				if(filter)
					filter->setPoints(edges[ids[i]]);
				grid.set(ids[i], computer->compute(edges[ids[i]]));
			}

			if(callbacks)
				callbacks->overallCallback(1.0f);

		}

	} // point

} // geotools

