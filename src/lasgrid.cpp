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

#include "lasgrid.hpp"
#include "lasutil.hpp"
#include "pointstream.hpp"
#include "simplegeom.hpp"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/Triangle.h"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::las;
using namespace geotools::geom;

namespace geotools {

	namespace las {

		namespace lasgrid_config {
			
			double defaultResolution = 2.0; 
			double defaultRadius = std::sqrt(g_sq(defaultResolution / 2.0) * 2.0); 
			bool defaultSnapToGrid = true;
			unsigned char defaultType = TYPE_MEAN;
			unsigned char defaultAttribute = ATT_HEIGHT;
			unsigned char defaultAngleLimit = 180;
			std::set<unsigned char> defaultClasses = {2}; 
			std::map<std::string, unsigned char> types = {
				{"Minimum", TYPE_MIN}, {"Maximum", TYPE_MAX}, {"Mean", TYPE_MEAN}, {"Density", TYPE_DENSITY},
				{"Sample Variance", TYPE_VARIANCE}, {"Sample Std. Dev.", TYPE_STDDEV}, {"Population Variance", TYPE_PVARIANCE},
				{"Population Std. Dev.", TYPE_PSTDDEV}, {"Count", TYPE_COUNT}, {"Quantile", TYPE_QUANTILE}, 
				{"Median", TYPE_MEDIAN}, {"Rugosity", TYPE_RUGOSITY}, {"Kurtosis", TYPE_KURTOSIS}, {"Skewness", TYPE_SKEW}
			};
			std::map<std::string, unsigned char> attributes = {
				{"Height", ATT_HEIGHT}, {"Intensity", ATT_INTENSITY}
			};

		} // config
		
		namespace lasgrid_util {

			/**
			 * Comparator for sorting doubles.
			 */
			int fcmp(const void * a, const void * b) {
				const double * aa = (const double *) a;
				const double * bb = (const double *) b;
				return (*aa > *bb) - (*bb > *aa);
			}

			/**
			 * Returns true if the point is within the radius associated with a cell's centroid.
			 * @param px The x coordinate of the point.
			 * @param py The y coordinate of the point.
			 * @param col The column of the cell of interest.
			 * @param row The row of the cell of interest.
			 * @param radius The radius around the cell's centroid.
			 * @param resolution The resolution of the output raster.
			 * @param bounds The bounds of the raster.
			 */
			bool inRadius(double px, double py, int col, int row, double radius,
					double resolution, Bounds &bounds){
				if(radius == 0.0) return true;
				// If a radius is given, extract the x and y of the current cell's centroid
				// and measure its distance (squared) from the point.
				double x = col * resolution + bounds[0] + resolution * 0.5;
				double y = row * resolution + bounds[1] + resolution * 0.5;
				// If the cell is outside the radius, ignore it.
				double r = sqrt(g_sq(x - px) + g_sq(y - py));
				return r <= radius;
			}

			CellStats::~CellStats() {}

			class CellDensity : public CellStats {
			private:
				double m_cellArea;
			public:
				CellDensity(double cellArea) :
					m_cellArea(cellArea) {
				}

				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return values.size() / m_cellArea;
				}
			};

			class CellMean : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double sum = 0.0;
					for(const std::unique_ptr<Pt> &v : values)
						sum += v->z;
					return sum / values.size();
				}
			};

			class CellCount : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					return values.size();
				}
			};

			class CellMedian : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					int i = 0;
					std::vector<double> v(values.size());
					for(const std::unique_ptr<Pt> &pt : values)
						v[i++] = pt->z;
					std::sort(v.begin(), v.end());
					unsigned int size = v.size();
					if(size % 2 == 0) {
						return (v[(int) size / 2] + v[(int) size / 2 - 1]) / 2.0;
					} else {
						return v[(int) size / 2];
					}
				}
			};

			class CellMin : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double min = G_DBL_MAX_POS;
					for(const std::unique_ptr<Pt> &v : values) {
						if(v->z < min)
							min = v->z;
					}
					return min;
				}
			};

			class CellMax : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double max = G_DBL_MAX_NEG;
					for(const std::unique_ptr<Pt> &v : values) {
						if(v->z > max)
							max = v->z;
					}
					return max;
				}
			};

			class CellSampleVariance : public CellStats {
			private:
				CellMean m_mean;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::unique_ptr<Pt> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / (values.size() - 1);
				}
			};

			class CellPopulationVariance : public CellStats {
			private:
				CellMean m_mean;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::unique_ptr<Pt> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / values.size();
				}
			};

			class CellSampleStdDev : public CellStats {
			private:
				CellSampleVariance m_variance;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return std::sqrt(m_variance.compute(values));
				}
			};

			class CellPopulationStdDev : public CellStats {
			private:
				CellPopulationVariance m_variance;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return std::sqrt(m_variance.compute(values));
				}
			};

			class CellSkewness : public CellStats {
			private:
				CellMean m_mean;
				CellSampleStdDev m_stdDev;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					// Fisher-Pearson
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::unique_ptr<Pt> &v : values)
						sum += std::pow(v->z - mean, 3.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 3.0);
				}
			};

			class CellKurtosis : public CellStats {
			private:
				CellMean m_mean;
				CellSampleStdDev m_stdDev;
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::unique_ptr<Pt> &v : values)
						sum += std::pow(v->z - mean, 4.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 4.0) - 3.0;
				}
			};

			class CellQuantile : public CellStats {
			private:
				unsigned int m_quantile;
				unsigned int m_quantiles;
			public:
				CellQuantile(unsigned char quantile, unsigned char quantiles) :
					m_quantile(quantile), m_quantiles(quantiles) {
				}
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					return 0;
				}
			};

			class CellRugosity : public CellStats {
			public:
				double compute(const std::list<std::unique_ptr<Pt> > &values) {
					if(values.size() == 0)
						return -9999.0;
					std::list<std::tuple<float, float, double> > coords;
					for(const std::unique_ptr<Pt> &v : values)
						coords.push_back(std::tuple<float, float, double>(v->x, v->y, v->z));
					std::unique_ptr<SFCGAL::TriangulatedSurface> geom = SimpleGeom::getDelaunayTriangles(coords);
					double tarea = 0;
					double parea = 0;
					for(size_t i = 0; i < geom->numGeometries(); ++i) {
						const SFCGAL::Triangle &t = geom->geometryN(i);
						tarea += SFCGAL::algorithm::area3D(t);
						parea += SFCGAL::algorithm::area(t);
					}
					return tarea / parea;
				}
			};

		} // util

		unsigned char LASGridConfig::parseAtt(const std::string &attStr) {
			if("intensity" == attStr) {
				return ATT_INTENSITY;
			} else if("height" == attStr) {
				return ATT_HEIGHT;
			} 
			return 0;
		}

		unsigned char LASGridConfig::parseType(const std::string &typeStr) {
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

		void LASGrid::setCallbacks(Callbacks *callbacks) {
			m_callbacks.reset(callbacks);
		}

		void LASGrid::checkConfig(const LASGridConfig &config) {
			if(config.resolution <= 0.0)
				g_argerr("Resolution must be > 0: " << config.resolution);
			if(config.radius <= 0.0)
				g_argerr("Radius invalid: " << config.radius);
			if(config.lasFiles.size() == 0)
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

			g_debug("Radius: " << config.radius);
			g_debug("Resolution: " << config.resolution);
			g_debug("Files: " << config.lasFiles.size());
			g_debug("Destination: " << config.dstFile);
			g_debug("Attribute: " << config.attribute);
			g_debug("Type: " << config.type);
			g_debug("Classes: " << config.classes.size());
			g_debug("Angle Limit: " << config.angleLimit);
		}

		void LASGrid::computeWorkBounds(const std::list<std::string> &files, const Bounds &bounds,
			std::set<std::string> &selectedFiles, Bounds &workBounds, unsigned long *pointCount) {
			g_debug(" -- computeWorkBounds - work bounds initial: " << workBounds.print());
			liblas::ReaderFactory rf;
			unsigned long count = 0;
			for(const std::string &file : files) {
				g_debug(" -- computeWorkBounds - checking file " << file);
				PointStream ps(file);
				count += ps.pointCount();
				if(bounds.intersects(ps.fileBounds(), 2)) {
					selectedFiles.insert(file);
					workBounds.extend(ps.fileBounds());
				}
			}
			*pointCount = count;
			g_debug(" -- computeWorkBounds - work bounds final: " << workBounds.print() << "; point count " << *pointCount);
		}

		geotools::las::lasgrid_util::CellStats* LASGrid::getComputer(const LASGridConfig &config) {
			using namespace geotools::las::lasgrid_util;
			switch(config.type) {
			case TYPE_MEAN: return new CellMean();
			case TYPE_MEDIAN: return new CellMedian();
			case TYPE_COUNT: return new CellCount();
			case TYPE_STDDEV: return new CellSampleStdDev();
			case TYPE_VARIANCE: return new CellSampleVariance();
			case TYPE_PSTDDEV: return new CellPopulationStdDev();
			case TYPE_PVARIANCE: return new CellPopulationVariance();
			case TYPE_DENSITY: return new CellDensity(g_sq(config.resolution));
			case TYPE_RUGOSITY: return new CellRugosity();
			case TYPE_MAX: return new CellMax();
			case TYPE_MIN: return new CellMin();
			case TYPE_KURTOSIS: return new CellKurtosis();
			case TYPE_SKEW: return new CellSkewness();
			case TYPE_QUANTILE: return new CellQuantile(config.quantile, config.quantiles);
			default:
				g_argerr("Invalid statistic type: " << config.type);
			}
		}

		void LASGrid::lasgrid(const LASGridConfig &config) {

			checkConfig(config);
			
			using namespace geotools::las::lasgrid_util;

			// Compute the work bounds, and store the list
			// of relevant files. Snap the bounds if necessary.
			std::set<std::string> files;
			unsigned long pointCount;
			Bounds workBounds;
			workBounds.collapse();
			computeWorkBounds(config.lasFiles, config.bounds, files, workBounds, &pointCount);
			if(config.snap) {
				workBounds.snap(config.resolution);
				g_debug(" -- lasgrid - snapped work bounds: " << workBounds.print());
			}

			// Prepare the grid
			// TODO: Only works with UTM north.
			Raster<float> grid(config.dstFile, 1, workBounds, config.resolution, 
				-config.resolution, -9999, config.hsrid);
			grid.fill(-9999.0);

			g_debug(" -- lasgrid - raster size: " << grid.cols() << ", " << grid.rows());

			std::vector<std::string> filesv(files.size());
			filesv.assign(files.begin(), files.end());

			g_debug(" -- lasgrid - " << filesv.size() << " files");

			CellStats *computer = getComputer(config);

			#pragma omp parallel
			{
			
				#pragma omp for
				for(unsigned int i = 0; i < filesv.size(); ++i) {
					const std::string &file = filesv[i];

					std::unordered_map<unsigned long, std::list<std::unique_ptr<Pt> > > values;

					LASPoint pt;
					PointStream ps(file);

					while(ps.next(pt)) {
						
						unsigned long idx = grid.toRow(pt.y) * grid.cols() + grid.toCol(pt.x);
						float x = (float) pt.x;
						float y = (float) pt.y;
						double val = 0;

						switch(config.attribute) {
						case ATT_INTENSITY: 
							val = (double) pt.intensity;
							break;
						default: 
							val = pt.z;
							break;
						}

						std::unique_ptr<Pt> pv(new Pt(x, y, val));
						values[idx].push_back(std::move(pv));
					}

					std::set<unsigned long> ids;
					for(const auto &it : values)
						ids.insert(it.first);

					g_debug(" -- lasgrid - computing results");
					for(const unsigned long &id : ids)
						grid.set(id, computer->compute(values[id]));
				}


			}

			delete computer;


//				if(m_callbacks)
//					m_callbacks->overallCallback((curFile + 0.5f) / files.size());

//						if(m_callbacks) 
//							m_callbacks->fileCallback((curPt + 1.0f) / pts.pointCount());
//			if(m_callbacks)
//				m_callbacks->overallCallback(1.0f);

		}

	} // las

} // geotools

