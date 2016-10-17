#ifndef __CELLSTATS_HPP__
#define __CELLSTATS_HPP__

#include <list>
#include <memory>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>


#include "pointstream.hpp"

#define GAP_IR 1
#define GAP_BLA 2
#define GAP_BLB 3
#define GAP_RR 4
#define GAP_FR 5

using namespace geotools::las;


typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Projection_traits_xy_3<K>  					Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> 					Delaunay;
typedef K::Point_3 											Point_3;
typedef K::Plane_3 											Plane_3;
typedef Delaunay::Finite_faces_iterator						Finite_faces_iterator;
typedef Delaunay::Face										Face;



namespace geotools {
	
	namespace point {

		namespace stats {

			class CellStatsFilter {
			protected:
				std::list<std::shared_ptr<LASPoint> > m_points;
				CellStatsFilter *m_chain;

				virtual bool keepImpl(const LASPoint &) const =0;
				virtual void init() {}

			public:
				CellStatsFilter() : m_chain(nullptr) {}

				CellStatsFilter* chain(CellStatsFilter *next) {
					m_chain = next;
					return next;
				}

				void setPoints(const std::list<std::shared_ptr<LASPoint> > &points) {
					m_points.clear();
					m_points.assign(points.begin(), points.end());
					init();
					if(m_chain)
						m_chain->setPoints(points);
				}

				virtual bool keep(const LASPoint &pt) const {
					return  keepImpl(pt) && (m_chain && m_chain->keep(pt));
				}

				virtual ~CellStatsFilter() {}
			};

			class ClassFilter : public CellStatsFilter {
			private:
				std::set<unsigned char> m_classes;

			protected:
				bool keepImpl(const LASPoint &pt) const {
					return m_classes.size() == 0 || (m_classes.find(pt.cls) != m_classes.end());
				}

			public:
				ClassFilter(const std::set<unsigned char> &classes) {
					m_classes.insert(classes.begin(), classes.end());
				}

				~ClassFilter() {}
			};

			class QuantileFilter : public CellStatsFilter {
			private:
				int m_quantiles;
				int m_from;
				int m_to;
				double m_min;
				double m_max;

			protected:
				bool keepImpl(const LASPoint &pt) const {
					return pt.z > m_min && pt.z <= m_max;
				}

			public:
				QuantileFilter(int quantiles, int from, int to) :
					m_quantiles(quantiles), m_from(from), m_to(to) {
				}

				void init() {
					unsigned int span = m_points.size() / m_quantiles;
					unsigned int start = m_from * span;
					unsigned int end = m_to * span + 1;
					auto it = m_points.begin(); //std::advance(m_points.begin(), start);
					if(it == m_points.end())
						g_argerr("Quantile start index out of bounds.");
					m_min = (*it)->z;
					//it = std::advance(it, end - start);
					if(it == m_points.end())
						g_argerr("Quantile end index out of bounds.");
					m_max = (*it)->z;
					g_debug(" -- quant init start: " << start << "; end: " << end << "; min: " << m_min << "; max: " << m_max);
					std::stringstream ss;
					for(const std::shared_ptr<LASPoint> &pt : m_points)
						ss << pt->z;
					g_debug(" --- points " << ss.str());
				}

				~QuantileFilter() {}
			};

			class CellStats {
			private:
				const CellStatsFilter *m_filter;

			public:
				CellStats() : m_filter(nullptr) {}

				void setFilter(const CellStatsFilter *filter) {
					m_filter = filter;
				}

				std::list<std::shared_ptr<LASPoint> > filtered(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(m_filter) {
						std::list<std::shared_ptr<LASPoint> > out;
						for(const std::shared_ptr<LASPoint> &pt : values) {
							if(m_filter->keep(*pt))
								out.push_back(std::move(pt));
						}
						return std::move(out);
					} else {
						return std::move(values);
					}				
				}

				virtual double compute(const std::list<std::shared_ptr<LASPoint> >&) =0;
				~CellStats() {
					if(m_filter)
						delete m_filter;
				}
			};

			class CellDensity : public CellStats {
			private:
				double m_cellArea;
			public:
				CellDensity(double cellArea) : CellStats(), 
					m_cellArea(cellArea) {
				}

				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return filtered(values).size() / m_cellArea;
				}
			};

			class CellMean : public CellStats {
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double sum = 0.0;
					for(const std::shared_ptr<LASPoint> &v : filtered(values))
						sum += v->z;
					return sum / values.size();
				}
			};

			class CellCount : public CellStats {
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					return values.size();
				}
			};

			class CellMedian : public CellStats {
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					int i = 0;
					std::vector<double> v(values.size());
					for(const std::shared_ptr<LASPoint> &pt : values)
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
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double min = G_DBL_MAX_POS;
					for(const std::shared_ptr<LASPoint> &v : values) {
						if(v->z < min)
							min = v->z;
					}
					return min;
				}
			};

			class CellMax : public CellStats {
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;
					double max = G_DBL_MAX_NEG;
					for(const std::shared_ptr<LASPoint> &v : values) {
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
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::shared_ptr<LASPoint> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / (values.size() - 1);
				}
			};

			class CellPopulationVariance : public CellStats {
			private:
				CellMean m_mean;
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() <= 1)
						return -9999.0;			
					double mean = m_mean.compute(values);
					double sum = 0;
					for(const std::shared_ptr<LASPoint> &v : values)
						sum += g_sq(g_abs(v->z - mean));
					return sum / values.size();
				}
			};

			class CellSampleStdDev : public CellStats {
			private:
				CellSampleVariance m_variance;
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					return std::sqrt(m_variance.compute(values));
				}
			};

			class CellPopulationStdDev : public CellStats {
			private:
				CellPopulationVariance m_variance;
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
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
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					// Fisher-Pearson
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::shared_ptr<LASPoint> &v : values)
						sum += std::pow(v->z - mean, 3.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 3.0);
				}
			};

			class CellKurtosis : public CellStats {
			private:
				CellMean m_mean;
				CellSampleStdDev m_stdDev;
			public:
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					if(values.size() == 0)
						return -9999.0;
					double mean = m_mean.compute(values);
					double sum = 0.0;
					unsigned int count = values.size();
					for(const std::shared_ptr<LASPoint> &v : values)
						sum += std::pow(v->z - mean, 4.0) / count;
					return sum / std::pow(m_stdDev.compute(values), 4.0) - 3.0;
				}
			};

			class CellQuantile : public CellStats {
			private:
				unsigned int m_quantile;
				unsigned int m_quantiles;
			public:
				CellQuantile(unsigned char quantile, unsigned char quantiles) : CellStats(), 
					m_quantile(quantile), m_quantiles(quantiles) {
				}
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					return 0;
				}
			};


			class CellRugosity : public CellStats { // TODO: Seems to be density-dependent
			private:
				static const int a = 2, b = 1, c = 0;

				double computeArea(const Face &face) {
					Point_3 p1 = face.vertex(0)->point();
					Point_3 p2 = face.vertex(1)->point();
					Point_3 p3 = face.vertex(2)->point();
					double sides[] {
						std::sqrt(std::pow(p1.x() - p2.x(), 2.0) + std::pow(p1.y() - p2.y(), 2.0) + std::pow(p1.z() - p2.z(), 2.0)),
						std::sqrt(std::pow(p2.x() - p3.x(), 2.0) + std::pow(p2.y() - p3.y(), 2.0) + std::pow(p2.z() - p3.z(), 2.0)),
						std::sqrt(std::pow(p3.x() - p1.x(), 2.0) + std::pow(p3.y() - p1.y(), 2.0) + std::pow(p3.z() - p1.z(), 2.0)),
					};
					double s = (sides[a] + sides[b] + sides[c]) / 2.0;
					return std::sqrt(s * (s - sides[a]) * (s - sides[b]) * (s - sides[c]));
				}

				double toPlane(const Point_3 &p, const Plane_3 &plane, const Point_3 &centroid) {
					return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
				}

			public:

				/**
				 * Using Du Preez, 2014 - Arc-Chord Ratio (ACR) Index.
				 */
				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {

					if(values.size() == 0)
						return -9999.0;
					
					std::list<Point_3> pts;
					for(const std::shared_ptr<LASPoint> &v : values)
						pts.push_back(Point_3(v->x, v->y, v->z));
					
					// Delaunay 3D surface area.
					double tarea = 0.0;
					Delaunay dt(pts.begin(), pts.end());
					for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
						tarea += computeArea(*it);

					// Convex hull and POBF
					std::list<Point_3> hull;
					Plane_3 plane;
					Point_3 centroid;
					CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), Gt());
					CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

					// POBF surface area.
					double parea = 0.0;
					std::list<Point_3> poly; // TODO: Is there a faster way?
					for(const Point_3 &pt : hull)
						poly.push_back(Point_3(pt.x(), pt.y(), toPlane(pt, plane, centroid)));
					Delaunay dt2(poly.begin(), poly.end());
					for(Finite_faces_iterator it = dt2.finite_faces_begin(); it != dt2.finite_faces_end(); ++it)
						parea += computeArea(*it);

					return tarea / parea;
				}
			};

			/**
			 * Adapted from:
			 * Hopkinson, C., & Chasmer, L. (2009). Testing LiDAR models of fractional cover across 
			 *	multiple forest ecozones. Remote Sensing of Environment, 113(1), 275–288. 
			 *	http://doi.org/10.1016/j.rse.2008.09.012
			 */
			class CellGapFraction : public CellStats {
			private:
				unsigned char m_type;

				double fcLidarBLa(const std::list<std::shared_ptr<LASPoint> > &values) {
					double gnd = 0.0, all = 0.0;
					for(const std::shared_ptr<LASPoint> &pt : values) {
						if(pt->ground())
							gnd += pt->intensity;
						all += pt->intensity; // TODO: This should perhaps be filtered by class to remove bogus points.
					}
					g_debug(" -- fcLidarBLa " << gnd << ", " << all << ", " << values.size());
					return all != 0.0 ? 1.0 - std::sqrt(gnd / all) : -9999.0;
				}

				double fcLidarBLb(const std::list<std::shared_ptr<LASPoint> > &values) {
					double gndSingle = 0.0, gndLast = 0.0, first = 0.0, single = 0.0, intermediate = 0.0, last = 0.0, total = 0.0;
					for(const std::shared_ptr<LASPoint> &pt : values) {
						if(pt->ground()) {
							if(pt->single())
								gndSingle += pt->intensity;
							if(pt->last())
								gndLast += pt->intensity;
						}
						if(pt->first())
							first += pt->intensity;
						if(pt->single())
							single += pt->intensity;
						if(pt->intermediate())
							intermediate += pt->intensity;
						if(pt->last())
							last += pt->intensity;
						total += pt->intensity; // TODO: This should perhaps be filtered by class to remove bogus points.
					}
					if(total == 0.0) return -9999.0;
					double denom = (first + single) / total + std::sqrt((intermediate + last) / total);
					if(denom == 0.0) return -9999.;
					return (gndSingle / total + std::sqrt(gndLast / total)) / denom;
				}

				double fcLidarIR(const std::list<std::shared_ptr<LASPoint> > &values) {
					double canopy = 0.0, total = 0.0;
					for(const std::shared_ptr<LASPoint> &pt : values) {
						if(!pt->ground())
							canopy += pt->intensity;
						total += pt->intensity;
					}
					return total != 0.0 ? canopy / total : -9999.0;
				}

				
				double fcLidarRR(const std::list<std::shared_ptr<LASPoint> > &values) {
					unsigned int canopy = 0, total = 0;
					for(const std::shared_ptr<LASPoint> &pt : values) {
						if(!pt->ground())
							++canopy;
						++total;
					}
					return total != 0.0 ? (double) canopy / total : -9999.0;
				}
				
				double fcLidarFR(const std::list<std::shared_ptr<LASPoint> > &values) {
					unsigned int canopy = 0, total = 0;
					for(const std::shared_ptr<LASPoint> &pt : values) {
						if(pt->first()) {
							if(!pt->ground())
								++canopy;
							++total;
						}
					}
					return total != 0.0 ? (double) canopy / total : -9999.0;
				}
			public:
				const static unsigned char IR = GAP_IR;
				const static unsigned char BLA = GAP_BLA;
				const static unsigned char BLB = GAP_BLB;
				const static unsigned char RR = GAP_RR;
				const static unsigned char FR = GAP_FR;

				CellGapFraction(unsigned char type) : CellStats(), 
					m_type(type) {
				}

				double compute(const std::list<std::shared_ptr<LASPoint> > &values) {
					switch(m_type) {
					case BLA:	return fcLidarBLa(values);
					case BLB:	return fcLidarBLb(values);
					case IR:	return fcLidarIR(values);
					case RR:	return fcLidarRR(values);
					case FR:	return fcLidarFR(values);
					default:
						g_argerr("Unknown Gap Fraction method: " << m_type);
					}
				}
			};

		}
	}
}

#endif