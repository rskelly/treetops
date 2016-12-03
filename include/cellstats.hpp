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

#include "laspoint.hpp"

#define GAP_IR 1
#define GAP_BLA 2
#define GAP_BLB 3
#define GAP_RR 4
#define GAP_FR 5
#define GAP_CCF 6
#define GAP_GAP 7

using namespace geotools::las;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Face Face;

namespace geotools {

    namespace point {

        namespace stats {

            class CellStatsFilter {
            protected:
                CellStatsFilter *m_chain;
                const std::list<LASPoint*> *m_points;

                virtual bool keepImpl(const LASPoint *) const = 0;

                virtual void init() {
                }

            public:

                CellStatsFilter() : m_chain(nullptr) {
                }

                CellStatsFilter* chain(CellStatsFilter *next) {
                    m_chain = next;
                    return next;
                }

                void setPoints(const std::list<LASPoint*> *points) {
                    m_points = points;
                    init();
                    if (m_chain)
                        m_chain->setPoints(points);
                }

                bool keep(const LASPoint *pt) const {
                    return keepImpl(pt) && (!m_chain || m_chain->keep(pt));
                }

                virtual ~CellStatsFilter() {
                }
            };

            class ClassFilter : public CellStatsFilter {
            private:
                std::set<unsigned char> m_classes;

            protected:

                bool keepImpl(const LASPoint *pt) const {
                    return m_classes.size() == 0 || (m_classes.find(pt->cls) != m_classes.end());
                }

            public:

                ClassFilter(const std::set<unsigned char> &classes) {
                    m_classes.insert(classes.begin(), classes.end());
                }

                ~ClassFilter() {
                }
            };

            class QuantileFilter : public CellStatsFilter {
            private:
                int m_quantiles;
                int m_from;
                int m_to;
                double m_min;
                double m_max;

            protected:

                bool keepImpl(const LASPoint *pt) const {
                    return pt->z > m_min && pt->z <= m_max;
                }

            public:

                QuantileFilter(int quantiles, int from, int to) :
                m_quantiles(quantiles), m_from(from), m_to(to) {
                }

                void init() {
                    auto it = m_points->begin(); //std::advance(m_points.begin(), start);
                    if (it == m_points->end())
                        g_argerr("Quantile start index out of bounds.");
                    m_min = (*it)->z;
                    if (it == m_points->end())
                        g_argerr("Quantile end index out of bounds.");
                    m_max = (*it)->z;
                }

                ~QuantileFilter() {
                }
            };

            class CellStats {
            private:
                CellStatsFilter *m_filter;

            public:

                CellStats() : m_filter(nullptr) {
                }

                void setFilter(CellStatsFilter *filter) {
                    m_filter = filter;
                }

                std::list<LASPoint*> filtered(const std::list<LASPoint*> &values) {
                    std::list<LASPoint*> filtered;
                    if (m_filter) {
                        m_filter->setPoints(&values);
                        for (LASPoint *pt : values) {
                            if (m_filter->keep(pt))
                                filtered.push_back(pt);
                        }
                    } else {
                        filtered.assign(values.begin(), values.end());
                    }
                    return filtered;
                }

                virtual void compute(const std::list<LASPoint*>&, double*) = 0;

                virtual int bands() const {
                    return 1;
                }
                
                ~CellStats() {
                }
            };

            class CellDensity : public CellStats {
            private:
                double m_cellArea;
            public:

                CellDensity(double area = 0.0) :
                    CellStats(),
                    m_cellArea(area) {
                }

                void setArea(double area) {
                    m_cellArea = area;
                }

                double area() {
                    return m_cellArea;
                }

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    g_debug(" -- density filtered " << filt.size());
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        result[0] = filt.size() / m_cellArea;
                    }
                }
            };

            class CellMean : public CellStats {
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        //g_debug(" -- mean: " << filt.size());
                        double sum = 0.0;
                        for (const LASPoint *v : filt)
                            sum += v->z;
                        result[0] =  sum / filt.size();
                    }
                }
            };

            class CellCount : public CellStats {
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    result[0] = filtered(values).size();
                }
            };

            class CellMedian : public CellStats {
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        int i = 0;
                        std::vector<double> v(filt.size());
                        for (const LASPoint *pt : filt)
                            v[i++] = pt->z;
                        std::sort(v.begin(), v.end());
                        unsigned int size = v.size();
                        if (size % 2 == 0) {
                            result[0] = (v[(int) size / 2] + v[(int) size / 2 - 1]) / 2.0;
                        } else {
                            result[0] = v[(int) size / 2];
                        }
                    }
                }
            };

            class CellMin : public CellStats {
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        double min = G_DBL_MAX_POS;
                        for (const LASPoint *v : filt) {
                            if (v->z < min)
                                min = v->z;
                        }
                        result[0] = min;
                    }
                }
            };

            class CellMax : public CellStats {
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        double max = G_DBL_MAX_NEG;
                        for (const LASPoint *v : filt) {
                            if (v->z > max)
                                max = v->z;
                        }
                        result[0] = max;
                    }
                }
            };

            class CellSampleVariance : public CellStats {
            private:
                CellMean m_mean;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        double mean[1];
                        m_mean.compute(filt, mean);
                        double sum = 0;
                        for (const LASPoint *v : filt)
                            sum += g_sq(g_abs(v->z - mean[0]));
                        result[0] = sum / (filt.size() - 1);
                    }
                }
            };

            class CellPopulationVariance : public CellStats {
            private:
                CellMean m_mean;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        double mean[1];
                        m_mean.compute(filt, mean);
                        double sum = 0;
                        for (const LASPoint *v : filt)
                            sum += g_sq(g_abs(v->z - mean[0]));
                        result[0] = sum / filt.size();
                    }
                }
            };

            class CellSampleStdDev : public CellStats {
            private:
                CellSampleVariance m_variance;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        double var[1];
                        m_variance.compute(filt, var);
                        result[0] =  std::sqrt(var[0]);
                    }
                }
            };

            class CellPopulationStdDev : public CellStats {
            private:
                CellPopulationVariance m_variance;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] =  -9999.0;
                    } else {
                        double var[1];
                        m_variance.compute(filt, var);
                        result[0] =  std::sqrt(var[0]);
                    }
                }
            };

            class CellSkewness : public CellStats {
            private:
                CellMean m_mean;
                CellSampleStdDev m_stdDev;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] =  -9999.0;
                    } else {
                        // Fisher-Pearson
                        double mean[1];
                        m_mean.compute(filt, mean);
                        double sum = 0.0;
                        unsigned int count = filt.size();
                        for (const LASPoint *v : filt)
                            sum += std::pow(v->z - mean[0], 3.0) / count;
                        double sd[1];
                        m_stdDev.compute(filt, sd);
                        result[0] =  sum / std::pow(sd[0], 3.0);
                    }
                }
            };

            class CellKurtosis : public CellStats {
            private:
                CellMean m_mean;
                CellSampleStdDev m_stdDev;
            public:

                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] =  -9999.0;
                    } else {
                        double mean[1];
                        m_mean.compute(filt, mean);
                        double sum = 0.0;
                        unsigned int count = values.size();
                        for (const LASPoint *v : filt)
                            sum += std::pow(v->z - mean[0], 4.0) / count;
                        double sd[1];
                        m_stdDev.compute(filt, sd);
                        result[0] =  sum / std::pow(sd[0], 4.0) - 3.0;
                    }
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

                void compute(const std::list<LASPoint*> &values, double *result) {
                    result[0] =  0;
                }
            };

            class CellRugosity : public CellStats { // TODO: Seems to be density-dependent
            private:
                CellDensity m_density;
                double m_avgDensity;

                double computeArea(const Face &face) {
                    Point_3 p1 = face.vertex(0)->point();
                    Point_3 p2 = face.vertex(1)->point();
                    Point_3 p3 = face.vertex(2)->point();
                    return computeArea(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p3.x(), p3.y(), p3.z());
                }

                double computeArea(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
                    double side0 = std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0) + std::pow(z1 - z2, 2.0));
                    double side1 = std::sqrt(std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0) + std::pow(z2 - z3, 2.0));
                    double side2 = std::sqrt(std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0) + std::pow(z3 - z1, 2.0));
                    double s = (side0 + side1 + side2) / 2.0;
                    return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
                }

                double toPlane(const Point_3 &p, const Plane_3 &plane, const Point_3 &centroid) {
                    return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
                }

                double densityFactor(const std::list<LASPoint*> &values) {
                    if (values.size() == 0 || m_avgDensity <= 0.0 || m_density.area() <= 0.0)
                        return 1.0;
                    double density[1];
                    m_density.compute(values, density);
                    return 1.0 / (2.49127261 + 9.01659384 * std::sqrt(density[0] * 32.65748276));
                }

                double polyArea(const std::list<Point_3> &hull, const Plane_3 &plane, const Point_3 &centroid) {
                    double area = 0.0;
                    auto it0 = hull.begin();
                    auto it1 = hull.begin();
                    it1++;
                    do {
                        double z0 = toPlane(*it0, plane, centroid);
                        double z1 = toPlane(*it1, plane, centroid);
                        area += computeArea(it0->x(), it0->y(), z0, it1->x(), it1->y(), z1, centroid.x(), centroid.y(), centroid.z());
                        it0++;
                        it1++;
                        if (it1 == hull.end())
                            it1 = hull.begin();
                    } while (it0 != hull.end());
                    return area;
                }

            public:

                CellRugosity(double cellArea = 0.0, double avgDensity = 0.0) :
                    m_avgDensity(avgDensity) {
                    if(cellArea <= 0.0)
                        g_argerr("The cell area must be greater than zero.");
                    m_density.setArea(cellArea);
                }

                /**
                 * Using Du Preez, 2014 - Arc-Chord Ratio (ACR) Index.
                 */
                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] =  -9999.0;
                    } else {
                        std::list<Point_3> pts;
                        for (const LASPoint *v : filt)
                            pts.push_back(Point_3(v->x, v->y, v->z));

                        // Delaunay 3D surface area.
                        double tarea = 0.0;
                        Delaunay dt(pts.begin(), pts.end());
                        for (Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
                            tarea += computeArea(*it);

                        // Convex hull and POBF
                        std::list<Point_3> hull;
                        Plane_3 plane;
                        Point_3 centroid;
                        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), Gt());
                        CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

                        // POBF surface area.
                        double parea = polyArea(hull, plane, centroid);
                        double df = m_cellArea > 0.0 ? densityFactor(filt) : 1.0;
                        result[0] =  (tarea / parea) * df;
                    }
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
                double m_threshold;
                
                void fcLidarBLa(const std::list<LASPoint*> &values, double *result) {
                    double gnd = 0.0;
                    double all = 0.0;
                    for (const LASPoint *pt : values) {
                        if (pt->ground())
                            gnd += pt->intensity;
                        if (pt->cls < 2) // TODO: This should perhaps be filtered by class to remove bogus points.
                            all += pt->intensity;
                    }
                    //g_debug(" -- fcLidarBLa " << gnd << ", " << all << ", " << values.size());
                    result[0] =  all != 0.0 ? 1.0 - std::sqrt(gnd / all) : -9999.0;
                }

                void fcLidarBLb(const std::list<LASPoint*> &values, double *result) {
                    double gndSingle = 0.0, gndLast = 0.0, first = 0.0, single = 0.0, intermediate = 0.0, last = 0.0, total = 0.0;
                    for (const LASPoint *pt : values) {
                        if (pt->ground()) {
                            if (pt->single())
                                gndSingle += pt->intensity;
                            if (pt->last())
                                gndLast += pt->intensity;
                        }
                        if (pt->first())
                            first += pt->intensity;
                        if (pt->single())
                            single += pt->intensity;
                        if (pt->intermediate())
                            intermediate += pt->intensity;
                        if (pt->last())
                            last += pt->intensity;
                        total += pt->intensity; // TODO: This should perhaps be filtered by class to remove bogus points.
                    }
                    if (total == 0.0) {
                        result[0] =  -9999.0;
                    } else {
                        double denom = (first + single) / total + std::sqrt((intermediate + last) / total);
                        if (denom == 0.0) {
                            result[0] =  -9999.;
                        } else {
                            result[0] =  (gndSingle / total + std::sqrt(gndLast / total)) / denom;
                        }
                    }
                }

                void fcLidarIR(const std::list<LASPoint*> &values, double *result) {
                    double canopy = 0.0, total = 0.0;
                    for (const LASPoint *pt : values) {
                        if (!pt->ground())
                            canopy += pt->intensity;
                        total += pt->intensity;
                    }
                    result[0] =  total != 0.0 ? canopy / total : -9999.0;
                }

                void fcLidarRR(const std::list<LASPoint*> &values, double *result) {
                    unsigned int canopy = 0, total = 0;
                    for (const LASPoint *pt : values) {
                        if (!pt->ground())
                            ++canopy;
                        ++total;
                    }
                    result[0] =  total != 0.0 ? (double) canopy / total : -9999.0;
                }

                void fcLidarFR(const std::list<LASPoint*> &values, double *result) {
                    unsigned int canopy = 0, total = 0;
                    for (const LASPoint *pt : values) {
                        if (pt->first()) {
                            if (!pt->ground())
                                ++canopy;
                            ++total;
                        }
                    }
                    result[0] =  total != 0.0 ? (double) canopy / total : -9999.0;
                }
                
                void ccf(const std::list<LASPoint*> &values, double *result) {
                    if(values.size() < 75) {
                        result[0] = -9999.0;
                    } else {
                        double maxZ = -9999.0;
                        for(const LASPoint *pt : values)
                            maxZ = g_max(maxZ, pt->z);
                        double htIncrement = (maxZ - m_threshold) / 20.0;
                        double curHeight = m_threshold;
                        for(int band = 0; band < 20; ++band) {
                            double count = 0;
                            for(const LASPoint *pt : values) {
                                if(pt->z > curHeight)
                                    ++count;
                            }
                            result[band] = (double) count / values.size();
                            curHeight += htIncrement;
                        }
                    }
                }
                
                void gap(const std::list<LASPoint*> &values, double *result) {
                    if(!values.size()) {
                        result[0] = -9999.0;
                    } else {
                        int cnt = 0;
                        for(const LASPoint *p : values) {
                            if(p->z > m_threshold)
                                ++cnt;
                        }
                        result[0] = 1.0 - ((double) cnt / values.size());
                    }
                }
                
            public:

                CellGapFraction(unsigned char type, double threshold) : CellStats(),
                    m_type(type),
                    m_threshold(threshold) {
                }

                void threshold(double t) {
                    m_threshold = t;
                }
                
                int bands() const {
                    switch(m_type) {
                        case GAP_CCF:
                            return 20;
                        default:
                            return 1;
                    }
                }
                
                void compute(const std::list<LASPoint*> &values, double *result) {
                    std::list<LASPoint*> filt = filtered(values);
                    if (!filt.size()) {
                        result[0] = -9999.0;
                    } else {
                        switch (m_type) {
                            case GAP_BLA: 
                                fcLidarBLa(filt, result);
                                break;
                            case GAP_BLB: 
                                fcLidarBLb(filt, result);
                                break;
                            case GAP_IR: 
                                fcLidarIR(filt, result);
                                break;
                            case GAP_RR: 
                                fcLidarRR(filt, result);
                                break;
                            case GAP_FR: 
                                fcLidarFR(filt, result);
                                break;
                            case GAP_CCF:
                                ccf(filt, result);
                                break;
                            case GAP_GAP:
                                gap(filt, result);
                                break;
                            default:
                                g_argerr("Unknown Gap Fraction method: " << m_type);
                        }
                    }
                }
            };

        }
    }
}

#endif
