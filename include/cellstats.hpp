#ifndef __CELLSTATS_HPP__
#define __CELLSTATS_HPP__

#include "laspoint.hpp"

#define GAP_IR 1
#define GAP_BLA 2
#define GAP_BLB 3
#define GAP_RR 4
#define GAP_FR 5
#define GAP_CCF 6
#define GAP_GAP 7

using namespace geotools::las;

namespace geotools {

    namespace point {

        namespace stats {

            class CellStatsFilter {
            protected:
                CellStatsFilter *m_chain;
                const std::list<LASPoint*> *m_points;

                virtual bool keepImpl(const LASPoint *) const = 0;

                virtual void init();

            public:
                CellStatsFilter();

                CellStatsFilter* chain(CellStatsFilter *next);

                void setPoints(const std::list<LASPoint*> *points);

                bool keep(const LASPoint *pt) const;

                virtual ~CellStatsFilter() = 0;
            };

            class ClassFilter : public CellStatsFilter {
            private:
                std::set<unsigned char> m_classes;

            protected:
                bool keepImpl(const LASPoint *pt) const;

                void init();
                
            public:
                ClassFilter(const std::set<unsigned char> &classes);

                ~ClassFilter();
            };

            class QuantileFilter : public CellStatsFilter {
            private:
                int m_quantiles;
                int m_from;
                int m_to;
                double m_min;
                double m_max;

            protected:
                bool keepImpl(const LASPoint *pt) const;

                void init();
                
            public:
                QuantileFilter(int quantiles, int from, int to);

                ~QuantileFilter();
            };

            class CellStats {
            private:
                CellStatsFilter *m_filter;

            public:
                CellStats();

                void setFilter(CellStatsFilter *filter);
                
                std::list<LASPoint*> filtered(const std::list<LASPoint*> &values);

                virtual void compute(const std::list<LASPoint*>&, double*);

                virtual int bands() const;
                
                ~CellStats();
            };

            class CellDensity : public CellStats {
            private:
                double m_cellArea;
            public:
                CellDensity(double area = 0.0);

                void compute(const std::list<LASPoint*>&, double*);
                
                void setArea(double area);

                double area();
            };

            class CellMean : public CellStats {
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellCount : public CellStats {
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellMedian : public CellStats {
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellMin : public CellStats {
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellMax : public CellStats {
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellSampleVariance : public CellStats {
            private:
                CellMean m_mean;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellPopulationVariance : public CellStats {
            private:
                CellMean m_mean;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellSampleStdDev : public CellStats {
            private:
                CellSampleVariance m_variance;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellPopulationStdDev : public CellStats {
            private:
                CellPopulationVariance m_variance;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellSkewness : public CellStats {
            private:
                CellMean m_mean;
                CellSampleStdDev m_stdDev;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellKurtosis : public CellStats {
            private:
                CellMean m_mean;
                CellSampleStdDev m_stdDev;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellCoV : public CellStats {
            private:
                CellMean m_mean;
                CellSampleStdDev m_stdDev;
            public:
                void compute(const std::list<LASPoint*>&, double*);
            };

            class CellQuantile : public CellStats {
            private:
                unsigned int m_quantile;
                unsigned int m_quantiles;
            public:
                CellQuantile(unsigned char quantile, unsigned char quantiles);
                
                void compute(const std::list<LASPoint*>&, double*);
            };

            // Using Du Preez, 2014 - Arc-Chord Ratio (ACR) Index.
            class CellRugosity : public CellStats { // TODO: Seems to be density-dependent
            private:
                CellDensity m_density;
                double m_avgDensity;

            public:
                CellRugosity(double cellArea = 0.0, double avgDensity = 0.0);

                void compute(const std::list<LASPoint*>&, double*);
            };

             // Adapted from:
             // Hopkinson, C., & Chasmer, L. (2009). Testing LiDAR models of fractional cover across 
             //	multiple forest ecozones. Remote Sensing of Environment, 113(1), 275–288. 
             //	http://doi.org/10.1016/j.rse.2008.09.012
            class CellGapFraction : public CellStats {
            private:
                unsigned char m_type;
                double m_threshold;

            public:
                CellGapFraction(unsigned char type, double threshold);

                void compute(const std::list<LASPoint*>&, double*);

                void threshold(double t);
                
                int bands() const;                
            };

        }
    }
}

#endif
