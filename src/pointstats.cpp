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
#include <cmath>
#include <thread>
#include <chrono>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "pointstats.hpp"
#include "lasutil.hpp"
#include "lasreader.hpp"
#include "laspoint.hpp"
#include "cellstats.hpp"

using namespace geotools::util;
using namespace geotools::raster;
using namespace geotools::point;
using namespace geotools::las;
using namespace geotools::point::stats;

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
        int idxa = ((int) (ab.miny() / m_rowSize)) * ((int) (ab.width() / m_colSize)) + ((int) (ab.minx() / m_colSize));
        int idxb = ((int) (bb.miny() / m_rowSize)) * ((int) (bb.width() / m_colSize)) + ((int) (bb.minx() / m_colSize));
        return idxa < idxb;
    }
};

namespace geotools {

    namespace point {

        namespace pointstats_config {

            std::map<std::string, uint8_t> types = {
                {"Minimum", TYPE_MIN},
                {"Maximum", TYPE_MAX},
                {"Mean", TYPE_MEAN},
                {"Density", TYPE_DENSITY},
                {"Sample Variance", TYPE_VARIANCE},
                {"Sample Std. Dev.", TYPE_STDDEV},
                {"Population Variance", TYPE_PVARIANCE},
                {"Population Std. Dev.", TYPE_PSTDDEV},
                {"Count", TYPE_COUNT},
                {"Quantile", TYPE_QUANTILE},
                {"Median", TYPE_MEDIAN},
                {"Rugosity", TYPE_RUGOSITY},
                {"Kurtosis", TYPE_KURTOSIS},
                {"Skewness", TYPE_SKEW},
                {"Gap Fraction", TYPE_GAP_FRACTION}, 
                {"CoV", TYPE_COV}
            };
            std::map<std::string, uint8_t> attributes = {
                {"Height", ATT_HEIGHT},
                {"Intensity", ATT_INTENSITY}
            };
            std::map<std::string, uint8_t> gapFractionTypes = {
                {"IR", GAP_IR},
                {"BLa", GAP_BLA},
                {"BLb", GAP_BLB},
                {"RR", GAP_RR},
                {"FR", GAP_FR},
                {"CCF", GAP_CCF},
                {"GAP", GAP_GAP}
            };
            
            std::map<std::string, uint8_t> snapModes = {
                {"None" , SNAP_NONE},
                {"Grid" , SNAP_GRID},
                {"Origin" , SNAP_ORIGIN}
            };
            
            std::map<std::string, uint8_t> areaModes = {
                {"Full Cell", AREA_CELL},
                {"Radius", AREA_RADIUS}
            };

        } // config

        PointStatsConfig::PointStatsConfig() : 
            snapMode(SNAP_NONE),
            normalize(false),
            resolutionX(10.0),
            resolutionY(-10.0),
            originX(0.0),
            originY(0.0),
            gapFractionType(GAP_GAP),
            gapThreshold(0.0),
            threads(1),
            hsrid(0),
            vsrid(0),
            attribute(ATT_HEIGHT),
            angleLimit(90),
            quantile(0),
            quantiles(1),
            quantileFilter(0),
            quantileFilterFrom(0),
            quantileFilterTo(0),
            areaMode(AREA_CELL),
            areaSize(0) {
        }
        
        uint8_t PointStatsConfig::parseAtt(const std::string &attStr) {
            if ("intensity" == attStr) {
                return ATT_INTENSITY;
            } else if ("height" == attStr) {
                return ATT_HEIGHT;
            }
            return 0;
        }

        std::vector<uint8_t> PointStatsConfig::parseTypes(const std::vector<std::string> &typeStrs) {
            std::vector<uint8_t> types;
            for (const std::string &typeStr : typeStrs) {
                if ("min" == typeStr) {
                    types.push_back(TYPE_MIN);
                } else if ("max" == typeStr) {
                    types.push_back(TYPE_MAX);
                } else if ("mean" == typeStr) {
                    types.push_back(TYPE_MEAN);
                } else if ("density" == typeStr) {
                    types.push_back(TYPE_DENSITY);
                } else if ("variance" == typeStr) {
                    types.push_back(TYPE_VARIANCE);
                } else if ("stddev" == typeStr) {
                    types.push_back(TYPE_STDDEV);
                } else if ("pvariance" == typeStr) {
                    types.push_back(TYPE_PVARIANCE);
                } else if ("pstddev" == typeStr) {
                    types.push_back(TYPE_PSTDDEV);
                } else if ("count" == typeStr) {
                    types.push_back(TYPE_COUNT);
                } else if ("median" == typeStr) {
                    types.push_back(TYPE_MEDIAN);
                } else if ("skew" == typeStr) {
                    types.push_back(TYPE_SKEW);
                } else if ("rugosity" == typeStr) {
                    types.push_back(TYPE_RUGOSITY);
                } else if ("kurtosis" == typeStr) {
                    types.push_back(TYPE_KURTOSIS);
                } else if ("gap" == typeStr) {
                    types.push_back(TYPE_GAP_FRACTION);
                }
            }
            return types;
        }

        uint8_t PointStatsConfig::parseGap(const std::string &gapStr) {
            if ("bla" == gapStr) {
                return GAP_BLA;
            } else if ("blb" == gapStr) {
                return GAP_BLB;
            } else if ("fr" == gapStr) {
                return GAP_FR;
            } else if ("rr" == gapStr) {
                return GAP_RR;
            } else if ("ir" == gapStr) {
                return GAP_IR;
            } else if("ccf" == gapStr) {
                return GAP_CCF;
            } else if("gap" == gapStr) {
                return GAP_GAP;
            }
            return 0;
        }

        void PointStats::checkConfig(const PointStatsConfig &config) {
            if (g_abs(config.resolutionX) == 0.0 || config.resolutionY == 0.0)
                g_argerr("Resolution must be != 0: " << config.resolutionX << ", " << config.resolutionY);
            if (!config.sourceFiles.size())
                g_argerr("At least one input file is required.");
            if (!config.dstFiles.size())
                g_argerr("At least one output file is required.");
            if (config.attribute == 0)
                g_argerr("An attribute is required.");
            if (!config.types.size())
                g_argerr("At least one valid type is required.");
            if (config.classes.size() == 0)
                g_warn("No classes given. Matching all classes.");
            if (config.angleLimit <= 0)
                g_argerr("Angle limit must be greater than zero.");
            if (config.dstFiles.size() != config.types.size())
                g_argerr("There should be one output file for each type.");
            if(config.areaMode != AREA_CELL && config.areaSize <= 0) 
                g_argerr("If area mode is not cell, a size must be given.");
          
            g_debug("Resolution: " << config.resolutionX << ", " << config.resolutionY);
            g_debug("Files: " << config.sourceFiles.size());
            g_debug("Destinations: " << config.dstFiles.size());
            g_debug("Attribute: " << config.attribute);
            g_debug("Types: " << config.types.size());
            g_debug("Classes: " << config.classes.size());
            g_debug("Angle Limit: " << config.angleLimit);
        }

        CellStats* PointStats::getComputer(const uint8_t &type, const PointStatsConfig &config) {
            using namespace geotools::point::stats;
            switch (type) {
                case TYPE_MEAN:      	return new CellMean();
                case TYPE_MEDIAN:    	return new CellMedian();
                case TYPE_COUNT:     	return new CellCount();
                case TYPE_STDDEV:    	return new CellSampleStdDev();
                case TYPE_VARIANCE:   	return new CellSampleVariance();
                case TYPE_PSTDDEV: 		return new CellPopulationStdDev();
                case TYPE_PVARIANCE: 	return new CellPopulationVariance();
                case TYPE_DENSITY: 		return new CellDensity(g_sq(config.resolutionX));
                case TYPE_RUGOSITY: 	return new CellRugosity(g_sq(config.resolutionX), 0.0);
                case TYPE_MAX: 			return new CellMax();
                case TYPE_MIN: 			return new CellMin();
                case TYPE_KURTOSIS: 	return new CellKurtosis();
                case TYPE_SKEW: 		return new CellSkewness();
                case TYPE_QUANTILE: 	return new CellQuantile(config.quantile, config.quantiles);
                case TYPE_GAP_FRACTION:	return new CellGapFraction(config.gapFractionType, config.gapThreshold);
                case TYPE_COV:			return new CellCoV();
                default:
                    g_argerr("Invalid statistic type: " << type);
            }
        }

        void _runner(PointStats *ps) {
            ps->runner();
        }

        PointStats::PointStats() :
            m_callbacks(nullptr),
			m_running(false), m_cancel(nullptr),
            m_finalizedCount(0), m_cellCount(0),
			m_cols(0), m_rows(0),
			m_resolutionX(0), m_resolutionY(0),
			m_tlx(0), m_tly(0) {
        }

        PointStats::~PointStats() {
        }

        void PointStats::runner() {
            uint64_t idx;
            std::list<LASPoint*> pts;
            while (true) {
                while(m_bq.empty() && m_running && !*m_cancel)
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                if(!m_running || *m_cancel)
                    break;
                {
                    // Grab an index.
                    std::unique_lock<std::mutex> lk(m_qmtx);
                    if(m_bq.empty())
                        continue;
                    idx = m_bq.front();
                    m_bq.pop();
                }
                
                {
                    // Get the relevant points from the cache
                    std::unique_lock<std::mutex> lk(m_cmtx);
                    pts.assign(m_cache[idx].begin(), m_cache[idx].end());
                    m_cache.erase(idx);
                    //g_debug(" -- cache size: " << m_cache.size());
                }

                // Process the points.
                if (!pts.empty()) {
                    // Get the cell's centre.
                	double x = m_bounds.toX(idx % m_cols, m_resolutionX) + m_resolutionX / 2.0;
                	double y = m_bounds.toY(idx / m_cols, m_resolutionY) + m_resolutionY / 2.0;
                    for (size_t i = 0; i < m_computers.size(); ++i) {
                        int bands = m_computers[i]->bands();
                        Buffer buf(sizeof(double) * bands);
                        m_computers[i]->compute(x, y, pts, (double *) buf.buf);
                        std::unique_lock<std::mutex> flk(*(m_mtx[i].get()));
                        int b = 0;
                        for(const std::unique_ptr<MemRaster<float> > &m : m_mem[i])
                            m->set(idx, ((double *) buf.buf)[b++]);
                    }
                    for (LASPoint *pt : pts)
                        delete pt;
                    pts.clear();

                    if(m_callbacks)
                        m_callbacks->overallCallback(((float) ++m_finalizedCount / m_cellCount) * 0.52f + 0.47f);

                }
            }
            g_debug(" -- exit thread");
        }
            
        bool _cancel = false;
        
        class InitCallback : public LASReaderCallback {
        public:
            const Callbacks *m_callbacks;
            InitCallback(const Callbacks *callbacks) :
                m_callbacks(callbacks) {
            }
            void status(float status) const {
                m_callbacks->overallCallback(0.01f + status * 0.45f);
            }
        };

        void PointStats::pointstats(const PointStatsConfig &config, const Callbacks *callbacks, bool *cancel) {

            checkConfig(config);

            // In no cancel flag is given, use the dummy.
            m_cancel = cancel == nullptr ? &_cancel : cancel;
            m_callbacks = callbacks;
            m_resolutionX = config.resolutionX;
            m_resolutionY = config.resolutionY;

            if(callbacks) {
                callbacks->overallCallback(0.01f);
                callbacks->statusCallback("Starting...");
            }

            FileSorter sorter(20.0, 20.0);
            std::vector<std::string> files(config.sourceFiles.begin(), config.sourceFiles.end());
            std::sort(files.begin(), files.end(), sorter);

            if(*cancel) return;            
            
            // Initialize the point stream
            LASMultiReader ps(config.sourceFiles, m_resolutionX, m_resolutionY, m_cancel);
            m_bounds = ps.bounds();
            switch(config.snapMode) {
                case SNAP_ORIGIN:
                    m_bounds.align(config.originX, config.originY, m_resolutionX, m_resolutionY);
                    break;
                case SNAP_GRID:
                    m_bounds.snap(m_resolutionX);
                    break;
            }
            if(callbacks)
                callbacks->statusCallback("Initializing point reader...");

            ps.setBounds(m_bounds);

            // Initialize the filter.
            LASFilter filter;
            filter.setClasses(config.classes);
            if(config.areaMode == AREA_RADIUS) {
                filter.setRadius(config.areaSize);
                filter.setGrid(config.resolutionX, config.resolutionY, config.originX, config.originY);
            }

            ps.setFilter(&filter);
            InitCallback initStatus(callbacks);
            ps.buildFinalizer(&initStatus);
            m_cellCount = ps.cellCount();

            m_cols = m_bounds.maxCol(m_resolutionX) + 1;
            m_rows = m_bounds.maxRow(m_resolutionY) + 1;

            if(callbacks) {
                callbacks->overallCallback(0.5f);
                callbacks->statusCallback("Configuring computers...");
            }

            // Create a computer, grid and mutex for each statistic.
            for (size_t i = 0; i < config.types.size(); ++i) {
                if(*cancel) return;
                // Create computers for each stat
                g_debug(" -- configuring computer " << (int) config.types[i]);
                std::unique_ptr<CellStats> cs(getComputer(config.types[i], config));
                int bands = cs->bands();
                m_computers.push_back(std::move(cs));

                // Create raster grid for each stat.
                g_debug(" -- configuring grids");
                std::vector<std::unique_ptr<MemRaster<float> > > rasters;
                for(int b = 0; b < bands; ++b) {
                    std::unique_ptr<MemRaster<float> > mr(new MemRaster<float>(m_cols, m_rows, true));
                    mr->fill(-9999.0);
                    mr->setNodata(-9999.0);
                    rasters.push_back(std::move(mr));
                }
                m_mem.push_back(std::move(rasters));

                // Create mutex for grid.
                g_debug(" -- creating mutexes");
                std::unique_ptr<std::mutex> m(new std::mutex());
                m_mtx.push_back(std::move(m));
            }

            if(callbacks) {
                callbacks->overallCallback(0.47f);
                callbacks->statusCallback("Processing points...");
            }

            // Initialize the thread group for runner.
            std::list<std::thread> threads;
            m_running = true;
            uint64_t finalIdx;
            bool final;
            
            // Start the runner threads.
            for (int i = 0; i < g_max(1, config.threads - 1); ++i) {
                std::thread t(_runner, this);
                threads.push_back(std::move(t));
            }

            // Begin streaming the points into the cache for processing.
            g_debug(" -- streaming points");
            LASPoint pt;
            m_finalizedCount = 0;
            while (ps.next(pt, &final, &finalIdx) && !*m_cancel) {
                {
                    uint64_t idx = m_bounds.toRow(pt.y, m_resolutionY) * m_cols
                    		+ m_bounds.toCol(pt.x, m_resolutionX);
                    std::unique_lock<std::mutex> lk(m_cmtx);
                    m_cache[idx].push_back(new LASPoint(pt));
                }
                if (final) {
                    {
                        std::unique_lock<std::mutex> lk(m_qmtx);
                        m_bq.push(finalIdx);
                    }
                }
                if(m_bq.size() > 5000) {
                    while(m_bq.size() > 100  && !*m_cancel)
                        std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
            }

            while(!m_bq.empty() &&  !*m_cancel)
                std::this_thread::sleep_for(std::chrono::milliseconds(1));

            g_debug(" -- done " << m_finalizedCount << " of " << m_cellCount);

            m_running = false;
            // Shut down and join the runners.
            for (std::thread &t : threads)
                t.join();

            // Write the grids to files.
            for (size_t i = 0; i < m_mem.size(); ++i) {
                if(*cancel) return;
                // Normalize the grids if desired.
                if (config.normalize) {
                    g_debug(" -- normalizing");
                    for(const std::unique_ptr<MemRaster<float> > &m : m_mem[i])
                        m->normalize();
                }
                // Prepare the grid
                // TODO: Only works with UTM north.
                int band = 1;
                Raster<float> grid(config.dstFiles[0].c_str(), m_mem[i].size(), m_bounds, m_resolutionX,
                    m_resolutionY, config.hsrid);
                for(const std::unique_ptr<MemRaster<float> > &m : m_mem[i]) {
                    if(*cancel) return;
                    // Write grid to file.
                    grid.setNodata(-9999.0, band);
                    grid.writeBlock(band, *(m.get()));
                    ++band;
                }
            }

            if (callbacks) {
                callbacks->overallCallback(1.0f);
                callbacks->statusCallback("Finished.");
            }

        }

    } // point

} // geotools

