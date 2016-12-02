#ifndef __LASGRID_HPP__
#define __LASGRID_HPP__

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <mutex>
#include <condition_variable>

#include "geotools.h"
#include "raster.hpp"
#include "util.hpp"
#include "laspoint.hpp"
#include "cellstats.hpp"

#define TYPE_MIN 1
#define TYPE_MAX 2
#define TYPE_MEAN 3
#define TYPE_DENSITY 4
#define TYPE_VARIANCE 7
#define TYPE_STDDEV 8
#define TYPE_PVARIANCE 12
#define TYPE_PSTDDEV 13
#define TYPE_COUNT 9
#define TYPE_QUANTILE 10
#define TYPE_MEDIAN 11
#define TYPE_SKEW 16
#define TYPE_KURTOSIS 15
#define TYPE_RUGOSITY 14
#define TYPE_GAP_FRACTION 17

#define SNAP_NONE 0
#define SNAP_GRID 1
#define SNAP_ORIGIN 2

#define LAS_EXT ".las"

#define ATT_HEIGHT 1
#define ATT_INTENSITY 2

namespace geotools {

    namespace point {

        namespace pointstats_config {

            extern std::map<std::string, uint8_t> types;
            extern std::map<std::string, uint8_t> attributes;
            extern std::map<std::string, uint8_t> gapFractionTypes;
            extern std::map<std::string, uint8_t> snapModes;
        }

        /**
         * Contains configuration values for running the lasgrid process.
         */
        class PointStatsConfig {
        public:
            std::vector<std::string> dstFiles;
            std::vector<std::string> sourceFiles;
            std::vector<uint8_t> types;
            std::set<uint8_t> classes;
            geotools::util::Bounds bounds;
            bool fill;
            uint8_t snapMode;
            bool rebuild;
            bool normalize;
            double resolution;
            double originX;
            double originY;
            double gapThreshold;
            uint32_t threads;
            uint16_t hsrid;
            uint16_t vsrid;
            uint8_t attribute;
            uint8_t angleLimit;
            uint8_t quantile;
            uint8_t quantiles;
            uint8_t gapFractionType;
            uint32_t quantileFilter;
            uint32_t quantileFilterFrom;
            uint32_t quantileFilterTo;

            PointStatsConfig();
            
            /**
             * Interpret the attribute and  return the constant int value.
             */
            uint8_t parseAtt(const std::string &attStr);

            /**
             * Interpret the output type and return the constant int value.
             */
            std::vector<uint8_t> parseTypes(const std::vector<std::string> &typeStrs);

            uint8_t parseGap(const std::string &typeStr);

            /**
             * Returns true if the classes set contains the class.
             */
            bool hasClass(uint8_t cls) const {
                return classes.find(cls) != classes.end();
            }

            bool hasClasses() const {
                return classes.size() > 0;
            }

            /**
             * Returns true if the quantile filter does not pass all points.
             */
            bool hasQuantileFilter() const {
                return quantileFilterFrom == 0 && quantileFilterTo == quantileFilter - 1;
            }
            
            bool check() {
                if(!dstFiles.size())
                    return false;
                if(!sourceFiles.size())
                    return false;
                if(!types.size())
                    return false;
                if(!classes.size())
                    g_warn("No classes specified; using all.");
                if(!hsrid)
                    g_warn("No horizontal spatial reference ID specified.");
                if(!vsrid)
                    g_warn("No vertical spatial reference ID specified.");
                if(!attribute)
                    return false;
                return true;
            }

        };

        class PointStats {
        private:
            std::mutex m_cmtx;
            std::mutex m_qmtx;
            std::condition_variable m_cdn;

            bool m_running;
            bool *m_cancel;
            std::unordered_map<size_t, std::list<geotools::las::LASPoint*> > m_cache;
            std::vector<std::unique_ptr<geotools::point::stats::CellStats> > m_computers;
            std::vector<std::vector<std::unique_ptr<geotools::raster::MemRaster<float> > > > m_mem;
            std::vector<std::unique_ptr<std::mutex> > m_mtx;
            std::queue<size_t> m_bq;
            std::queue<size_t> m_idxq;

            /**
             * Check the configuration for validity. 
             * Throw an exception if it's invalid or absent.
             */
            void checkConfig(const PointStatsConfig &config);

            geotools::point::stats::CellStats* getComputer(const uint8_t &type, const PointStatsConfig &config);

            /**
             * Compute the working bounds and the selection of files
             * to include in the working set.
             */
            void computeWorkBounds(const std::vector<std::string> &files, const Bounds &bounds,
                    std::set<std::string> &selectedFiles, Bounds &workBounds, unsigned long *pointCount);

        public:

            PointStats();

            ~PointStats();

            /**
             * Runs the compute loop. Inteded to be used by a thread.
             */
            void runner();

            /**
             * Runs the read loop. Inteded to be used by a thread.
             */
            void reader();

            /**
             * Execute the gridding process.
             */
            void pointstats(const PointStatsConfig &config, const Callbacks *callbacks = nullptr, bool *cancel = nullptr);
        };

    } // las

} // geotools

#endif
