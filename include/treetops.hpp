#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>

#include "geotools.hpp"
#include "util.hpp"

namespace geotools {

    namespace treetops {

        namespace config {

            // Contains configuration information for performing tree top extraction.
            class G_DLL_EXPORT TreetopsConfig {
            public:

                // Defines the boundaries of work to be performed. Every step of the process
                // will be confined, including smoothing and searching.
                geotools::util::Bounds bounds;

                // A spatial reference ID for the output files.
                int srid;

                // If true, build the index on the tops table. Can be slow.
                bool buildIndex;

                // The cache size (B) for the sqlite database. A performance optimization.
                int tableCacheSize;

                // The cache size (B) for rows when reading the raster.
                int rowCacheSize;

                // The number of threads to use in execution.
                uint8_t threads;

                // Set to true to perform smoothing on the input raster.
                // This will force a check that the smoothing params are valild.
                // A smoothed filename must be provided for output, plus a sigma and window size.
                bool doSmoothing;

                // The size of the smoothing window >=3; an odd number.
                int smoothWindowSize;

                // The std. deviation used for generating the Gaussian kernel.
                // 0 < n <= 1.
                double smoothSigma;

                std::string smoothOriginalCHM;
                std::string smoothSmoothedCHM;

                // If true, treetop location will be performed.
                bool doTops;

                // The minimum height of a pixel that will be considered for inclusion as a 
                // tree top.
                double topsMinHeight;

                // The size of the top finding window >=3; an odd number.
                int topsWindowSize;

                std::string topsOriginalCHM;
                std::string topsSmoothedCHM;
                std::string topsTreetopsDatabase;

                // Set to true to delineate crowns.
                bool doCrowns;

                // The maximum crown radius.
                double crownsRadius;

                // The maximum height of a crown as a fraction of top height.
                double crownsHeightFraction;

                // The maximum height of pixels to consider for inclusion.
                double crownsMinHeight;

                // The input raster -- ideally the same one used for tops.
                std::string crownsSmoothedCHM;

                // The treetops database file.
                std::string crownsTreetopsDatabase;

                std::string crownsCrownsRaster;
                std::string crownsCrownsDatabase;

                // Build a TreetopsConfig with defaults.
                TreetopsConfig();

                void checkSmoothing() const;

                void checkTops() const;

                void checkCrowns() const;

                void checkMerge() const;

                // Check the validity of the configuration.
                void check() const;

                // Returns true if any function can be successfully run.
                bool canRun() const;
            };

        } // config

        namespace util {

            // A simple class for maintaining information about a tree top.
            class Top {
            public:
                uint64_t id;
                double x, y, z;
                double sx, sy, sz; // The smoothed z value.
                int col, row;
				int scol, srow; // Smoothed peak

                Top(uint64_t id, double x, double y, double z, double sx, double sy, double sz, int col, int row, int scol, int srow);

                Top();
            };

        } // util

        class G_DLL_EXPORT Treetops {
        private:
            geotools::util::Callbacks *m_callbacks;

        public:

            void setCallbacks(geotools::util::Callbacks *callbacks);

            // A convenience method for smoothing the input raster before using it to generate crowns
            // or treetops.
            void smooth(const geotools::treetops::config::TreetopsConfig &config, bool *cancel = nullptr);

            // Locates tree top points on a canopy height model.
            void treetops(const geotools::treetops::config::TreetopsConfig &config, bool *cancel = nullptr);

            // Performs tree crown delineation using a (preferrably smoothed) input raster and a
            // vector file (sqlite) containing tree tops as seeds. Output is an integer raster with 
            // cell values representing tree top IDs, and an optional vector which is the polygonized 
            // version of the raster. The table should have been generated using the treetops() method
            // to ensure that its structure is correct.
            void treecrowns(const geotools::treetops::config::TreetopsConfig &config, bool *cancel = nullptr);

            // Using predefined criteria, merge disparate tree crowns into single crowns,
            // potentially with more than one top.
            void merge(const geotools::treetops::config::TreetopsConfig &config, bool *cancel = nullptr);

        };

    } // trees

} // geotools

#endif
