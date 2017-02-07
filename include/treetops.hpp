#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>
#include <vector>

#include "geotools.hpp"
#include "util.hpp"
#include "db.hpp"

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

                // The method to use for smoothing; e.g. Gaussian, Median
                std::string smoothMethod;

                // The path to the original CHM.
                std::string smoothOriginalCHM;

                // The path to the smoothed CHM.
                std::string smoothSmoothedCHM;

                // The raster band to smooth.
                int smoothOriginalCHMBand;

                // The driver to use for creating the smoothed CHM.
                std::string smoothSmoothedCHMDriver;

                // If true, treetop location will be performed.
                bool doTops;

                // For pixels equal or above each height (double) use the given
                // window size to detect maxima. Previously-detected maxima will be
                // obliterated if a new maximum is found whose window encompases
                // the previous one.
                std::vector<std::pair<double, uint8_t> > topsThresholds;

                // The path to the smoothed CHM.
                std::string topsSmoothedCHM;

                // The path to the treetops database.
                std::string topsTreetopsDatabase;

                // The river to use for creating the database.
                std::string topsTreetopsDatabaseDriver;

                // The max proportion of pixels in a given kernel that are allowed
                // to be null. Any kernel with this many or greater is ignored.
                double topsMaxNulls;

                // Set to true to delineate crowns.
                bool doCrowns;

                // The maximum crown radius.
                double crownsRadius;

                // The maximum height of a crown as a fraction of top height.
                double crownsHeightFraction;

                // The maximum height of pixels to consider for inclusion.
                double crownsMinHeight;

                // If true, the heights of the treetops will be updated
                // using values from the original CHM within the bounds
                // of the crowns
                bool crownsUpdateHeights;

                // The input raster -- ideally the same one used for tops.
                std::string crownsSmoothedCHM;

                // The original CHM from which the smoothed raster was created.
                std::string crownsOriginalCHM;

                // The treetops database file.
                std::string crownsTreetopsDatabase;

                // The path to the crowns raster.
                std::string crownsCrownsRaster;

                // If true, a crowns database will be produced
                bool crownsDoDatabase;

                // The path to the crowns database.
                std::string crownsCrownsDatabase;

                // The driver to use for the raster.
                std::string crownsCrownsRasterDriver;

                // The driver to use for the database.
                std::string crownsCrownsDatabaseDriver;

                // Build a TreetopsConfig with defaults.
                TreetopsConfig();

                // Check that the settings are appropriate for a smoothing
                // job, throw an exception otherwise.
                void checkSmoothing() const;

                // Check that the settings are appropriate for a treetops
                // job, throw an exception otherwise.
                void checkTops() const;

                // Check that the settings are appropriate for a crowns
                // job, throw an exception otherwise.
                void checkCrowns() const;

                // Check that the settings are appropriate for a merge
                // job, throw an exception otherwise.
                void checkMerge() const;

                // Check the validity of the configuration.
                void check() const;

                // Returns true if any function can be successfully run.
                bool canRun() const;

                // Returns the thresholds as a comma-delimited list.
                std::string thresholds() const;

                // Parses a comma-delimited list of thresholds into an internal list.
                void parseThresholds(const std::string &str);

            };

        } // config

        namespace util {

            // A simple class for maintaining information about a tree top.
            class Top {
            public:
            	uint64_t fid;		// The geomid.
                uint64_t id;		// The ID of this top.
                uint64_t parentID;	// The ID of this top's parent.
                double ox, oy, oz; 	// Original x, y, z value
                double sx, sy, sz; 	// Smoothed x, y, z value.
                int32_t sc, sr;    	// Smoothed col, row.
                Top(uint64_t id, uint64_t parentId, 
                    double ox, double oy, double oz, 
                    double sx, double sy, double sz, 
                    int32_t sc, int32_t sr);

                Top();
            };

            using namespace geotools::db;

            // A subclass of DB with specific methods for managing treetops.
            class TTDB : public DB {
            private:

            	// Returns a map of fields/types to use in the constructor.
                static std::unordered_map<std::string, FieldType> fields();

            public:

                // Build a new database.
                TTDB(const std::string &file, const std::string &layer, const std::string &driver, int srid = 0, bool replace = false);

                // Open an existing database
                TTDB(const std::string &file, const std::string &layer);

                // Add a single treetop
                void addTop(const std::unique_ptr<Top> &top);

                // Add a list of treetops.
                void addTops(const std::list<std::unique_ptr<Top> > &tops);

                // Return the treetops within the given bounds
                void getTops(std::list<std::unique_ptr<Top> > &tops, const geotools::util::Bounds &bounds);

                // Update a treetop.
                void updateTop(const std::unique_ptr<Top> &top);

                // Update a list of treetops
                void updateTops(std::list<std::unique_ptr<Top> > &tops);

            };

        } // util

        class G_DLL_EXPORT Treetops {
        private:
            geotools::util::Callbacks *m_callbacks;

            void updateOriginalCHMHeights(const geotools::treetops::config::TreetopsConfig &config,
            		bool *cancel, float start, float end);

            void delineateCrowns(const geotools::treetops::config::TreetopsConfig &config,
            		bool *cancel, float start, float end);

            void polygonizeCrowns(const geotools::treetops::config::TreetopsConfig &config,
            		bool *cancel, float start, float end);

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
