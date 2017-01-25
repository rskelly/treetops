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

                std::string smoothOriginalCHM;
                std::string smoothSmoothedCHM;

                // If true, treetop location will be performed.
                bool doTops;

                // For pixels equal or above each height (float) use the given
                // window size to detect maxima. Previously-detected maxima will be
                // obliterated if a new maximum is found whose window encompases
                // the previous one.
                std::vector<std::pair<float, uint8_t> > topsThresholds;

                std::string topsOriginalCHM;
                std::string topsSmoothedCHM;
                std::string topsTreetopsDatabase;

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

                // The input raster -- ideally the same one used for tops.
                std::string crownsSmoothedCHM;

                // The original CHM from which the smoothed raster was created.
                std::string crownsOriginalCHM;

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

                std::string thresholds() const;

                void parseThresholds(const std::string &str);

            };

        } // config

        namespace util {

            // A simple class for maintaining information about a tree top.
            class Top {
            public:
                uint64_t id;
                uint64_t parentID;
                double ox, oy, oz; // Original x, y, z value
                double sx, sy, sz; // Smoothed x, y, z value.
                int32_t sc, sr;    // Smoothed col, row.
                Top(uint64_t id, uint64_t parentId, 
                    double ox, double oy, double oz, 
                    double sx, double sy, double sz, 
                    int32_t sc, int32_t sr);

                Top();
            };

            using namespace geotools::db;

            class TTDB : public DB {
            public:

                TTDB(const std::string &file, const std::string &layer,
                    const std::unordered_map<std::string, FieldType> &fields,
                    GeomType type, int srid = 0, bool replace = false) : 
                    DB(file, layer, fields, type, srid, replace) {}

                TTDB(const std::string &file, const std::string &layer) :
                    DB(file, layer) {}

                void addTop(Top *top) {
                    if(m_type != GeomType::GTPoint)
                        g_runerr("This dataset is not a point dataset.");
                    OGRFeature feat(m_fdef);
                    feat.SetField("id", (GIntBig) top->id);
                    feat.SetField("parentId", (GIntBig) top->parentID);
                    feat.SetField("originalX", top->ox);
                    feat.SetField("originalY", top->oy);
                    feat.SetField("originalZ", top->oz);
                    feat.SetField("smoothedX", top->sx);
                    feat.SetField("smoothedY", top->sy);
                    feat.SetField("smoothedZ", top->sz);
                    feat.SetField("smoothedCol", top->sc);
                    feat.SetField("smoothedRow", top->sr);
                    OGRPoint geom(top->sx, top->sy, top->sz);
                    feat.SetGeometry(&geom);
                    if(CPLE_None != m_layer->CreateFeature(&feat))
                        g_runerr("Failed to add feature to " << m_file << ".");
                }

                void addTops(std::list<Top*> &tops) {
                    for(Top *t : tops)
                        addTop(t);
                }

                void getTops(std::vector<Top*> &tops, const geotools::util::Bounds &bounds) {
                    if(m_type != GeomType::GTPoint)
                        g_runerr("This dataset is not a point dataset.");
                    m_layer->SetSpatialFilterRect(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy());
                    OGRFeature *feat;
                    while((feat = m_layer->GetNextFeature())) {
                        Top *t = new Top(
                            feat->GetFieldAsInteger64("id"),
                            feat->GetFieldAsInteger64("parentID"),
                            feat->GetFieldAsDouble("originalX"),
                            feat->GetFieldAsDouble("originalY"),
                            feat->GetFieldAsDouble("originalZ"),
                            feat->GetFieldAsDouble("smoothedX"),
                            feat->GetFieldAsDouble("smoothedY"),
                            feat->GetFieldAsDouble("smoothedZ"),
                            feat->GetFieldAsInteger("smoothedCol"),
                            feat->GetFieldAsInteger("smoothedRow")
                        );
                        tops.push_back(t);
                        OGRFeature::DestroyFeature(feat);
                    }
                    m_layer->SetSpatialFilter(NULL);
                }

                void updateTop(Top *top) {
                    std::stringstream ss;
                    ss << "\"id\"=" << top->id;
                    m_layer->SetAttributeFilter(ss.str().c_str());
                    OGRFeature *feat = m_layer->GetNextFeature();
                    if(!feat)
                        g_runerr("Failed to find feature with ID: id=" << top->id);
                    feat->SetField("parentID", (GIntBig) top->parentID);
                    feat->SetField("originalX", top->ox);
                    feat->SetField("originalY", top->oy);
                    feat->SetField("originalZ", top->oy);
                    feat->SetField("smoothedX", top->sx);
                    feat->SetField("smoothedY", top->sy);
                    feat->SetField("smoothedZ", top->sz);
                    feat->SetField("smoothedCol", top->sc);
                    feat->SetField("smoothedRow", top->sr);
                    if(CPLE_None != m_layer->SetFeature(feat))
                        g_runerr("Failed to save feature: " << top->id << ".");
                    OGRFeature::DestroyFeature(feat);
                }

                void updateTops(std::vector<Top*> &tops) {
                    for(Top *top : tops)
                        updateTop(top);
                }

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
