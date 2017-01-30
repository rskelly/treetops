#ifndef __SQLITE_HPP__
#define __SQLITE_HPP__

#include <sqlite3.h>
#include <spatialite.h>

#include <sstream>
#include <iomanip>
#include <vector>

#include <ogr_core.h>
#include <ogr_feature.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "geotools.hpp"
#include "util.hpp"

using namespace geotools::util;

namespace geotools {

    namespace db {

        enum GeomType {
            GTUnknown = 0,
            GTPoint = 1,
            GTLine = 2,
            GTPolygon = 3,
            GTMultiPoint = 4,
            GTMultiLine = 5,
            GTMultiPolygon = 6
        };

        enum FieldType {
            FTUnknown = 0,
            FTInt = 1,
            FTDouble = 2,
            FTString = 3,
            FTBlob = 4
        };


        // Provides some easy database read/write methods.
        class DB {
        protected:
            GeomType m_type;
            int m_srid;
            std::string m_file;
            std::string m_layerName;
            std::string m_driver;
            std::string m_geomName;
            std::unordered_map<std::string, FieldType> m_fieldTypes;
            GDALDataset *m_ds;
            OGRLayer *m_layer;
            OGRFeatureDefn *m_fdef;

        public:

            DB(const std::string &file, const std::string &layer, const std::string &driver,
            		const std::unordered_map<std::string, FieldType> &fields,
            		GeomType type, int srid = 0, bool replace = false);

            DB(const std::string &file, const std::string &layer);

            ~DB();

            // Returns a map with file extensions for keys, and a list of driver names
            // as values.
            static std::map<std::string, std::set<std::string> > extensions();

            // Returns a map with driver short names as keys and long names for values.
            static std::map<std::string, std::string> drivers();

            // Returns a vector driver for the given filename.
            static std::string getDriverForFilename(const std::string &filename);

            void clear();

            void addPoint(double x, double y, double z, const std::unordered_map<std::string, std::string> &fields);

            void addPoints(std::vector<geotools::util::Point*> &points);

            void setCacheSize(size_t size);

            void dropGeomIndex();

            void createGeomIndex();

            void getPoints(std::vector<geotools::util::Point*> &points,
                    const geotools::util::Bounds &bounds);

            void getPoints(std::vector<geotools::util::Point*> &points,
                    int count, int offset = 0);

            void getNearestPoints(const geotools::util::Point &target, int count,
            		std::vector<std::unique_ptr<geotools::util::Point> > &points);

            uint64_t getGeomCount() const;

            void updateFeature(const std::string &idField, uint64_t id, std::unordered_map<std::string, void*> &values);

            void deleteFeature(const std::string &idField, uint64_t id);

            void execute(std::string &sql);

            void begin();

            void rollback();

            void commit();

            int srid() const;

        };

    } // db

} // geotools

#endif

