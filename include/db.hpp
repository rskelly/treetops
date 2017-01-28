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


        namespace util {

            GeomType geomType(OGRwkbGeometryType type);

            OGRwkbGeometryType geomType(GeomType type);

            FieldType fieldType(OGRFieldType type);

            OGRFieldType fieldType(FieldType type);

            void addFields(OGRFeature *feat, geotools::util::Point *pt, OGRFeatureDefn *ftdef);

            void updateFields(OGRFeature *feat, std::unordered_map<std::string, void*> &values, OGRFeatureDefn *ftdef);

        }

        using namespace geotools::db::util;

        // Provides some easy database read/write methods.
        class DB {
        protected:
            GeomType m_type;
            int m_srid;
            std::string m_file;
            std::string m_layerName;
            std::string m_geomName;
            std::unordered_map<std::string, FieldType> m_fieldTypes;
            GDALDataset *m_ds;
            OGRLayer *m_layer;
            OGRFeatureDefn *m_fdef;

        public:

            DB(const std::string &file, const std::string &layer,
            		const std::unordered_map<std::string, FieldType> &fields,
            		GeomType type, int srid = 0, bool replace = false) :
                m_type(type),
                m_srid(srid),
                m_file(file),
				m_layerName(layer),
				m_fieldTypes(fields),
				m_ds(nullptr),
				m_layer(nullptr),
				m_fdef(nullptr) {

                if(replace && Util::exists(file))
                    Util::rm(file);

                GDALAllRegister();
                std::string drvName = "SQLite";//Util::getDriverForFilename(m_file);
                GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(drvName.c_str());
                if(!drv)
                    g_runerr("Driver not found for " << m_file << " (" << drvName << ")");
                m_ds = drv->Create(m_file.c_str(), 0, 0, 0, GDT_Unknown, NULL);
                if(!m_ds)
                    g_runerr("Failed to create data set for " << m_file);
            	if(m_layerName.empty())
            		m_layerName = "data";
            	char **options = nullptr;
            	OGRSpatialReference sr;
            	sr.importFromEPSG(m_srid);
				m_layer = m_ds->CreateLayer(m_layerName.c_str(), &sr, geomType(m_type), options);
				if(!m_layer)
					g_runerr("Failed to create layer, " << m_layerName << ".");
            	m_fdef = m_layer->GetLayerDefn();
            	for(const auto &it : m_fieldTypes) {
            		OGRFieldDefn def(it.first.c_str(), fieldType(it.second));
            		m_fdef->AddFieldDefn(&def);
            	}
                OGRGeomFieldDefn *gdef = m_layer->GetLayerDefn()->GetGeomFieldDefn(0);
                m_geomName = std::string(gdef->GetNameRef());
            }

            DB(const std::string &file, const std::string &layer) :
                m_type(GeomType::GTUnknown),
                m_srid(0),
                m_file(file),
				m_layerName(layer),
				m_ds(nullptr),
				m_layer(nullptr),
				m_fdef(nullptr) {

                GDALAllRegister();
            	m_ds = (GDALDataset *) GDALOpenEx(m_file.c_str(), GDAL_OF_VECTOR|GDAL_OF_UPDATE, NULL, NULL, NULL);
                if(!m_ds)
                    g_runerr("Failed to open data set for " << m_file);
            	if(m_layerName.empty()) {
            		m_layer = m_ds->GetLayer(1);
                	if(!m_layer)
                		g_runerr("No layer was found on this data set.");
            	} else {
            		m_layer = m_ds->GetLayerByName(m_layerName.c_str());
                	if(!m_layer)
                		g_runerr("No layer, " << m_layerName << " was found on this data set.");
            	}
            	m_type = geomType(m_layer->GetGeomType());
            	OGRGeomFieldDefn *gdef = m_layer->GetLayerDefn()->GetGeomFieldDefn(0);
            	m_geomName = std::string(gdef->GetNameRef());
            	m_fdef = m_layer->GetLayerDefn();
            	for(int i = 0; i < m_fdef->GetFieldCount(); ++i) {
            		OGRFieldDefn *def = m_fdef->GetFieldDefn(i);
            		m_fieldTypes[std::string(def->GetNameRef())] = fieldType(def->GetType());
            	}
            }

            ~DB() {
            	GDALClose(m_ds);
            }

            void clear();

            void addPoint(double x, double y, double z, const std::unordered_map<std::string, std::string> &fields) {
            	if(m_type != GeomType::GTPoint)
            		g_runerr("This dataset is not a point dataset.");
            	OGRFeature feat(m_fdef);
            	for(const auto &it : fields) {
            		switch(m_fieldTypes[it.first]) {
            		//case FieldType::Blob:
            		case FieldType::FTDouble:
                		feat.SetField(it.first.c_str(), atof(it.second.c_str()));
						break;
            		case FieldType::FTInt:
                		feat.SetField(it.first.c_str(), atoi(it.second.c_str()));
						break;
            		case FieldType::FTString:
                		feat.SetField(it.first.c_str(), it.second.c_str());
						break;
            		default:
            			g_runerr("Field type " << m_fieldTypes[it.first] << " is not currently supported.");
            		}
            	}
            	OGRPoint geom(x, y, z);
            	feat.SetGeometry(&geom);
                if(CPLE_None != m_layer->CreateFeature(&feat))
                    g_runerr("Failed to add feature to " << m_file << ".");
            }

            void addPoints(std::vector<geotools::util::Point*> &points) {
                for (const geotools::util::Point *pt : points)
                    addPoint(pt->x, pt->y, pt->z, pt->fields);
            }


            void setCacheSize(size_t size) {
            }

            void dropGeomIndex() {
            }

            void createGeomIndex() {
            }

            void getPoints(std::vector<geotools::util::Point*> &points,
                    const geotools::util::Bounds &bounds) {
            	if(m_type != GeomType::GTPoint)
            		g_runerr("This dataset is not a point dataset.");
            	m_layer->SetSpatialFilterRect(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy());
            	OGRFeature *feat;
            	while((feat = m_layer->GetNextFeature())) {
            		OGRPoint *opt = (OGRPoint *) feat->GetGeometryRef();
            		geotools::util::Point *pt = new Point(opt->getX(), opt->getY(), opt->getZ());
            		addFields(feat, pt, m_fdef);
            		points.push_back(pt);
                    OGRFeature::DestroyFeature(feat);
            	}
            	m_layer->SetSpatialFilter(NULL);
            }

            void getPoints(std::vector<geotools::util::Point*> &points,
                    int count, int offset = 0) {
            	if(m_type != GeomType::GTPoint)
            		g_runerr("This dataset is not a point dataset.");
            	long total = m_layer->GetFeatureCount(1);
            	if(offset > total) return;
            	if(count > total - offset) count = total - offset;
            	for(long i = offset; i < offset + count; ++i) {
                	OGRFeature *feat = m_layer->GetFeature(i);
					OGRPoint *opt = (OGRPoint *) feat->GetGeometryRef();
					geotools::util::Point *pt = new Point(opt->getX(), opt->getY(), opt->getZ());
            		addFields(feat, pt, m_fdef);
					points.push_back(pt);
                    OGRFeature::DestroyFeature(feat);
            	}
            }

            /*
            void getNearestPoints(const geotools::util::Point &target, int count,
            		std::vector<std::unique_ptr<geotools::util::Point> > &points) {
            	if(!m_type == GeomType::Point)
            		g_runerr("This dataset is not a point dataset.");
            	double rad = 1.0;
            	std::unordered_set<geotools::util::Point*>
            	while(points.size() < count) {
					m_layer->SetSpatialFilterRect(target.x - rad, target.y - rad, target.x + rad, target.y + rad);
					OGRFeature *feat;
					int gidx = this->m_fdef->GetGeomFieldIndex(m_geomName.c_str());
					while((feat = m_layer->GetNextFeature())) {
						OGRPoint *opt = feat->GetGeometryRef();
						geotools::util::Point *pt = new Point(opt->getX(), opt->getY(), opt->getZ());
						SQLite::addFields(feat, pt);
						points.push_back(pt);
					}
					rad *= 2;
            	}
            	m_layer->SetSpatialFilter(NULL);
            }
			*/

            uint64_t getGeomCount() {
            	return m_layer->GetFeatureCount(1);
            }

            void updateFeature(const std::string &idField, uint64_t id, std::unordered_map<std::string, void*> &values) {
            	std::stringstream ss;
            	ss << "\"" << idField << "\"=" << id;
            	m_layer->SetAttributeFilter(ss.str().c_str());
            	OGRFeature *feat;
            	if((feat = m_layer->GetNextFeature())) {
            	   updateFields(feat, values, m_fdef);
                   if(CPLE_None != m_layer->SetFeature(feat))
                        g_runerr("Failed to update feature on " << m_file << ".");
            	} else {
            		g_runerr("Failed to find feature with ID: " << idField << "=" << id);
            	}
            }

            void deleteFeature(const std::string &idField, uint64_t id) {
                std::stringstream ss;
                ss << "\"" << idField << "\"=" << id;
                m_layer->SetAttributeFilter(ss.str().c_str());
                OGRFeature *feat;
                while((feat = m_layer->GetNextFeature())) {
                    if(CPLE_None != m_layer->DeleteFeature(feat->GetFID()))
                        g_runerr("Failed to delete feature with ID: " << idField << "=" << id << ". Deletion may not be implemented for this file type.");
                }
            }

            /*
            void execute(std::string &sql) {
                char *err;
                if(SQLITE_OK != sqlite3_exec(m_db, sql.c_str(), NULL, NULL, &err))
                    handleError("Failed to update: ", err);
            }
			*/

            void begin() {
                if(CPLE_None != m_layer->StartTransaction())
                    g_warn("Failed to start transaction.");
            }

            void rollback() {
                if(CPLE_None != m_layer->RollbackTransaction())
                    g_warn("Failed to roll back transaction.");
            }

            void commit() {
                if(CPLE_None != m_layer->CommitTransaction())
                    g_warn("Failed to commit transaction.");
                m_layer->SyncToDisk();
            }

            int srid() {
                return m_srid;
            }

        };

    } // db

} // geotools

#endif

