#include "db.hpp"

using namespace geotools::db;

namespace geotools {

	namespace db {

		namespace util {

			GeomType geomType(OGRwkbGeometryType type) {
				switch(type) {
				case wkbPoint: return GeomType::GTPoint;
				case wkbLineString: return GeomType::GTLine;
				case wkbPolygon: return GeomType::GTPolygon;
				case wkbMultiPoint: return GeomType::GTMultiPoint;
				case wkbMultiLineString: return GeomType::GTMultiLine;
				case wkbMultiPolygon: return GeomType::GTMultiPolygon;
				default: return GeomType::GTUnknown;
				}
			}

			OGRwkbGeometryType geomType(GeomType type) {
				switch(type) {
				case GeomType::GTPoint: return wkbPoint;
				case GeomType::GTLine: return wkbLineString;
				case GeomType::GTPolygon: return wkbPolygon;
				case GeomType::GTMultiPoint: return wkbMultiPoint;
				case GeomType::GTMultiLine: return wkbMultiLineString;
				case GeomType::GTMultiPolygon: return wkbMultiPolygon;
				default: return wkbUnknown;
				}
			}

			FieldType fieldType(OGRFieldType type) {
				switch(type) {
				case OFTInteger:
				case OFTInteger64: return FieldType::FTInt;
				case OFTString: return FieldType::FTString;
				case OFTReal: return FieldType::FTDouble;
				case OFTBinary: return FieldType::FTBlob;
				default: return FieldType::FTUnknown;
				}
			}

			OGRFieldType fieldType(FieldType type) {
				switch(type) {
				case FieldType::FTInt: return OFTInteger64;
				case FieldType::FTString: return OFTString;
				case FieldType::FTDouble: return OFTReal;
				case FieldType::FTBlob: return OFTBinary;
				default:
					g_argerr("Unknown or unimplemented type: " << type);
				}
			}

			void addFields(OGRFeature *feat, geotools::util::Point *pt, OGRFeatureDefn *ftdef) {
				for(int i = 0; i < feat->GetFieldCount(); ++i) {
					OGRFieldDefn *fdef = ftdef->GetFieldDefn(i);
					OGRField *fld = feat->GetRawFieldRef(i);
					switch(fdef->GetType()) {
					case OFTInteger:
					case OFTInteger64:
						pt->fields[fdef->GetNameRef()] = fld->Integer64;
						break;
					case OFTReal:
						pt->fields[fdef->GetNameRef()] = fld->Real;
						break;
					case OFTString:
						pt->fields[fdef->GetNameRef()] = fld->String;
						break;
					case OFTBinary:
						g_runerr("Binary fields not implemented.");
						break;
					default:
						g_runerr("Field type " << fdef->GetType() << " not implemented.");
						break;
					}
				}
			}

			void updateFields(OGRFeature *feat, std::unordered_map<std::string, void*> &values, OGRFeatureDefn *ftdef) {
				for(int i = 0; i < feat->GetFieldCount(); ++i) {
					OGRFieldDefn *fdef = ftdef->GetFieldDefn(i);
					OGRField *fld = feat->GetRawFieldRef(i);
					std::string *v;
					switch(fdef->GetType()) {
					case OFTInteger:
					case OFTInteger64:
						fld->Integer64 = *((int *) values[fdef->GetNameRef()]);
					std::cerr << fdef->GetNameRef() << ", " << fld->Integer64 << ", " << values[fdef->GetNameRef()] << "\n";
						break;
					case OFTReal:
						fld->Real = *((double *) values[fdef->GetNameRef()]);
					std::cerr << fdef->GetNameRef() << ", " << fld->Real << ", " << values[fdef->GetNameRef()] << "\n";
						break;
					case OFTString:
						v = (std::string *) values[fdef->GetNameRef()];
						std::strcpy(fld->String, v->c_str());
					std::cerr << fdef->GetNameRef() << ", " << fld->String << ", " << values[fdef->GetNameRef()] << "\n";
						break;
					case OFTBinary:
						g_runerr("Binary fields not implemented.");
						break;
					default:
						g_runerr("Field type " << fdef->GetType() << " not implemented.");
						break;
					}
				}
			}

			bool isRast(GDALDriver *drv) {
				const char* cc = drv->GetMetadataItem(GDAL_DCAP_RASTER);
				return cc != NULL && std::strncmp(cc, "YES", 3) == 0;
			}

		} // util
	} // db
} // geotools

using namespace geotools::db::util;

DB::DB(const std::string &file, const std::string &layer, const std::string &driver,
		const std::unordered_map<std::string, FieldType> &fields,
		GeomType type, int srid, bool replace) :
    m_type(type),
    m_srid(srid),
    m_file(file),
	m_layerName(layer),
	m_driver(driver),
	m_fieldTypes(fields),
	m_ds(nullptr),
	m_layer(nullptr),
	m_fdef(nullptr) {

    if(replace && Util::exists(file))
        Util::rm(file);

    GDALAllRegister();

    GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(m_driver.c_str());
    if(!drv)
        g_runerr("Driver not found for " << m_file << " (" << driver << ")");

    m_ds = drv->Create(m_file.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if(!m_ds)
        g_runerr("Failed to create data set for " << m_file);

	if(m_layerName.empty())
		m_layerName = "data";

	OGRSpatialReference sr;
	sr.importFromEPSG(m_srid);

	char **options = nullptr;

	m_layer = m_ds->CreateLayer(m_layerName.c_str(), &sr, geomType(m_type), options);
	if(!m_layer)
		g_runerr("Failed to create layer, " << m_layerName << ".");

	for(const auto &it : m_fieldTypes) {
		OGRFieldDefn def(it.first.c_str(), fieldType(it.second));
		m_layer->CreateField(&def);
	}
	m_fdef = m_layer->GetLayerDefn();

    OGRGeomFieldDefn *gdef = m_layer->GetLayerDefn()->GetGeomFieldDefn(0);
    m_geomName = std::string(gdef->GetNameRef());
}

DB::DB(const std::string &file, const std::string &layer) :
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

    if(!m_layerName.empty())
		m_layer = m_ds->GetLayerByName(m_layerName.c_str());
    if(!m_layer)
		m_layer = m_ds->GetLayer(0);
	if(!m_layer) {
		if(m_layerName.empty()) {
			g_runerr("No layer, " << m_layerName << " was found on this data set, and no default was available.");
		} else {
			g_runerr("No layer was found on this data set.");
		}
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

DB::~DB() {
	#pragma omp critical(__db_destruct)
	{
		if(m_ds)
			GDALClose(m_ds);
		m_ds = nullptr;
	}
}

std::map<std::string, std::set<std::string> > DB::extensions() {
	GDALAllRegister();
	std::map<std::string, std::set<std::string> > extensions;
	GDALDriverManager* mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver* drv = mgr->GetDriver(i);
		if(!isRast(drv)) {
			const char* desc = drv->GetDescription();
			if(desc != NULL) {
				const char *ext = drv->GetMetadataItem(GDAL_DMD_EXTENSION);
				if(ext != NULL ) {
					std::list<std::string> lst;
					Util::splitString(std::string(ext), lst);
					for(const std::string &item : lst)
						extensions[desc].insert(Util::lower(item));
				}
			}
		}
	}
	return extensions;
}

std::map<std::string, std::string> DB::drivers() {
	GDALAllRegister();
	std::map<std::string, std::string> drivers;
	GDALDriverManager *mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver *drv = mgr->GetDriver(i);
		if(!isRast(drv)) {
			const char* name = drv->GetMetadataItem(GDAL_DMD_LONGNAME);
			const char* desc = drv->GetDescription();
			if(name != NULL && desc != NULL)
				drivers[desc] = name;
		}
	}
	return drivers;
}

std::string DB::getDriverForFilename(const std::string &filename) {
	std::string ext = Util::extension(filename);
	std::map<std::string, std::set<std::string> > drivers = extensions();
	std::string result;
	for(const auto &it : drivers) {
		if(it.second.find(ext) != it.second.end())
			result = it.first;
	}
	return result;
}

void DB::clear() {
	g_runerr("Not implemented.");
}

void DB::addPoint(double x, double y, double z, const std::unordered_map<std::string, std::string> &fields) {
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

void DB::addPoints(std::vector<geotools::util::Point*> &points) {
    for (const geotools::util::Point *pt : points)
        addPoint(pt->x, pt->y, pt->z, pt->fields);
}


void DB::setCacheSize(size_t size) {
	g_runerr("Not implemented.");
}

void DB::dropGeomIndex() {
	g_runerr("Not implemented.");
}

void DB::createGeomIndex() {
	g_runerr("Not implemented.");
}

void DB::getPoints(std::vector<geotools::util::Point*> &points,
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

void DB::getPoints(std::vector<geotools::util::Point*> &points,
        int count, int offset) {
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

void DB::getNearestPoints(const geotools::util::Point &target, int count,
		std::vector<std::unique_ptr<geotools::util::Point> > &points) {
	g_runerr("Not implemented.");
}

uint64_t DB::getGeomCount() const {
	return m_layer->GetFeatureCount(1);
}

void DB::updateFeature(const std::string &idField, uint64_t id, std::unordered_map<std::string, void*> &values) {
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

void DB::deleteFeature(const std::string &idField, uint64_t id) {
    std::stringstream ss;
    ss << "\"" << idField << "\"=" << id;
    m_layer->SetAttributeFilter(ss.str().c_str());
    OGRFeature *feat;
    while((feat = m_layer->GetNextFeature())) {
        if(CPLE_None != m_layer->DeleteFeature(feat->GetFID()))
            g_runerr("Failed to delete feature with ID: " << idField << "=" << id << ". Deletion may not be implemented for this file type.");
    }
}

void DB::execute(std::string &sql) {
	g_runerr("Not implemented.");
}

void DB::begin() {
    if(CPLE_None != m_layer->StartTransaction())
        g_warn("Failed to start transaction.");
}

void DB::rollback() {
    if(CPLE_None != m_layer->RollbackTransaction())
        g_warn("Failed to roll back transaction.");
}

void DB::commit() {
    if(CPLE_None != m_layer->CommitTransaction())
        g_warn("Failed to commit transaction.");
    m_layer->SyncToDisk();
}

int DB::srid() const {
    return m_srid;
}

