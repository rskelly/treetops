#include "db.hpp"

using namespace geotools::db;

GeomType geotools::db::util::geomType(OGRwkbGeometryType type) {
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

OGRwkbGeometryType geotools::db::util::geomType(GeomType type) {
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

FieldType geotools::db::util::fieldType(OGRFieldType type) {
    switch(type) {
    case OFTInteger:
    case OFTInteger64: return FieldType::FTInt;
    case OFTString: return FieldType::FTString;
    case OFTReal: return FieldType::FTDouble;
    case OFTBinary: return FieldType::FTBlob;
    default: return FieldType::FTUnknown;
    }
}

OGRFieldType geotools::db::util::fieldType(FieldType type) {
    switch(type) {
    case FieldType::FTInt: return OFTInteger64;
    case FieldType::FTString: return OFTString;
    case FieldType::FTDouble: return OFTReal;
    case FieldType::FTBlob: return OFTBinary;
    default:
        g_argerr("Unknown or unimplemented type: " << type);
    }
}

void geotools::db::util::addFields(OGRFeature *feat, geotools::util::Point *pt, OGRFeatureDefn *ftdef) {
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

void geotools::db::util::updateFields(OGRFeature *feat, std::unordered_map<std::string, void*> &values, OGRFeatureDefn *ftdef) {
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
