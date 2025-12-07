/*
 * treetops
 * 
 * This library provides methods for isolating tree tops and crowns from a
 * LiDAR (or other) canopy height model (CHM.) The output from this program is 
 * ideal for use with the spectral extraction module (spectral) or other 
 * analysis.
 * 
 * The usual sequence for producing useful output is:
 * 1) Smooth the original CHM. The CHM should possibly have pit- or 
 *    spike-removal applied to it first. Smoothing is a Gaussian kernel with
 *    configurable sigma and window size.
 * 2) Locate treetops. This is performed on the smoothed CHM and uses a 
 *    maximum-value kernel to locate local maxima.
 * 3) Delineate crowns. Uses the tree tops as seeds in a region-growing 
 *    algorithm that creates tree crown boundaries with configurable limits.
 * 4) Optional: Locates the "actual" tree top height from the original CHM.
 *              Polygonizes the crowns; transfers data from treetops table to
 *              crowns table.
 *              Cleans up polygons.
 *
 *  Created on: May 3, 2016
 *      Author: Rob Skelly 
 *       Email: rob@dijital.ca
 */

#include "geo.hpp"

#ifdef _WIN32
#include <Windows.h>
#include <io.h>
typedef long ssize_t;
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#define O_BINARY 0	// Not available in linux because would make no sense (text == binary).
#endif

#include <queue>
#include <iostream>
#include <fstream>
#include <atomic>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <vector>

#include <sqlite3.h>
#include <json/json.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>

#include "ds/simple_interval_tree.hpp"
#include "ds/mqtree.hpp"
#include "crowns.hpp"
#include "grid.hpp"
#include "util.hpp"

using namespace dijital;
using namespace dijital::grid;
using namespace dijital::util;
using namespace dijital::ds;
using namespace dijital::crowns;
using namespace dijital::crowns::util;
using namespace dijital::crowns::config;

namespace {

	/**
	 * Returns the extension for a database driver name. Anything other than
	 * 'sqlite' or 'esri shapefile' (case irrelevant) raises an argument error.
	 *
	 * \param driver The driver name.
	 */
	std::string dbExt(const std::string& driver) {
		std::string d = lowercase(driver);
		if(d == "sqlite") {
			return ".sqlite";
		} else if(d== "esri shapefile") {
			return ".shp";
		} else {
			g_argerr("Unknown database driver: " << driver);
		}
	}

	/**
	 * Returns the extension for a raster driver name. Anything other than
	 * 'gtiff' or 'hfa' (case irrelevant) raises an argument error.
	 *
	 * \param driver The driver name.
	 */
	std::string rastExt(const std::string& driver) {
		std::string d = lowercase(driver);
		if(d == "gtiff") {
			return ".tif";
		} else if(d == "hfa") {
			return ".img";
		} else {
			g_argerr("Unknown raster driver: " << driver);
		}
	}

	/**
	 * Replaces the file extension of the given filename with a new extension. The
	 * period is included in the new extension.
	 *
	 * \param filename The original filename.
	 * \param ext The new extension.
	 */
	void replaceExt(std::string& filename, const std::string& ext) {
		std::string ext0 = extension(filename);
		filename.replace(filename.find(ext0), ext0.size(), ext);
	}

	/**
	 * Updates the given strings with filenames appropriate to each output type,
	 * derived from the first argument, the source file name. If the source file
	 * doesn't exist, raises an argument error.
	 *
	 * \param filename The template filename.
	 * \param smoothed The name of the smoothed file.
	 * \param smoothedDriver The name of the driver for the smoothed file.
	 * \param topsDb The name of the treetops database file.
	 * \param topsDriver The name of the driver for the treetops database.
	 * \param crownsRast The name of the crowns raster file.
	 * \param crownsRasterDriver The name of the crowns raster driver.
	 * \param crownsDb The name of the crowns database file.
	 * \param crownsDbDriver The name of the crowns database driver.
	 * \param settings The name of the settings file.
	 */
	void formatFilenames(const std::string& filename,
			std::string& smoothed, const std::string& smoothedDriver,
			std::string& topsDb, const std::string& topsDriver,
			std::string& crownsRast, const std::string& crownsRastDriver,
			std::string& crownsDb, const std::string& crownsDbDriver,
			std::string& settings) {

		if(!isfile(filename))
			g_argerr("The given CHM doesn't exist: " << filename)

		std::string dir = parent(filename);
		std::string ext = extension(filename);
		std::string base = basename(filename);
		std::string tpl = base.substr(0, base.find(ext));

		smoothed = join(dir, tpl + "_smoothed" + rastExt(smoothedDriver));
		topsDb = join(dir, tpl + "_tops" + dbExt(topsDriver));
		crownsRast = join(dir, tpl + "_crowns" + rastExt(crownsRastDriver));
		crownsDb = join(dir, tpl + "_crowns" + dbExt(crownsDbDriver));
		settings = join(dir, tpl + "_settings.txt");
	}

	/**
	 * \brief Creates and configures a GDAL dataset to contain manage the crowns database.
	 *
	 * \param[inout] ds A reference to the GDALDataset pointer.
	 * \param[inout] layer A reference to the OGRLayer pointer.
	 * \param filename The database filename. Will be deleted and recreated if it exists.
	 * \param driver The database driver.
	 * \param layerName The name of the table. Manifestation depends on the file format.
	 * \param projection A WKT projection string.
	 * \param type The geometry type.
	 * \param idField The name of the identifier field. Must be represented in the fields parameter.
	 * \param fields A list of field names and types to create on the dataset.
	 */
	void makeCrownDataset(GDALDataset*& ds, OGRLayer*& layer,
			const std::string& filename, const std::string& driver,
			const std::string& layerName,
			const std::string& projection, OGRwkbGeometryType type,
			const std::string& idField, const std::vector<std::pair<std::string, OGRFieldType>>& fields) {

		ds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, nullptr, nullptr, nullptr);

		std::string drvl = lowercase(driver);

		if(!ds) {
			GDALAllRegister();
			// Get the vector driver.
			GDALDriverManager* dm = GetGDALDriverManager();
			GDALDriver* drv = dm->GetDriverByName(driver.c_str());
			if(!drv) {
				int dcnt = dm->GetDriverCount();
				std::string ext2(lowercase(extension(filename).substr(1)));
				for(int i = 0; i < dcnt; ++i) {
					GDALDriver* d = dm->GetDriver(i);
					const char* ext = GDALGetMetadataItem(d, GDAL_DMD_EXTENSION, "");
					if(!ext)
						continue;
					std::string ext1(ext);
					if(lowercase(ext1) == ext2) {
						drv = d;
						break;
					}
				}
			}
			if(!drv)
				g_runerr("Failed to find driver for " << driver << ".");

			// Create an output dataset for the polygons.
			char** dopts = NULL;
			if(drvl == "spatialite" || drvl == "sqlite") {
				drvl = "sqlite";
				//dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
			}

			ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
			CSLDestroy(dopts);
			if(!ds)
				g_runerr("Failed to create dataset " << filename << ".");

		}

		// If a projection is given, create a spatial reference object for it.
		OGRSpatialReference* sr = nullptr;
		if(!projection.empty())
			sr = new OGRSpatialReference(projection.c_str());

		// Create the layer.
		char** lopts = NULL;
		if(drvl == "sqlite") {
			lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
		} else if(drvl == "esri shapefile") {
			lopts = CSLSetNameValue(lopts, "2GB_LIMIT", "YES");
		}

		layer = ds->CreateLayer(layerName.c_str(), sr, type, lopts);
		CSLDestroy(lopts);

		if(!layer)
			g_runerr("Failed to create layer " << layerName << ".");

		// Create the field set for the database.
		for(const std::pair<std::string, OGRFieldType>& field : fields) {
			OGRFieldDefn dfn(field.first.c_str(), field.second);
			layer->CreateField(&dfn, field.first == idField); // Flag as the key if it matches idField.
		}

		// Dispose of the spatial reference object.
		if(sr)
			sr->Release();
	}

	/**
	 * \brief Set the z coordinate on a multipolygon, recursively.
	 * \param gctx The GEOS context handle.
	 * \param geom The geometry.
	 * \param z The height to set.
	 * \return A new geometry. Input and output are reposibility of the caller.
	 */
	GEOSGeometry* setGeomZ(GEOSContextHandle_t gctx, const GEOSGeometry* geom, float z) {
		int type = GEOSGeomTypeId_r(gctx, geom);
		if(type == GEOS_MULTIPOLYGON) {
			std::vector<GEOSGeometry*> geoms;
			for(int i = 0; i < GEOSGetNumGeometries_r(gctx, geom); ++i)
				geoms.push_back(setGeomZ(gctx, GEOSGetGeometryN_r(gctx, geom, i), z));
			return GEOSGeom_createCollection_r(gctx, GEOS_MULTIPOLYGON, geoms.data(), geoms.size());
		} else if(type == GEOS_POLYGON) {
			const GEOSGeometry* lr = GEOSGetExteriorRing_r(gctx, geom);
			GEOSGeometry* lr0 = setGeomZ(gctx, lr, z);
			std::vector<GEOSGeometry*> irs;
			for(int i = 0; i < GEOSGetNumInteriorRings_r(gctx, geom); ++i) {
				const GEOSGeometry* ir = GEOSGetInteriorRingN_r(gctx, geom, i);
				irs.push_back(setGeomZ(gctx, ir, z));
			}
			return GEOSGeom_createPolygon_r(gctx, lr0, irs.data(), irs.size());
		} else if(type == GEOS_LINEARRING) {
			const GEOSCoordSequence* cs = GEOSGeom_getCoordSeq_r(gctx, geom);
			unsigned int size;
			GEOSCoordSeq_getSize_r(gctx, cs, &size);
			GEOSCoordSequence* cs0 = GEOSCoordSeq_create_r(gctx, size, 3);
			double x, y;
			for(unsigned i = 0; i < size; ++i) {
				GEOSCoordSeq_getX_r(gctx, cs, i, &x);
				GEOSCoordSeq_getY_r(gctx, cs, i, &y);
				GEOSCoordSeq_setX_r(gctx, cs0, i, x);
				GEOSCoordSeq_setY_r(gctx, cs0, i, y);
				GEOSCoordSeq_setZ_r(gctx, cs0, i, z);
			}
			return GEOSGeom_createLinearRing_r(gctx, cs0);
		} else{
			g_runerr("Unexpected geom type: " << type);
		}
	}

	bool hasMember(Json::Value& v, const std::string& name) {
		for(const auto& id : v.getMemberNames()) {
			if(name == id)
				return true;
		}
		return false;
	}

	template<class T>
	T readJSON(Json::Value& v, const std::string& name, T def) {
		try {
			if(hasMember(v, name))
				return v[name].as<T>();
		} catch(...) {}
		return def;
	}

} // anon


// CrownDB implementation

// Size for read/write buffer.
static const size_t CDBBUFSIZE = sizeof(size_t) * 3 + sizeof(Top);

CrownDB::CrownDB() :
		m_tfd(0), m_cfd(0),
		m_toffset(0), m_coffset(0),
		m_tidx(0), m_count(0) {

	m_gctx = GEOS_init_r();
	m_rdr = GEOSWKBReader_create_r(m_gctx);
	m_wtr = GEOSWKBWriter_create_r(m_gctx);

	m_tfile = tmpfile("tt");
	m_cfile = tmpfile("tt");
	if((m_tfd = open(m_tfile.c_str(), O_CREAT|O_RDWR|O_BINARY, 0777)) <= 0)
		g_runerr("Failed to open tops DB file at " << m_tfile);
	if((m_cfd = open(m_cfile.c_str(), O_CREAT|O_RDWR|O_BINARY, 0777)) <= 0)
		g_runerr("Failed to open tops DB file at " << m_cfile);
}

GEOSContextHandle_t CrownDB::gctx() {
	return m_gctx;
}

bool CrownDB::next(Top& top) {
	static CrownGeom geom;
	return next(top, geom, false);
}

bool CrownDB::next(Top& top, std::vector<unsigned char>& crownwkb, bool loadCrown) {
	size_t id;
	size_t plen;
	size_t ppos;
	size_t pos = 0;
	std::vector<unsigned char> buf(CDBBUFSIZE);

	// Loop until a next one is found or the file runs out.
	while(true) {
		// Go to start of line at offset.
		if (-1 == lseek(m_tfd, m_tidx * CDBBUFSIZE, SEEK_SET)) {
			g_warn("Failed to seek in tops DB file: " << strerror(errno));
			return false;
		}
		// Read the line into the buffer.
		int r;
		if ((r = read(m_tfd, buf.data(), CDBBUFSIZE)) < (ssize_t)CDBBUFSIZE) {
			// The file may be finished; not an error.
			if(r != 0)
				g_warn("Failed to read tops DB file (this may not be an error): " << r << "; " << strerror(errno));
			return false;
		}
		// Read the values.
		std::memcpy(&id, buf.data(), sizeof(size_t));			pos += sizeof(size_t);
		std::memcpy(&top, buf.data() + pos, sizeof(Top));		pos += sizeof(Top);
		std::memcpy(&plen, buf.data() + pos, sizeof(size_t));	pos += sizeof(size_t);
		// If there's a poly, the plen will be > 0. Read it. Otherwise crown is null.
		if(loadCrown && plen > 0) {
			std::memcpy(&ppos, buf.data() + pos, sizeof(size_t));
			crownwkb.resize(plen);

			// Go to the poly.
			if(-1 == lseek(m_cfd, ppos, SEEK_SET)) {
				g_warn("Failed to seek in crowns DB file: " << strerror(errno));
				return false;
			}
			if((r = read(m_cfd, crownwkb.data(), plen)) < (ssize_t) plen) {
				g_warn("Failed to read " << plen << " bytes from crowns DB file (this is an error): " << r << "; " << strerror(errno));
				return false;
			}

		}

		++m_tidx;
		return true;
	}
}

bool CrownDB::next(Top& top, CrownGeom& crown, bool loadCrown) {
	static std::vector<unsigned char> cbuf;
	bool ret = next(top, cbuf, loadCrown);
	if(ret && loadCrown) {
		if(!crown.gctx)
			crown.gctx = m_gctx;
		if(!(crown.crown = GEOSWKBReader_read_r(crown.gctx, m_rdr, cbuf.data(), cbuf.size())))
			return false;
	}
	return ret;
}

void CrownDB::reset() {
	m_tidx = 0;
}

bool CrownDB::insert(const Top& top) {
	CrownGeom crown;
	return insert(top, crown);
}

bool CrownDB::insert(const Top& top, CrownGeom& crown) {

	// If it exists, update it.
	if(update(top, crown))
		return true;

	size_t pos = 0;
	size_t plen = 0;
	size_t ppos = 0;
	std::vector<unsigned char> buf(CDBBUFSIZE);

	// Write the ID, top.
	std::memcpy(buf.data(), &top.id, sizeof(size_t));		pos += sizeof(size_t);
	std::memcpy(buf.data() + pos, &top, sizeof(Top));		pos += sizeof(Top);

	// Write the geometry.
	if(crown.crown) {
		if(!crown.gctx)
			crown.gctx = m_gctx;
		if(!writeCrown(plen, ppos, crown))
			return false;
	}

	// Write back the crown length and position.
	std::memcpy(buf.data() + pos, &plen, sizeof(size_t));	pos += sizeof(size_t);
	std::memcpy(buf.data() + pos, &ppos, sizeof(size_t));

	// Write the Top record.
	{
		std::lock_guard<std::mutex> lk(m_tmtx);
		if(-1 == lseek(m_tfd, m_toffset, SEEK_SET))
			return false;
		if(write(m_tfd, buf.data(), CDBBUFSIZE) < (ssize_t) CDBBUFSIZE)
			return false;
	}

	{
		std::lock_guard<std::mutex> lk(m_imtx);
		m_idx[top.id] = m_toffset;
		m_toffset += CDBBUFSIZE;
		++m_count;
	}

	return true;
}

bool CrownDB::find(int id, Top& top, CrownGeom& crown, bool withCrown) {

	if(id <= 0)
		return false;

	size_t offset = -1;
	{
		std::lock_guard<std::mutex> lk(m_imtx);
		if(m_idx.find(id) != m_idx.end())
			offset = m_idx.at(id);
	}

	if(offset == (size_t) -1)
		return false;

	std::vector<unsigned char> buf(CDBBUFSIZE);

	{
		std::lock_guard<std::mutex> lk(m_tmtx);
		// Go to start of line at offset.
		if(-1 == lseek(m_tfd, offset, SEEK_SET))
			return false;
		// Read the size.
		if(read(m_tfd, buf.data(), CDBBUFSIZE) < (ssize_t) CDBBUFSIZE)
			return false;
	}

	size_t id_;
	size_t plen;
	size_t ppos;
	size_t pos = 0;

	// Read the values.
	std::memcpy(&id_, buf.data(), sizeof(size_t));			pos += sizeof(size_t);
	std::memcpy(&top, buf.data() + pos, sizeof(Top));			pos += sizeof(Top);
	std::memcpy(&plen, buf.data() + pos, sizeof(size_t));		pos += sizeof(size_t);

	// If there's a poly, the plen will be > 0. Read it. Otherwise crown is null.
	if(withCrown && plen > 0) {
		if(!crown.gctx)
			crown.gctx = m_gctx;
		std::memcpy(&ppos, buf.data() + pos, sizeof(size_t));
		std::vector<unsigned char> pbuf(plen);

		std::lock_guard<std::mutex> lk(m_cmtx);
		// Go to the poly.
		if(-1 == lseek(m_cfd, ppos, SEEK_SET))
			return false;
		if(read(m_cfd, pbuf.data(), plen) < (ssize_t) plen)
			return false;
		if(!(crown.crown = GEOSWKBReader_read_r(crown.gctx, m_rdr, pbuf.data(), plen)))
			return false;
	} 

	return true;
}

int CrownDB::count() {
	return m_count;
}

bool CrownDB::writeCrown(size_t& plen, size_t& ppos, CrownGeom& crown) {
	bool ret = false;
	if(crown.crown) {
		if(!crown.gctx)
			crown.gctx = m_gctx;

		// Write the geometry to a char buffer.
		size_t plen0;
		unsigned char* poly = GEOSWKBWriter_write_r(crown.gctx, m_wtr, crown.crown, &plen0);

		// If the polygon is a different size, have to write it at the end of the crown file.
		if(plen0 > plen) {
			ppos = m_coffset;
			m_coffset += plen0;
		}
		plen = plen0;

		if(plen) {
			std::lock_guard<std::mutex> lk(m_cmtx);
			// Seek and write.
			if (lseek(m_cfd, ppos, SEEK_SET) == (ssize_t)ppos) {
				if (write(m_cfd, poly, plen) == (ssize_t)plen)
					ret = true;
			}
		}

		GEOSFree_r(crown.gctx, poly);
	}
	return ret;
}

bool CrownDB::update(const Top& top) {
	CrownGeom crown;
	return update(top, crown);
}

bool CrownDB::update(const Top& top, CrownGeom& crown) {

	size_t offset = -1;

	if(top.id > 0) {
		std::lock_guard<std::mutex> lk(m_imtx);
		if(m_idx.find(top.id) != m_idx.end())
			offset = m_idx.at(top.id);
	}

	if(offset == (size_t) -1)
		return false;

	std::vector<unsigned char> buf(CDBBUFSIZE);
	{
		std::lock_guard<std::mutex> lk(m_tmtx);
		if(-1 == lseek(m_tfd, offset, SEEK_SET))
			return false;
		if(read(m_tfd, buf.data(), CDBBUFSIZE) < (ssize_t) CDBBUFSIZE)
			return false;
	}

	// Get the poly length and position first.
	size_t plen;
	size_t ppos;
	size_t pos = sizeof(size_t) + sizeof(Top); // Skip the ID and Top.
	std::memcpy(&plen, buf.data() + pos, sizeof(size_t));		pos += sizeof(size_t);
	std::memcpy(&ppos, buf.data() + pos, sizeof(size_t));

	// Write the Top.
	pos = sizeof(size_t);
	std::memcpy(buf.data() + pos, &top, sizeof(Top));			pos += sizeof(Top);

	if(crown.crown) {
		if(!crown.gctx)
			crown.gctx = m_gctx;
		// Write the crown, update length and position in the crown file.
		if(!writeCrown(plen, ppos, crown))
			return false;
	} else {
		// If there's no crown, zero out the position/length.
		ppos = 0;
		plen = 0;
	}

	// Write the original pos/size back.
	std::memcpy(buf.data() + pos, &plen, sizeof(size_t));		pos += sizeof(size_t);
	std::memcpy(buf.data() + pos, &ppos, sizeof(size_t));

	{
		std::lock_guard<std::mutex> lk(m_tmtx);
		// Seek to position and overwrite the ID as zero.
		if(-1 == lseek(m_tfd, offset, SEEK_SET))
			return false;
		if(write(m_tfd, buf.data(), CDBBUFSIZE) < (ssize_t) CDBBUFSIZE)
			return false;
	}

	return true;
}

void CrownDB::saveTops(const std::string& filename, const std::string& driver,
		const std::string& layerName, const std::string& projection) {

	rem(filename);

	std::string idField = "tree_id";
	std::vector<std::pair<std::string, OGRFieldType> > fields = {
		{idField, OFTInteger},
		{"parent_id", OFTInteger},
		{"ground_z", OFTReal},
		{"smoothed_z", OFTReal},
		{"tree_ht", OFTReal},
		{"top_z", OFTReal}
	};

	// Create the output dataset
	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();

	// Create the output dataset
	GDALDataset* ds;
	OGRLayer* layer;

	makeCrownDataset(ds, layer, filename, driver, layerName,
			projection,	wkbPoint25D,
			idField, fields);

	// Start a transaction on the layer.
	if(OGRERR_NONE != layer->StartTransaction())
		g_runerr("Failed to start transaction.");

	Top top;
	GEOSGeometry* geom = nullptr;
	int count = 0;

	reset();

	CrownGeom crown(geom, gctx);

	while(next(top, crown, false)) {
		float x, y, z, th, gz;
		if(top.oz > 0) {
			x = top.ox;
			y = top.oy;
			z = top.oz;
			gz = top.groundZ;
			th = z * 0.5 + gz;
		} else {
			x = top.sx;
			y = top.sy;
			z = top.sz;
			gz = 0;
			th = 0;
		}
		// Note: Fields added in same order as defined.
	    OGRFeature feat(layer->GetLayerDefn());
	    feat.SetField(0, (GIntBig) top.id);
	    feat.SetField(1, (GIntBig) top.parentID);
	    feat.SetField(2, gz);
	    feat.SetField(3, top.sz);
	    feat.SetField(4, th);
	    feat.SetField(5, z);
		OGRPoint geom(x, y, z);
		feat.SetGeometry(&geom);
	    if(OGRERR_NONE != layer->CreateFeature(&feat))
	    	throw CrownDBSaveException("Failed to add feature to " + filename);

	    Monitor::get().status((float) count++ / m_count);
	}

	// Commit and release the layer -- GDAL will take care of it. But close the dataset so that can happen.
	if(OGRERR_NONE != layer->CommitTransaction())
		g_runerr("Failed to commit transaction.");
	layer->Dereference();
	GDALClose(ds);

	OGRGeometry::freeGEOSContext(gctx);

	Monitor::get().status(1.0f);
}

bool CrownDB::checkSaveCrowns(const std::string& driver) {
	size_t gb = 2L * 1024L * 1024L * 1024L;
	// This is an empirically-derived estimate for the filesize.
	// Crowns are variable-length and this may still fail in the save step.
	if(lowercase(driver) == "esri shapefile" && 2 * crownsDBSize() > gb)
		return false;
	return true;
}

bool CrownDB::checkSaveTops(const std::string& driver) {
	size_t gb = 2L * 1024L * 1024L * 1024L;
	// This is an empirically-derived estimate for the filesize with a buffer.
	// Unlike crowns, Tops are fixed-length.
	if(lowercase(driver) == "esri shapefile" && 0.35 * topsDBSize() > gb)
		return false;
	return true;
}

void CrownDB::saveCrowns(const std::string& filename, const std::string& driver,
		const std::string& layerName, const std::string& projection) {

	rem(filename);

	std::string idField = "tree_id";
	std::vector<std::pair<std::string, OGRFieldType> > fields = {
		{idField, OFTInteger},
		{"parent_id", OFTInteger},
		{"ground_z", OFTReal},
		{"smoothed_z", OFTReal},
		{"tree_ht", OFTReal},
		{"top_z", OFTReal}
	};

	// Create the output dataset
	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();

	// Create the output dataset
	GDALDataset* ds;
	OGRLayer* layer;

	makeCrownDataset(ds, layer, filename, driver, layerName,
			projection,	wkbMultiPolygon25D,
			idField, fields);

	reset();

	// Start a transaction on the layer.
	if(OGRERR_NONE != layer->StartTransaction())
		g_runerr("Failed to start transaction.");

	Top top;
	GEOSGeometry* geom = nullptr;
	int count = 0;

	reset();

	CrownGeom crown(geom, gctx);
	
	while(next(top, crown, true)) {
		if(crown.crown == nullptr) {
			// g_warn("Null geometry.");
			continue;
		}
		float z, th, gz;
		if(top.oz > 0) {
			z = top.oz;
			gz = top.groundZ;
			th = z * 0.5 + gz;
		} else {
			z = top.sz;
			gz = 0;
			th = 0;
		}

		// Set the z-component on the geometry and all its children.
		GEOSGeometry* geom0 = setGeomZ(gctx, crown.crown, th);
		GEOSGeom_destroy_r(gctx, crown.crown);

		// Note: Fields added in same order as defined.
	    OGRFeature feat(layer->GetLayerDefn());
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) geom0);
	    feat.SetField(0, (GIntBig) top.id);
	    feat.SetField(1, (GIntBig) top.parentID);
	    feat.SetField(2, gz);
	    feat.SetField(3, top.sz);
	    feat.SetField(4, th);
	    feat.SetField(5, z);
		feat.SetGeometry(ogeom);
		OGRErr res = layer->CreateFeature(&feat);
		OGRGeometryFactory::destroyGeometry(ogeom);
		GEOSGeom_destroy_r(gctx, geom0);
	    if(res != OGRERR_NONE)
	    	throw CrownDBSaveException("Failed to add feature to " + filename);

	    Monitor::get().status((float) count++ / m_count);
	}

	// Commit and release the layer -- GDAL will take care of it. But close the dataset so that can happen.
	if(OGRERR_NONE != layer->CommitTransaction())
		g_runerr("Failed to commit transaction.");
	layer->Dereference();
	GDALClose(ds);

	OGRGeometry::freeGEOSContext(gctx);

	Monitor::get().status(1.0f);
}

size_t CrownDB::crownsDBSize() {
	return filesize(m_cfile);
}

size_t CrownDB::topsDBSize() {
	return filesize(m_tfile);
}

CrownDB::~CrownDB() {
	close(m_tfd);
	close(m_cfd);
	rem(m_tfile);
	rem(m_cfile);
	GEOSWKBReader_destroy_r(m_gctx, m_rdr);
	GEOSWKBWriter_destroy_r(m_gctx, m_wtr);
	GEOS_finish_r(m_gctx);
}

// CrownsAppConfig implementation

CrownsAppConfig::CrownsAppConfig() :
	m_srid(0),
	m_buildIndex(false),
	m_tableCacheSize(1024 * 1024),
	m_rowCacheSize(24 * 1024 * 1024),
	m_threads(1),
	m_doSmoothing(false),
	m_doTops(false),
	m_doCrowns(false),
	m_originalCHMBand(1),
	m_smoothWindowSize(3),
	m_smoothSigma(0.8),
	m_topsMaxNulls(0.2),
	m_crownsUpdateHeights(true),
	m_crownsDoDatabase(false),
	m_crownsRemoveHoles(false),
	m_crownsRemoveDangles(false),
	m_crownsKeepSmoothed(false),
	m_listener(nullptr),
	m_active(false),
	m_locked(false),
	m_smoothExisted(false) {

	reset();
}

void CrownsAppConfig::loadFromJSON(const std::string& config_file) {
	// Load the config.
	std::ifstream f(config_file);
	Json::Value root;
	f >> root;
	m_srid = readJSON(root["srid"], "value", 4326);
	m_originalCHMBand = readJSON(root["originalCHMBand"], "value", 1);
	m_originalCHM = readJSON(root["originalCHM"], "value", "");
	m_smoothedCHM = readJSON(root["smoothedCHM"], "value", "");
	m_topsDatabase = readJSON(root["topsDatabase"], "value", "");
	m_topsDatabaseDriver = readJSON(root["topsDatabaseDriver"], "value", "");
	m_crownsRaster = readJSON(root["crownsRaster"], "value", "");
	m_crownsRasterDriver = readJSON(root["crownsRasterDriver"], "value", "");
	m_crownsDatabase = readJSON(root["crownsDatabase"], "value", "");
	m_crownsDatabaseDriver = readJSON(root["crownsDatabaseDriver"], "value", "");

	m_buildIndex = readJSON(root["buildIndex"], "value", false);
	m_tableCacheSize = readJSON(root["tableCacheSize"], "value", 1024 * 1024);
	m_rowCacheSize = readJSON(root["rowCacheSize"], "value", 24 * 1024 * 1024);
	m_threads = readJSON(root["threads"], "value", 1);
	m_doSmoothing = readJSON(root["doSmoothing"], "value", true);
	m_doTops = readJSON(root["doTops"], "value", false);
	m_doCrowns = readJSON(root["doCrowns"], "value", false);
	m_smoothWindowSize = readJSON(root["smoothingWindowSize"], "value", 3);
	m_smoothSigma = readJSON(root["smoothSigma"], "value", 0.8);
	m_topsMaxNulls = readJSON(root["topsMaxNulls"], "value", 0.2);
	m_crownsUpdateHeights = readJSON(root["crownsUpdateHeights"], "value", true);
	m_crownsDoDatabase = readJSON(root["crownsDoDatabase"], "value", false);
	m_crownsRemoveHoles = readJSON(root["crownsRemoveHoles"], "value", false);
	m_crownsRemoveDangles = readJSON(root["crownsRemoveDangles"], "value", false);
	m_crownsKeepSmoothed = readJSON(root["crownsKeepSmoothed"], "value", false);
	for(Json::Value& tt : root["topsThresholds"]) {
		float t = readJSON(tt["threshold"], "value", 4);
		int w = readJSON(tt["window"], "value", 3);
		m_topsThresholds.emplace_back(t, w);
	}
	for(Json::Value& ct : root["crownsThresholds"]) {
		float t = readJSON(ct["threshold"], "value", 4.0);
		float f = readJSON(ct["fraction"], "value", 0.65);
		float r = readJSON(ct["radius"], "value", 15.0);
		m_crownsThresholds.emplace_back(t, f, r);
	}
}

void CrownsAppConfig::reset() {
	m_crownDb.reset(new CrownDB());
}

CrownDB& CrownsAppConfig::crownDB() {
	return *m_crownDb;
}

Band<float>& CrownsAppConfig::chm() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_chm.get())
			m_chm.reset(new Band<float>(originalCHM(), originalCHMBand() - 1, false, true));
	}
	return *m_chm;
}

Band<uint32_t>& CrownsAppConfig::crowns() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_crowns.get()) {
			GridProps pr(chm().props());
			pr.setBands(1);
			pr.setDriver(crownsRasterDriver());
			pr.setWritable(true);
			pr.setDataType(DataType::UInt32);
			pr.setNoData(0);
			m_crowns.reset(new Band<uint32_t>(crownsRaster(), pr));
		}
	}
	return *m_crowns;
}

Band<float>& CrownsAppConfig::smoothed() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_smoothed.get()) {
			if(!isfile(smoothedCHM())) {
				GridProps pr(chm().props());
				pr.setBands(1);
				pr.setDriver(smoothedCHMDriver());
				pr.setWritable(true);
				pr.setDataType(DataType::Float32);
				pr.setNoData(-9999);
				m_smoothed.reset(new Band<float>(smoothedCHM(), pr));
				m_smoothExisted = false;
			} else {
				m_smoothed.reset(new Band<float>(smoothedCHM(), 0, true, true));
				m_smoothExisted = true;
			}
		}
	}
	return *m_smoothed;
}

bool CrownsAppConfig::smoothExisted() const {
	return m_smoothExisted;
}

void CrownsAppConfig::setSmoothExisted(bool v) {
	m_smoothExisted = v;
}

Band<uint32_t>& CrownsAppConfig::topsWindow() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_topsWindow.get()) {
			GridProps pr = chm().props();
			pr.setDataType(DataType::Byte);
			pr.setNoData(0);
			pr.setWritable(true);
			m_topsWindow.reset(new Band<uint32_t>(pr, true));
		}
	}
	return *m_topsWindow;
}

Band<uint32_t>& CrownsAppConfig::topsID() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_topsID.get()) {
			GridProps pr = chm().props();
			pr.setDataType(DataType::Byte);
			pr.setNoData(0);
			pr.setWritable(true);
			m_topsID.reset(new Band<uint32_t>(pr, true));
		}
	}
	return *m_topsID;
}

Band<uint32_t>& CrownsAppConfig::topsParent() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_topsParent.get()) {
			GridProps pr = chm().props();
			pr.setDataType(DataType::Byte);
			pr.setNoData(0);
			pr.setWritable(true);
			m_topsParent.reset(new Band<uint32_t>(pr, true));
		}
	}
	return *m_topsParent;
}

dijital::ds::mqtree<dijital::crowns::util::Top>& CrownsAppConfig::topsTree() {
	static std::mutex mtx;
	{
		std::lock_guard<std::mutex> lk(mtx);
		if(!m_topsTree.get())
			m_topsTree.reset(new dijital::ds::mqtree<dijital::crowns::util::Top>());
	}
	return *m_topsTree;
}

void CrownsAppConfig::flush() {
	if(m_chm.get())
		m_chm->flush();
	if(m_crowns.get())
		m_crowns->flush();
	if(m_smoothed.get())
		m_smoothed->flush();
}

void CrownsAppConfig::destroy() {
	flush();
	m_chm.reset(nullptr);
	m_crowns.reset(nullptr);
	m_smoothed.reset(nullptr);
	m_topsTree.reset(nullptr);
	m_topsParent.reset(nullptr);
	m_topsID.reset(nullptr);
	m_topsWindow.reset(nullptr);
}

void CrownsAppConfig::lock() {
	m_locked = true;
}

void CrownsAppConfig::unlock() {
	m_locked = false;
}

void CrownsAppConfig::setListener(CrownsAppConfigListener* listener) {
	m_listener = listener;
}

CrownsAppConfigListener* CrownsAppConfig::listener() {
	return m_listener;
}

bool CrownsAppConfig::setActive(bool active) {
	bool orig = m_active;
	m_active = active;
	return orig;
}

bool CrownsAppConfig::active() const {
	return m_active;
}

void CrownsAppConfig::update(long field) {
	if(m_active && m_listener)
		m_listener->configUpdate(*this, field);
}

void CrownsAppConfig::setSettings(const std::string& filename) {
	if(m_locked) return;
	m_settings = filename;
	update(TTSettingsFile);
}

const std::string& CrownsAppConfig::settings() const {
	return m_settings;
}

void CrownsAppConfig::setBounds(const Bounds<float>& bounds) {
	if(m_locked) return;
	m_bounds = bounds;
	update(BoundsField);
}

const Bounds<float>& CrownsAppConfig::bounds() const {
	return m_bounds;
}

std::string CrownsAppConfig::projection() {
	std::string projection;
	if(!m_originalCHM.empty() && isfile(m_originalCHM)){
		projection = chm().props().projection();
	} else if(projection.empty() && !m_smoothedCHM.empty() && isfile(m_smoothedCHM)) {
		projection = smoothed().props().projection();
	} else if(projection.empty() && !m_crownsRaster.empty() && isfile(m_crownsRaster)) {
		projection = crowns().props().projection();
	} else {
		g_warn("No file was found with a SRS to use for output.")
	}
	return projection;
}

void CrownsAppConfig::setBuildIndex(bool build) {
	if(m_locked) return;
	m_buildIndex = build;
	update(BuildIndex);
}

bool CrownsAppConfig::buildIndex() const {
	return m_buildIndex;
}

void CrownsAppConfig::setTableCacheSize(int size) {
	if(m_locked) return;
	m_tableCacheSize = size;
	update(TableCacheSize);
}

int CrownsAppConfig::tableCacheSize() const {
	return m_tableCacheSize;
}

void CrownsAppConfig::setRowCacheSize(int size) {
	if(m_locked) return;
	m_rowCacheSize = size;
	update(RowCacheSize);
}

int CrownsAppConfig::rowCacheSize() const {
	return m_rowCacheSize;
}

void CrownsAppConfig::setThreads(int threads) {
	if(m_locked) return;
	m_threads = threads;
	update(Threads);
}

int CrownsAppConfig::threads() const {
	return m_threads;
}

void CrownsAppConfig::setDoSmoothing(bool smoothing) {
	if(m_locked) return;
	m_doSmoothing = smoothing;
	update(DoSmoothing);
}

bool CrownsAppConfig::doSmoothing() const {
	return m_doSmoothing;
}

void CrownsAppConfig::setSmoothWindowSize(int size) {
	if(m_locked) return;
	m_smoothWindowSize = size;
	update(SmoothWindowSize);
}

int CrownsAppConfig::smoothWindowSize() const {
	return m_smoothWindowSize;
}

void CrownsAppConfig::setSmoothSigma(float sigma) {
	if(m_locked) return;
	m_smoothSigma = sigma;
	update(SmoothSigma);
}

float CrownsAppConfig::smoothSigma() const {
	return m_smoothSigma;
}

void CrownsAppConfig::setOriginalCHM(const std::string& filename, bool formatAll) {
	if(m_locked) return;
	m_originalCHM = filename;
	if(formatAll) {
		formatFilenames(filename,
				m_smoothedCHM, m_smoothedCHMDriver,
				m_topsDatabase, m_topsDatabaseDriver,
				m_crownsRaster, m_crownsRasterDriver,
				m_crownsDatabase, m_crownsDatabaseDriver,
				m_settings);
		update(OriginalCHM|SmoothedCHM|CrownsDatabase|CrownsRaster|CrownsDatabase|SettingsFile);
	} else {
		update(OriginalCHM);
	}
}

const std::string& CrownsAppConfig::originalCHM() const {
	return m_originalCHM;
}

void CrownsAppConfig::setSmoothedCHM(const std::string& filename) {
	if(m_locked) return;
	m_smoothedCHM = filename;
	update(SmoothedCHM);
}

const std::string& CrownsAppConfig::smoothedCHM() const {
	return m_smoothedCHM;
}

void CrownsAppConfig::setOriginalCHMBand(int band) {
	if(m_locked) return;
	m_originalCHMBand = band;
	update(OriginalCHMBand);
}

int CrownsAppConfig::originalCHMBand() const {
	return m_originalCHMBand;
}

void CrownsAppConfig::setSmoothedCHMDriver(const std::string& driver) {
	if(m_locked) return;
	m_smoothedCHMDriver = driver;
	if(!m_smoothedCHM.empty()) {
		replaceExt(m_smoothedCHM, rastExt(driver));
		update(SmoothedCHMDriver|SmoothedCHM);
	} else{
		update(SmoothedCHMDriver);
	}
}

const std::string& CrownsAppConfig::smoothedCHMDriver() const {
	return m_smoothedCHMDriver;
}

void CrownsAppConfig::setDoTops(bool tops) {
	if(m_locked) return;
	m_doTops = tops;
	update(DoTops);
}

bool CrownsAppConfig::doTops() const {
	return m_doTops;
}

void CrownsAppConfig::setTopsThresholds(const std::vector<TopThreshold>& thresholds) {
	if(m_locked) return;
	m_topsThresholds.assign(thresholds.begin(), thresholds.end());
	update(TopsThresholds);
}

const std::vector<TopThreshold>& CrownsAppConfig::topsThresholds() const {
	return m_topsThresholds;
}

void CrownsAppConfig::setTopsDatabase(const std::string& filename) {
	if(m_locked) return;
	m_topsDatabase = filename;
	update(TopsDatabase);
}

const std::string& CrownsAppConfig::topsDatabase() const {
	return m_topsDatabase;
}

void CrownsAppConfig::setTopsDatabaseDriver(const std::string& driver, bool forced) {
	if(m_locked) return;
	uint64_t msg = TopsDatabaseDriver;
	if(forced) {
		msg |= TopsDBFormatChanged;
		m_topsDatabaseDriver = "SQLite";
	} else {
		m_topsDatabaseDriver = driver;
	}
	if(!m_topsDatabase.empty())
		msg |= CrownsDatabase;
	if(forced || !m_topsDatabase.empty())
		replaceExt(m_topsDatabase, dbExt(driver));
	if(msg)
		update(msg);
}

const std::string& CrownsAppConfig::topsDatabaseDriver() const {
	return m_topsDatabaseDriver;
}

void CrownsAppConfig::setTopsMaxNulls(float nulls) {
	if(m_locked) return;
	m_topsMaxNulls = nulls;
	update(TopsMaxNulls);
}

float CrownsAppConfig::topsMaxNulls() const {
	return m_topsMaxNulls;
}

void CrownsAppConfig::setDoCrowns(bool crowns) {
	if(m_locked) return;
	m_doCrowns = crowns;
	update(DoCrowns);
}

bool CrownsAppConfig::doCrowns() const {
	return m_doCrowns;
}

void CrownsAppConfig::setCrownsThresholds(const std::vector<CrownThreshold>& thresholds) {
	if(m_locked) return;
	m_crownsThresholds.assign(thresholds.begin(), thresholds.end());
	update(CrownsThresholds);
}

const std::vector<CrownThreshold> CrownsAppConfig::crownsThresholds() const {
	return m_crownsThresholds;
}

void CrownsAppConfig::setCrownsUpdateHeights(bool updateHeights) {
	if(m_locked) return;
	m_crownsUpdateHeights = updateHeights;
	update(CrownsUpdateHeights);
}

bool CrownsAppConfig::crownsUpdateHeights() const {
	return m_crownsUpdateHeights;
}

void CrownsAppConfig::setCrownsRaster(const std::string& filename) {
	if(m_locked) return;
	m_crownsRaster = filename;
	update(CrownsRaster);
}

const std::string& CrownsAppConfig::crownsRaster() const {
	return m_crownsRaster;
}

void CrownsAppConfig::setCrownsDoDatabase(bool doDatabase) {
	if(m_locked) return;
	m_crownsDoDatabase = doDatabase;
	update(CrownsDoDatabase);
}

bool CrownsAppConfig::crownsDoDatabase() const {
	return m_crownsDoDatabase;
}

void CrownsAppConfig::setCrownsDatabase(const std::string& filename) {
	if(m_locked) return;
	m_crownsDatabase = filename;
	update(CrownsDatabase);
}

const std::string& CrownsAppConfig::crownsDatabase() const {
	return m_crownsDatabase;
}

void CrownsAppConfig::setCrownsRasterDriver(const std::string& driver) {
	if(m_locked) return;
	m_crownsRasterDriver = driver;
	if(!m_crownsRaster.empty()) {
		replaceExt(m_crownsRaster, rastExt(driver));
		update(CrownsRasterDriver|CrownsRaster);
	} else {
		update(CrownsRasterDriver);
	}
}

const std::string& CrownsAppConfig::crownsRasterDriver() const {
	return m_crownsRasterDriver;
}

void CrownsAppConfig::setCrownsDatabaseDriver(const std::string& driver, bool forced) {
	if(m_locked) return;
	uint64_t msg = CrownsDatabaseDriver;
	if(forced) {
		m_crownsDatabaseDriver = "SQLite";
		msg |= CrownsDBFormatChanged;
	} else {
		m_crownsDatabaseDriver = driver;
	}
	if(!m_crownsDatabase.empty())
		msg |= CrownsDatabase;
	if(forced || !m_crownsDatabase.empty())
		replaceExt(m_crownsDatabase, dbExt(driver));
	if(msg)
		update(msg);
}

const std::string& CrownsAppConfig::crownsDatabaseDriver() const {
	return m_crownsDatabaseDriver;
}

void CrownsAppConfig::setCrownsRemoveHoles(bool removeHoles) {
	if(m_locked) return;
	m_crownsRemoveHoles = removeHoles;
	update(CrownsRemoveHoles);
}

bool CrownsAppConfig::crownsRemoveHoles() const {
	return m_crownsRemoveHoles;
}

void CrownsAppConfig::setCrownsRemoveDangles(bool removeDangles) {
	if(m_locked) return;
	m_crownsRemoveDangles = removeDangles;
	update(CrownsRemoveDangles);
}

bool CrownsAppConfig::crownsRemoveDangles() const {
	return m_crownsRemoveDangles;
}

bool CrownsAppConfig::crownsKeepSmoothed() const {
	return m_crownsKeepSmoothed;
}

void CrownsAppConfig::setCrownsKeepSmoothed(bool keep) {
	if(m_locked) return;
	m_crownsKeepSmoothed = keep;
	update(CrownsKeepSmoothed);
}

void CrownsAppConfig::checkSmoothing() const {
	if (!doSmoothing())
		g_argerr("Not configured to perform smoothing.");
	if (originalCHM().empty())
		g_argerr("Smoothing: CHM filename must not be empty.");
	if (originalCHMBand() < 1)
		g_argerr("Smoothing: Band must be 1 or larger.");
	if (smoothedCHM().empty())
		g_argerr("Smoothing: Output filename must not be empty.");
	if (smoothedCHMDriver().empty())
		g_argerr("Smoothing: Output driver must not be empty.");
	if (smoothSigma() <= 0 || smoothSigma() > 100)
		g_argerr("Smoothing: Std. deviation must be 0 < n <= 100. " << smoothSigma() << " given.");
	if (smoothWindowSize() % 2 == 0 || smoothWindowSize() < 3)
		g_argerr("Smoothing: The window must be odd and >=3.");
}

void CrownsAppConfig::checkTops() const {
	if (!doTops())
		g_argerr("Not configured to find tops.");
	if (smoothedCHM().empty())
		g_argerr("Tops: Smoothed CHM filename must not be empty.");
	if (topsDatabase().empty())
		g_argerr("Tops: Database filename must not be empty.");
	if (topsDatabaseDriver().empty())
		g_argerr("Tops: Database driver must not be empty.");
	if (topsThresholds().empty()) {
		g_argerr("Tops: At least one threshold must be configured.");
	} else {
		float lastHeight;
		int lastWindow = 0;
		for(const TopThreshold& t : topsThresholds()) {
			if(t.threshold < 0.0)
				g_argerr("Threshold heights below zero are not allowed.");
			if(t.window % 2 == 0 || t.window < 3)
				g_argerr("Window size must be odd and >= 3.");
			if(lastWindow) {
				if(lastWindow >= t.window)
					g_argerr("Each window must be larger than the previous one.");
				if(lastHeight >= t.threshold)
					g_argerr("Each height must be larger than the previous one.");
			}
			lastWindow = t.window;
			lastHeight = t.threshold;
		}
	}
}

void CrownsAppConfig::checkCrowns() const {
	if (!doCrowns())
		g_argerr("Not configured to find crowns.");
	if (!doTops())
		g_argerr("Cannot run crowns without also running tops.");
	if (crownsRaster().empty())
		g_argerr("Crowns: Output raster filename must not be empty.");
	if (crownsRasterDriver().empty())
		g_argerr("Crowns: Output raster driver must not be empty.");
	if(crownsDoDatabase() && !crownsDatabase().empty() && crownsDatabaseDriver().empty())
		g_argerr("Crowns: If database file is given, driver must also be given.");
	if (topsDatabase().empty())
		g_argerr("Crowns: Crowns database filename must not be empty.");
	if (smoothedCHM().empty())
		g_argerr("Crowns: Smoothed CHM filename must not be empty.");
	if(crownsThresholds().empty()) {
		g_argerr("Crowns: At least one threshold must be configured.");
	} else {
		float lastHeight;
		int lastWindow = 0;
		for(const CrownThreshold& c : crownsThresholds()) {
			if(c.threshold < 0.0)
				g_argerr("Crown Thresholds: Threshold heights below zero are not allowed.");
			if(c.fraction <= 0.0)
				g_argerr("Crown Thresholds: Fraction height must be > 0.");
			if(c.radius <= 0.0)
				g_argerr("Crown Thresholds: Radius must be > 0.");
			if(lastWindow) {
				if(lastWindow >= c.radius)
					g_argerr("Crown Thresholds: Each radius must be larger than the previous one.");
				if(lastHeight >= c.threshold)
					g_argerr("Crown Thresholds: Each height must be larger than the previous one.");
			}
			lastWindow = c.radius;
			lastHeight = c.threshold;
		}
	}
}

void CrownsAppConfig::checkMerge() const {
	g_runerr("Not implemented.");
}

void CrownsAppConfig::check() const {
	if (doSmoothing())
		checkSmoothing();
	if (doTops())
		checkTops();
	if (doCrowns())
		checkCrowns();
}

bool CrownsAppConfig::canRun() const {
	try {
		check();
		return true;
	} catch (const std::exception& ex) {
		g_warn("Exception in canRun: " << ex.what());
	}
	return false;
}

std::string CrownsAppConfig::topsThresholdsList() const {
	std::vector<std::string> p(topsThresholds().size() * 2);
	char buf[256];
	int i = 0;
	for(const TopThreshold& t : topsThresholds()) {
		std::sprintf(buf, "%.2f", t.threshold);
		p[i++] = std::string(buf);
		std::sprintf(buf, "%u", t.window);
		p[i++] = std::string(buf);
	}
	return join(p.begin(), p.end(), ",");
}

void CrownsAppConfig::parseTopsThresholds(const std::string& str) {
	std::stringstream ss(str);
	std::string item1, item2;
	std::vector<TopThreshold> config;
	while(std::getline(ss, item1, ',') && std::getline(ss, item2, ',')) {
		float height = atof(item1.c_str());
		int window = atoi(item2.c_str());
		config.emplace_back(height, window);
	}
	setTopsThresholds(config);
}

std::string CrownsAppConfig::crownsThresholdsList() const {
	std::vector<std::string> p(crownsThresholds().size() * 3);
	char buf[256];
	int i = 0;
	for(const CrownThreshold& c : crownsThresholds()) {
		std::sprintf(buf, "%.2f", c.threshold);
		p[i++] = std::string(buf);
		std::sprintf(buf, "%.2f", c.fraction);
		p[i++] = std::string(buf);
		std::sprintf(buf, "%.2f", c.radius);
		p[i++] = std::string(buf);
	}
	return join(p.begin(), p.end(), ",");
}

void CrownsAppConfig::parseCrownsThresholds(const std::string& str) {
	std::stringstream ss(str);
	std::string item1, item2, item3;
	std::vector<CrownThreshold> config;
	while(std::getline(ss, item1, ',') && std::getline(ss, item2, ',') && std::getline(ss, item3, ',')) {
		float height = atof(item1.c_str());
		float frac = atof(item2.c_str());
		float radius = atof(item3.c_str());
		config.emplace_back(height, frac, radius);
	}
	setCrownsThresholds(config);
}

CrownsAppConfig::~CrownsAppConfig() {
	destroy();
}

// Crowns implementation

Crowns::Crowns() {}

CrownsAppConfig& Crowns::config() {
	return m_config;
}

void Crowns::smooth() {

	m_config.checkSmoothing();
	Monitor::get().status(0.01f, "Smoothing...");

	Band<float>& chm = m_config.chm();

	if(Monitor::get().canceled())
		return;

	Band<float>& smoothed = m_config.smoothed();
	smoothed.fill(smoothed.props().nodata());

	if(Monitor::get().canceled())
		return;

	chm.smooth(smoothed, m_config.smoothSigma(), m_config.smoothWindowSize());

	m_config.flush();
	Monitor::get().status(1.0f, "Smoothing: Done.");
}

void Crowns::run() {
	if(m_config.doTops())
		runTops();
	if(m_config.doCrowns())
		runCrowns();
}

void Crowns::runTops() {

	m_config.checkTops();
	Monitor::get().status(0.01f, "Crowns: Preparing...");

	// Initialize input rasters.
	Band<float>& smoothed = m_config.smoothed();

	if(Monitor::get().canceled())
		return;

	// Get the grid containing top windows and zero it.
	Band<uint32_t>& topsWindowGrid = m_config.topsWindow();
	topsWindowGrid.fill(0);

	if(Monitor::get().canceled())
		return;

	// Get the grid containing tree top IDs and zero it.
	Band<uint32_t>& topsIDGrid = m_config.topsID();
	topsIDGrid.fill(0);

	if(Monitor::get().canceled())
		return;

	Monitor::get().status(0.02f, "Crowns: Finding tops...");

	if (Monitor::get().canceled())
		return;

	// Find the tops in the smoothed raster.
	stage1(m_config);

	Monitor::get().status(0.31f, "Crowns: Finding parent tops...");

	if (Monitor::get().canceled())
		return;

	// Get the grid containing tree top parent tops and zero it.
	Band<uint32_t>& topsParentGrid = m_config.topsParent();
	topsParentGrid.fill(0);

	if(Monitor::get().canceled())
		return;

	// Find parent tops.
	stage2(m_config);

	if(Monitor::get().canceled())
		return;

	// Finally, scrape up all the tops and put them in the DB.
	Monitor::get().status(0.61f, "Crowns: Saving tops...");

	if (Monitor::get().canceled())
		return;

	// Must initialize the tree with the raster bounds.
	mqtree<Top>& topsTree = m_config.topsTree();
	{
		const GridProps& props = smoothed.props();
		const Bounds<double>& bounds = props.bounds();
		topsTree.init(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy());
		topsTree.reset();
	}

	// Create treetop objects.
	stage3(m_config);

	if (Monitor::get().canceled())
		return;

	// Save tops, but only if crowns aren't going to be performed also.
	// At this stage only the smoothed heights/positions are known.
	CrownDB& cdb = m_config.crownDB();
	Top top;
	// Add the tree tops to the internal database.
	topsTree.reset();
	while(topsTree.next(top) && !Monitor::get().canceled())
		cdb.insert(top);

	if(Monitor::get().canceled())
		return;

	// If not doing crowns, save the tops now.
	if(!m_config.doCrowns()) {

		// Check if the save is likely to succeed. Change the driver if not.
		if(!cdb.checkSaveTops(m_config.topsDatabaseDriver()))
			m_config.setCrownsDatabaseDriver("SQLite", true);

		try {
			cdb.saveTops(
				m_config.topsDatabase(),
				m_config.topsDatabaseDriver(),
				"tops",
				m_config.projection()
			);
		} catch(const CrownDBSaveException& ex) {
			// If it fails on save, change the driver and try again. If it fails once more,
			// will rethrow and exit.
			m_config.setCrownsDatabaseDriver("SQLite", true);
			cdb.saveTops(
				m_config.topsDatabase(),
				m_config.topsDatabaseDriver(),
				"tops",
				m_config.projection()
			);
		} catch(const std::exception& ex) {
			throw ex;
		}

	}
	topsTree.reset();

	m_config.flush();

	Monitor::get().status(1.0f, "Crowns: Done.");

}

void Crowns::runCrowns() {

	if(!m_config.doTops())
		g_runerr("Tops are required for running crowns.")

	m_config.checkCrowns();

	Monitor::get().status(0.01f, "Crowns: Preparing...");

	if(Monitor::get().canceled())
		return;

	// Delineate crowns on the smoothed raster.
	delineateCrowns(m_config);

	if(Monitor::get().canceled())
		return;

	// Save the crowns file.
	m_config.crowns().flush();

	if(Monitor::get().canceled())
		return;

	// Update the heights of tops from the original CHM.
	updateOriginalCHMHeights(m_config);

	if(Monitor::get().canceled())
		return;

	// Create the tops database.
	Monitor::get().status(0.0f, "Crowns: Saving tops...");

	// Check if the save is likely to succeed. Change the driver if not.
	if(!m_config.crownDB().checkSaveTops(m_config.topsDatabaseDriver()))
		m_config.setCrownsDatabaseDriver("SQLite", true);

	try{
		m_config.crownDB().saveTops(
			m_config.topsDatabase(),
			m_config.topsDatabaseDriver(),
			"tops",
			m_config.projection()
		);
	} catch(const CrownDBSaveException& ex) {
		// If it fails on save, change the driver and try again. If it fails once more,
		// will rethrow and exit.
		m_config.setCrownsDatabaseDriver("SQLite", true);
		m_config.crownDB().saveTops(
			m_config.topsDatabase(),
			m_config.topsDatabaseDriver(),
			"tops",
			m_config.projection()
		);
	} catch(const std::exception& ex) {
		throw ex;
	}

	if(Monitor::get().canceled())
		return;

	// Update the values in the crowns table from the tops table.
	if(m_config.crownsDoDatabase()) {
		polygonizeCrowns(m_config);
		Monitor::get().status(0.0f, "Crowns: Saving crowns...");

		if(Monitor::get().canceled())
			return;

		// Check if the save is likely to succeed. Change the driver if not.
		if(!m_config.crownDB().checkSaveCrowns(m_config.crownsDatabaseDriver()))
			m_config.setCrownsDatabaseDriver("SQLite", true);

		try {
			m_config.crownDB().saveCrowns(
				m_config.crownsDatabase(),
				m_config.crownsDatabaseDriver(),
				"crowns",
				m_config.projection()
			);
		} catch(const CrownDBSaveException& ex) {
			// If it fails on save, change the driver and try again. If it fails once more,
			// will rethrow and exit.
			m_config.setCrownsDatabaseDriver("SQLite", true);
			m_config.crownDB().saveCrowns(
				m_config.crownsDatabase(),
				m_config.crownsDatabaseDriver(),
				"crowns",
				m_config.projection()
			);
		} catch(const std::exception& ex) {
			throw ex;
		}
	}

	m_config.destroy();

	Monitor::get().status(1.0, "Crowns: Done.");

}

Crowns::~Crowns() {
}


