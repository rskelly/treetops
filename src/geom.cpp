#include <fcntl.h>

#include "treetops.hpp"
#include "geom.hpp"
#include "util.hpp"

using namespace tt::util;
using namespace tt::util::vec;
using namespace tt::vec;
using namespace tt::data;
using namespace tt::config;


static const int CDBBUFSIZE = sizeof(int) * 3 + sizeof(tt::data::Treetop);


CrownDB::CrownDB() :
		m_tfd(0), m_cfd(0),
		m_toffset(0), m_coffset(0),
		m_tidx(0), m_count(0) {

	m_gctx = GEOS_init_r();
	m_rdr = GEOSWKBReader_create_r(m_gctx);
	m_wtr = GEOSWKBWriter_create_r(m_gctx);

	if(!tmpfile("tt", "", m_tfile) || (m_tfd = open(m_tfile.c_str(), O_CREAT|O_RDWR, 0777)) <= 0)
		_runerr("Failed to open tops DB file: " << m_tfile);
    
	if(!tmpfile("tt", "", m_cfile) || (m_cfd = open(m_cfile.c_str(), O_CREAT|O_RDWR, 0777)) <= 0)
		_runerr("Failed to open crowns DB file:" << m_cfile);
}

GEOSContextHandle_t CrownDB::gctx() {
	return m_gctx;
}

bool CrownDB::next(Treetop& top) {
	static GEOSGeometry* crown;
	return next(top, crown, nullptr, false);
}

bool CrownDB::next(Treetop& top, std::vector<unsigned char>& crownwkb, bool loadCrown) {
	int id;
	int plen;
	int ppos;
	int pos = 0;
	std::vector<unsigned char> buf(CDBBUFSIZE);

	// Loop until a next one is found or the file runs out.
	while(true) {
		// Go to start of line at offset.
		if (-1 == lseek(m_tfd, m_tidx * CDBBUFSIZE, SEEK_SET)) {
			_warn("Failed to seek in tops DB file: " << strerror(errno));
			return false;
		}
		// Read the line into the buffer.
		int r;
		if ((r = read(m_tfd, buf.data(), CDBBUFSIZE)) < (ssize_t)CDBBUFSIZE) {
			// The file may be finished; not an error.
			if(r != 0)
				_warn("Failed to read tops DB file (this may not be an error): " << r << "; " << strerror(errno));
			return false;
		}
		// Read the values.
		std::memcpy(&id, buf.data(), sizeof(int));			
		pos += sizeof(int);
		std::memcpy(&top, buf.data() + pos, sizeof(Treetop));		
		pos += sizeof(Treetop);
		std::memcpy(&plen, buf.data() + pos, sizeof(int));	
		pos += sizeof(int);
		// If there's a poly, the plen will be > 0. Read it. Otherwise crown is null.
		if(loadCrown && plen > 0) {
			std::memcpy(&ppos, buf.data() + pos, sizeof(int));
			crownwkb.resize(plen);

			// Go to the poly.
			if(-1 == lseek(m_cfd, ppos, SEEK_SET)) {
				_warn("Failed to seek in crowns DB file: " << strerror(errno));
				return false;
			}
			if((r = read(m_cfd, crownwkb.data(), plen)) < (ssize_t) plen) {
				_warn("Failed to read " << plen << " bytes from crowns DB file (this is an error): " << r << "; " << strerror(errno));
				return false;
			}

		}

		++m_tidx;
		return true;
	}
}

bool CrownDB::next(Treetop& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx, bool loadCrown) {
	static std::vector<unsigned char> cbuf;
	bool ret = next(top, cbuf, loadCrown);
	if(ret && loadCrown) {
		if(!gctx)
			gctx = m_gctx;
		if(!(crown = GEOSWKBReader_read_r(gctx, m_rdr, cbuf.data(), cbuf.size())))
			return false;
	}
	return ret;
}

void CrownDB::reset() {
	m_tidx = 0;
}

bool CrownDB::insert(const Treetop& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx) {

	// If it exists, update it.
	if(update(top, crown, gctx))
		return true;

	int pos = 0;
	int plen = 0;
	int ppos = 0;
	std::vector<unsigned char> buf(CDBBUFSIZE);

	// Write the ID, top.
	std::memcpy(buf.data(), &top.id, sizeof(int));		
	pos += sizeof(int);
	std::memcpy(buf.data() + pos, &top, sizeof(Treetop));		
	pos += sizeof(Treetop);

	// Write the geometry.
	if(crown) {
		if(!gctx)
			gctx = m_gctx;
		if(!writeCrown(plen, ppos, crown, gctx))
			return false;
	}

	// Write back the crown length and position.
	std::memcpy(buf.data() + pos, &plen, sizeof(int));	
	pos += sizeof(int);
	std::memcpy(buf.data() + pos, &ppos, sizeof(int));

	// Write the Treetop record.
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

bool CrownDB::find(int id, Treetop& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx, bool withCrown) {

	if(id <= 0)
		return false;

	int offset = -1;
	{
		std::lock_guard<std::mutex> lk(m_imtx);
		if(m_idx.find(id) != m_idx.end())
			offset = m_idx.at(id);
	}

	if(offset == (int) -1)
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

	int id_;
	int plen;
	int ppos;
	int pos = 0;

	// Read the values.
	std::memcpy(&id_, buf.data(), sizeof(int));			
	pos += sizeof(int);
	std::memcpy(&top, buf.data() + pos, sizeof(Treetop));			
	pos += sizeof(Treetop);
	std::memcpy(&plen, buf.data() + pos, sizeof(int));		
	pos += sizeof(int);

	// If there's a poly, the plen will be > 0. Read it. Otherwise crown is null.
	if(withCrown && plen > 0) {
		if(!gctx)
			gctx = m_gctx;
		std::memcpy(&ppos, buf.data() + pos, sizeof(int));
		std::vector<unsigned char> pbuf(plen);

		std::lock_guard<std::mutex> lk(m_cmtx);
		// Go to the poly.
		if(-1 == lseek(m_cfd, ppos, SEEK_SET))
			return false;
		if(read(m_cfd, pbuf.data(), plen) < (ssize_t) plen)
			return false;
		if(!(crown = GEOSWKBReader_read_r(gctx, m_rdr, pbuf.data(), plen)))
			return false;
	} else {
		crown = nullptr;
	}

	return true;
}

int CrownDB::count() {
	return m_count;
}

bool CrownDB::writeCrown(int& plen, int& ppos, const GEOSGeometry* crown, GEOSContextHandle_t gctx) {
	bool ret = false;
	if(crown) {
		if(!gctx)
			gctx = m_gctx;

		// Write the geometry to a char buffer.
		size_t plen0;
		unsigned char* poly = GEOSWKBWriter_write_r(gctx, m_wtr, crown, &plen0);

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

		GEOSFree_r(gctx, poly);
	}
	return ret;
}

bool CrownDB::update(const Treetop& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx) {

	int offset = -1;

	if(top.id > 0) {
		std::lock_guard<std::mutex> lk(m_imtx);
		if(m_idx.find(top.id) != m_idx.end())
			offset = m_idx.at(top.id);
	}

	if(offset == (int) -1)
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
	int plen;
	int ppos;
	int pos = sizeof(int) + sizeof(Treetop); // Skip the ID and Treetop.
	std::memcpy(&plen, buf.data() + pos, sizeof(int));		
	pos += sizeof(int);
	std::memcpy(&ppos, buf.data() + pos, sizeof(int));

	// Write the Treetop.
	pos = sizeof(int);
	std::memcpy(buf.data() + pos, &top, sizeof(Treetop));			
	pos += sizeof(Treetop);

	if(crown) {
		if(!gctx)
			gctx = m_gctx;
		// Write the crown, update length and position in the crown file.
		if(!writeCrown(plen, ppos, crown, gctx))
			return false;
	} else {
		// If there's no crown, zero out the position/length.
		ppos = 0;
		plen = 0;
	}

	// Write the original pos/size back.
	std::memcpy(buf.data() + pos, &plen, sizeof(int));		
	pos += sizeof(int);
	std::memcpy(buf.data() + pos, &ppos, sizeof(int));

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
		_runerr("Failed to start transaction.");

	Treetop top;
	GEOSGeometry* geom;
	int count = 0;

	reset();

	while(next(top, geom, gctx, false)) {
		double x, y, z, th, gz;
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
	    feat.SetField(1, (GIntBig) top.parentId);
	    feat.SetField(2, gz);
	    feat.SetField(3, top.sz);
	    feat.SetField(4, th);
	    feat.SetField(5, z);
		OGRPoint geom(x, y, z);
		feat.SetGeometry(&geom);
	    if(OGRERR_NONE != layer->CreateFeature(&feat))
	    	_runerr("Failed to add feature to " << filename);

	}

	// Commit and release the layer -- GDAL will take care of it. But close the dataset so that can happen.
	if(OGRERR_NONE != layer->CommitTransaction())
		_runerr("Failed to commit transaction.");
	layer->Dereference();
	GDALClose(ds);

	OGRGeometry::freeGEOSContext(gctx);

}

bool CrownDB::checkSaveCrowns(const std::string& driver) {
	int gb = 2L * 1024L * 1024L * 1024L;
	// This is an empirically-derived estimate for the filesize.
	// Crowns are variable-length and this may still fail in the save step.
	if(lowercase(driver) == "esri shapefile" && 2 * crownsDBSize() > gb)
		return false;
	return true;
}

bool CrownDB::checkSaveTops(const std::string& driver) {
	int gb = 2L * 1024L * 1024L * 1024L;
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
		_runerr("Failed to start transaction.");

	Treetop top;
	GEOSGeometry* geom;
	int count = 0;

	reset();

	while(next(top, geom, gctx, true)) {
		if(geom == nullptr) {
			// _warn("Null geometry.");
			continue;
		}
		double z, th, gz;
		if(top.oz > 0) {
			z = top.oz;				// The "true" height from the unsmoothed raster.
			gz = top.groundZ;		// The height of ground, if a terrain model is given.
			th = z * 0.5 + gz;		// ?
		} else {
			z = top.sz;
			gz = 0;
			th = 0;
		}

		// Set the z-component on the geometry and all its children.
		GEOSGeometry* geom0 = setGeomZ(gctx, geom, th);
		GEOSGeom_destroy_r(gctx, geom);

		// Note: Fields added in same order as defined.
	    OGRFeature feat(layer->GetLayerDefn());
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) geom0);
	    feat.SetField(0, (GIntBig) top.id);
	    feat.SetField(1, (GIntBig) top.parentId);
	    feat.SetField(2, gz);
	    feat.SetField(3, top.sz);
	    feat.SetField(4, th);
	    feat.SetField(5, z);
		feat.SetGeometry(ogeom);
		OGRErr res = layer->CreateFeature(&feat);
		OGRGeometryFactory::destroyGeometry(ogeom);
		GEOSGeom_destroy_r(gctx, geom0);
	    if(res != OGRERR_NONE)
	    	_runerr("Failed to add feature to " << filename);

	}

	// Commit and release the layer -- GDAL will take care of it. But close the dataset so that can happen.
	if(OGRERR_NONE != layer->CommitTransaction())
		_runerr("Failed to commit transaction.");
	layer->Dereference();
	GDALClose(ds);

	OGRGeometry::freeGEOSContext(gctx);

}

int CrownDB::crownsDBSize() {
	return filesize(m_cfile);
}

int CrownDB::topsDBSize() {
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
