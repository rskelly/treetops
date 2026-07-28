/*
 * contrem_util.cpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */


#ifdef _WIN32
#include <io.h>
#include <winbase.h>
#include <processthreadsapi.h>
#include <minwinbase.h>
#include <sysinfoapi.h>
#else
#include <sys/time.h>
#include <glob.h>
#endif

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <sstream>
#include <regex>
#include <fstream>
#include <random>
#include <filesystem>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "util.hpp"


using namespace tt::util;

namespace fs = std::filesystem;

namespace {

	const char pathsep = std::filesystem::path::preferred_separator;

	const std::string defaultChars = "abcdefghijklmnaoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";

	//https://codereview.stackexchange.com/a/118957
	std::string randomString(size_t len = 15) {
		std::mt19937_64 gen { std::random_device()() };
	    std::uniform_int_distribution<size_t> dist { 0, defaultChars.length()-1 };
	    std::string ret;
	    std::generate_n(std::back_inserter(ret), len, [&] {
	    	return defaultChars[dist(gen)];
	    });
	    return ret;
	}

} // anon


FileType tt::util::getFileType(const std::string& filename) {
	std::string ext;
	{
		size_t p = filename.find('.');
		if(p < std::string::npos) {
			std::string ext0 = filename.substr(p);
			std::transform(ext0.begin(), ext0.end(), std::back_inserter(ext), ::tolower);
		}
	}
	if(ext == ".csv" || ext == ".txt") {
		return FileType::CSV;
	} else if(ext == ".roi") {
		return FileType::ROI;
	} else {
		GDALAllRegister();
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, 0, 0, 0));
		if(ds) {
			std::string drv(ds->GetDriverName());
			FileType type = FileType::Unknown;
			if(drv == "GTiff") {
				type = FileType::GTiff;
			} else if(drv == "ENVI") {
				type = FileType::ENVI;
			} else if(drv == "ESRI Shapefile") {
				type = FileType::SHP;
			} else if(drv == "SQLite") {
				type = FileType::SQLITE;
			}
			GDALClose(ds);
			return type;
		}
	}
	return FileType::Unknown;
}

std::string tt::util::fileTypeAsString(FileType type) {
	switch(type) {
	case FileType::GTiff: return "GTiff";
	case FileType::ENVI: return "ENVI";
	case FileType::ROI: return "ENVI ROI";
	case FileType::SHP: return "Shapefile";
	case FileType::CSV: return "CSV";
	default: return "";
	}
}

FileType tt::util::fileTypeFromString(const std::string& type) {
	if(type == "GTiff") {
		return FileType::GTiff;
	} else if(type == "ENVI") {
		return FileType::ENVI;
	} else if(type == "ENVI ROI") {
		return FileType::ROI;
	} else if(type == "Shapefile" || type == "ESRI Shapefile") {
		return FileType::SHP;
	} else if(type == "CSV") {
		return FileType::CSV;
	} else {
		return FileType::Unknown;
	}
}

std::string tt::util::normMethodAsString(NormMethod method) {
	switch(method) {
	case NormMethod::ConvexHull:
		return "Convex Hull";
	case NormMethod::ConvexHullLongestSeg:
		return "Convex Hull, Longest Segment";
	case NormMethod::Line:
		return "Line";
	case NormMethod::Unknown:
	default:
		return "Unknown";
	}
}

NormMethod tt::util::normMethodFromString(const std::string& method) {
	if(method == normMethodAsString(NormMethod::ConvexHull)) {
		return NormMethod::ConvexHull;
	} else if(method == normMethodAsString(NormMethod::ConvexHullLongestSeg)) {
		return NormMethod::ConvexHullLongestSeg;
	} else if(method == normMethodAsString(NormMethod::Line)) {
		return NormMethod::Line;
	} else {
		return NormMethod::Unknown;
	}
}

bool tt::util::isnonzero(const double& v) {
	return v != 0;
}

bool tt::util::checkValidInputFiles(const std::vector<std::string>& files) {
	if(files.empty()) {
		_warn("Input file list is empty.");
		return false;
	}
	int invalid = 0;
	for(const std::string& f : files) {
		if(f.empty()) {
			_warn("File name is empty.");
			++invalid;
		}
		if(!isfile(f)) {
			_warn(f << " is invalid or doesn't exist.");
			++invalid;
		}
	}
	return invalid == 0;
}

bool tt::util::safeToWrite(const std::string& path, bool force) {
	if(path.empty() || isdir(path))
		return false;
	if(isfile(path))
		return force;
	return true;
}

bool tt::util::exists(const std::string& path) {
	return isdir(path) || isfile(path);
}

bool tt::util::isdir(const std::string& path) {
	return !path.empty() && fs::is_directory(path);
}

bool tt::util::isfile(const std::string& path) {
	return !path.empty() && fs::is_regular_file(path);
}

bool tt::util::rem(const std::string& dir) {
	try {
		fs::path p(dir);
		fs::remove_all(p);
	} catch (const std::exception& ex) {
		_warn(ex.what());
		return false;
	}
	return true;
}

std::vector<std::string> tt::util::glob(const std::string& path) {
	std::vector<std::string> files;
#ifdef _WIN32
	std::string p = parent(path);
	WIN32_FIND_DATA data;
	HANDLE fh = FindFirstFile(path.c_str(), &data);
	if(fh != INVALID_HANDLE_VALUE) {
		files.push_back(join(p, data.cFileName));
		while (FindNextFileA(fh, &data)) {
			if(fh != INVALID_HANDLE_VALUE)
				files.push_back(join(p, data.cFileName));
		}
		FindClose(fh);
	}
#else
	//https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(path.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(return_value != 0) {
        globfree(&glob_result);
        _runerr("glob() failed with return_value " << return_value);
    }

    // collect all the filenames into a std::list<std::string>
    for(size_t i = 0; i < glob_result.gl_pathc; ++i)
        files.push_back(glob_result.gl_pathv[i]);

    // cleanup
    globfree(&glob_result);
#endif

    return files;
}

std::string tt::util::parent(const std::string& path) {
	std::string _p = path;

	while(_p.size() > 0 && _p.back() == pathsep)
		_p = _p.substr(0, _p.size() - 1);

	if(_p.empty())
		return _p;

	size_t a = _p.find_last_of(pathsep);
	if(a == std::string::npos)
		return "";

	return _p.substr(0, a);

}

size_t tt::util::filesize(const std::string& f) {
	fs::path p(f);
	return fs::file_size(p);
}

bool tt::util::rename(const std::string& from, const std::string& to) {
	if (isdir(to))
		_runerr(to << " is a directory.")
	return fs::copy_file(from, to);
}

std::string tt::util::join(const std::string& a, const std::string& b) {
	if(b.empty()) {
		return a.empty() ? "" : a;
	} else if(a.empty()) {
		return b.empty() ? "" : b;
	}
	std::string _a, _b;
	for(size_t i = a.size() - 1; i < std::string::npos; --i) {
		if(a[i] != pathsep) {
			_a = a.substr(0, i + 1);
			break;
		}
	}
	for(size_t i = 0; i < b.size(); ++i) {
		if(b[i] != pathsep) {
			_b = b.substr(i);
			break;
		}
	}
	return _a + pathsep + _b;
}

std::string tt::util::basename(const std::string& path) {
	std::string _p = path;

	while(_p.size() > 0 && _p.back() == pathsep)
		_p = _p.substr(0, _p.size() - 1);

	if(_p.empty())
		return _p;

	size_t a = _p.find_last_of(pathsep);
	size_t b = _p.find_last_of('.');
	return _p.substr(a + 1, b - a - 1);
}

std::string tt::util::extension(const std::string& path) {
	size_t pos = path.find_last_of('.');
	if(pos < std::string::npos)
		return path.substr(pos, std::string::npos);
	return path;
}

std::string tt::util::gettmpdir() {
	fs::path p = fs::temp_directory_path();
	return p.string();
}

int tt::util::pid() {
#ifdef _WIN32
	return GetCurrentProcessId();
#else
	return getpid();
#endif
}

bool tt::util::tmpdir(const std::string& prefix, const std::string& dir, std::string& result) {
	// Assemble the target directory, check and attempt to create if needed.
	std::string tdir = join(gettmpdir(), dir);
	if(!isdir(tdir)) {
		if(!makedir(tdir))
			_runerr("Failed to make target dir: " << tdir);
	}
	int tries = 16;
	while(--tries) {
		std::stringstream ss;
		ss << join(tdir, prefix) << '_' << pid() << '_' << randomString();
		std::string path = ss.str();
		if(makedir(path)) {
			result = path;
			return true;
		}
	}
	_runerr("Failed to make temporary directory.");
	return false;
}

bool tt::util::tmpfile(const std::string& prefix, const std::string& dir, std::string& result) {
	// Assemble the target directory, check and attempt to create if needed.
	std::string tdir = join(gettmpdir(), dir);
	if(!isdir(tdir)) {
		if(!makedir(tdir))
			_runerr("Failed to make target directory: " << tdir);
	}
	int tries = 16;
	while(--tries) {
		// Assemble the file path.
		std::stringstream ss;
		ss << join(tdir, prefix) << '_' << pid() << '_' << randomString();
		std::string path = ss.str();
		if(!isfile(path)) {
			result = path;
			return true;
		}
	}
	_runerr("Failed to create non-extant filename.");
	return false;
}

bool tt::util::makedir(const std::string& path) {
	bool res = fs::create_directories(path);
	if (res)
		fs::permissions(path, fs::perms::owner_all | fs::perms::group_all);
	return res;
}

std::string tt::util::sanitize(const std::string& str) {
	std::regex repl("([^0-9A-Za-z]+)");
	std::stringstream ss;
	std::regex_replace(std::ostreambuf_iterator<char>(ss), str.begin(), str.end(), repl, "_");
	return ss.str();
}

std::string tt::util::lowercase(const std::string& str) {
	std::string out;
	std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
	return out;
}

int tt::util::gdalTypeSize(GDALDataType type) {
	switch(type) {
	case GDT_Float32:
	case GDT_Int32:
	case GDT_UInt32: 	return 4;
	case GDT_Int16:
	case GDT_UInt16: 	return 2;
	case GDT_Float64: 	return 8;
	case GDT_Byte: 		return 1;
	default:
		throw std::runtime_error("Unknown GDAL data type: " + std::to_string((int) type));
	}
}

std::string tt::util::projectionFromSRID(int srid) {
	std::string out;
	OGRSpatialReference osr;
	if((OGRERR_NONE == osr.importFromEPSG(srid))) {
		char* wkt;
		osr.exportToWkt(&wkt);
		out = wkt;
		CPLFree(wkt);
	}
	return out;
}

TmpFile::TmpFile(size_t size) :
	fd(0), size(0) {
	if(!tmpfile("geo_util", "", filename))
		_runerr("Failed to create temp file.");
	fd = ::open(filename.c_str(), O_CREAT|O_RDWR, 0777);
	if(fd <= 0)
		_runerr("Failed to open temp file: " << strerror(errno));
	resize(size);
}

void TmpFile::resize(size_t newSize) {
	if(newSize > size) {
		if(fd <= 0)
			throw std::runtime_error(std::string("File is not open: ") + strerror(errno));

		if((size_t) lseek(fd, newSize - 1, SEEK_SET) != (newSize - 1))
			throw std::runtime_error(std::string("Failed to create temporary file for mapping.") + strerror(errno));

		if(write(fd, "", 1) < 1)
			throw std::runtime_error(std::string("Failed to create temporary file for mapping.") + strerror(errno));

		size = newSize;
	}
}

void TmpFile::close() {
	::close(fd);
}

TmpFile::~TmpFile() {
	::close(fd);
	rem(filename);
}

constexpr uint32_t MAX_UINT32 = 32768;

uint64_t tt::util::morton(uint32_t x, uint32_t y) {
	if(x >= MAX_UINT32)
		_runerr("x coordinate is too large for bit shuffling: " << x)
	if(y >= MAX_UINT32)
		_runerr("y coordinate is too large for bit shuffling: " << y)
	uint32_t xa = bitsplit2(x);
	uint32_t ya = bitsplit2(y);
	return xa | (ya << 1);
}

double tt::util::random(double min, double max) {
	return min + ((double) rand() / RAND_MAX) * (max - min);
}

uint64_t tt::util::microtime() {
#ifdef _WIN32
	FILETIME ft;
	GetSystemTimePreciseAsFileTime(&ft);
	return ((uint64_t) ft.dwHighDateTime << 32) | ft.dwLowDateTime;
#else
	struct timeval t;
	gettimeofday(&t, NULL);
	return (uint64_t) t.tv_sec * 1000000 + t.tv_usec;
#endif
}


using namespace tt::util::csv;

int CSVValue::asInt() const {
	if(t == Int) {
		return i;
	} else {
		return (int) d;
	}
}

double CSVValue::asDouble() const {
	if(t == Double) {
		return d;
	} else {
		return (double) i;
	}
}

const std::string& CSVValue::asString() const {
	return s;
}


bool CSV::isdouble(const std::string& s) {
	if(s == "inf" || s == "-inf" || s == "NaN")
		return true;
	for(size_t i = 0; i < s.size(); ++i) {
		if(!std::isdigit(s[i], std::locale()) && s[i] != '.' && s[i] != '+' && s[i] != '-' && s[i] != 'e')
			return false;
	}
	return true;
}

bool CSV::isint(const std::string& s) {
	for(size_t i = 0; i < s.size(); ++i) {
		if(!std::isdigit(s[i], std::locale()))
			return false;
	}
	return true;
}

CSV::CSV(const std::string& file, bool header) {
	if(!file.empty())
		load(file, header);
}

void CSV::load(const std::string& file, bool header) {
	std::ifstream in(file);
	std::string line;
	std::string cell;
	std::vector<std::string> names;
	std::vector<CSVType> types;
	std::vector<std::vector<std::string>> values;
	int colCount = 0;
	bool doNames = true;
	while(std::getline(in, line)) {
		std::stringstream ss(line);
		if(header) {
			while(std::getline(ss, cell, ','))
				names.push_back(cell);
			colCount = names.size();
			values.resize(colCount);
			header = false;
			doNames = false;
		} else {
			int idx = 0;
			while(std::getline(ss, cell, ',')) {
				if(idx == colCount)
					break;
				if(doNames) {
					names.push_back("col_" + std::to_string(++colCount));
					values.resize(colCount);
				}
				values[idx++].push_back(cell);
			}
			doNames = false;
		}
	}
	types.resize(names.size());
	for(size_t i = 0; i < names.size(); ++i) {
		int t = 2;
		for(size_t j = 0; j < values[i].size(); ++j) {
			if(t == 2 && !isint(values[i][j])) {
				--t;
			} else if(t == 1 && !isdouble(values[i][j])) {
				--t;
				break;
			}
		}
		switch(t) {
		case 2: types[i] = Int; break;
		case 1: types[i] = Double; break;
		default: types[i] = String; break;
		}
	}
	m_values.resize(names.size());
	for(size_t i = 0; i < names.size(); ++i) {
		m_values[i].name = names[i];
		m_values[i].type = types[i];
		m_values[i].values.resize(values[i].size());
		for(size_t j = 0; j < values[i].size(); ++j) {
			m_values[i].values[j].t = types[i];
			switch(types[i]) {
			case Double:
				m_values[i].values[j].d = atof(values[i][j].c_str());
				break;
			case Int:
				m_values[i].values[j].i = atoi(values[i][j].c_str());
				break;
			default:
				m_values[i].values[j].s = values[i][j];
				break;
			}
		}
		values[i].clear();
	}
}

std::vector<std::string> CSV::columnNames() const {
	std::vector<std::string> names;
	for(const CSVColumn& c : m_values)
		names.push_back(c.name);
	return names;
}

std::vector<CSVValue> CSV::row(size_t idx) const {
	if(idx >= m_values[0].values.size())
		throw std::runtime_error("Index is too large.");
	std::vector<CSVValue> row;
	for(size_t i = 0; i < m_values.size(); ++i)
		row.push_back(m_values[i].values[idx]);
	return row;
}

std::vector<CSVValue> CSV::column(const std::string& name) const {
	for(size_t i = 0; i < m_values.size(); ++i) {
		if(m_values[i].name == name)
			return m_values[i].values;
	}
	throw std::runtime_error("No column named " + name);
}

CSVType CSV::columnType(const std::string& name) const {
	for(size_t i = 0; i < m_values.size(); ++i) {
		if(m_values[i].name == name)
			return m_values[i].type;
	}
	throw std::runtime_error("No column named " + name);
}

std::vector<CSVValue> CSV::column(size_t i) const {
	if(i < m_values.size())
		return m_values[i].values;
	throw std::runtime_error("No column with index " + i);
}

CSVType CSV::columnType(size_t i) const {
	if(i < m_values.size())
		return m_values[i].type;
	throw std::runtime_error("No column with index " + i);
}

void Stopwatch::start() {
	m_begin = std::chrono::steady_clock::now();
}

void Stopwatch::reset() {
	start();
}

std::string Stopwatch::time() {
	int t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()- m_begin).count();
	std::stringstream ss;
	ss << std::setfill('0');
	ss << std::setw(2) << (t / 3600) << ':';
	ss << std::setw(2) << (t / 60) % 60 << ':';
	ss << std::setw(2) << (t % 60);
	return ss.str();
}

std::string tt::util::getDriverFromPath(const std::string& path) {
	std::string ext = lowercase(extension(path));
	if(ext == ".tif" || ext == ".tiff") {
		return "GTiff";
	} else if(ext == ".shp") {
		return "ESRI Shapefile";
	} else if(ext == ".sqlite") {
		return "Spatialite";
	} else {
		return "";
	}
}

void tt::util::saveGrid(const std::string& file, const std::vector<double> grid,
		int cols, int rows, double minx, double miny, double xres, double yres,
		const std::string& proj) {
	GDALAllRegister();
	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* drv = dm->GetDriverByName("GTiff");
	GDALDataset* ds = drv->Create(file.c_str(), cols, rows, 1, GDT_Float32, 0);
	double trans[] = {minx, xres, 0, miny, 0, yres};
	ds->SetGeoTransform(trans);
	ds->SetProjection(proj.c_str());
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, (void*) grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
		std::cerr << "Failed to write to raster.\n";
	ds->GetRasterBand(1)->SetNoDataValue(-9999.0);
	GDALClose(ds);
}



/**
 * Waits for finalized collections of polygon parts or unioning.
 *
 * \param callback The callback functor, which accepts the id, geometry and context pointer.
 * \param pc The PolygonContext pointer.
 */
void tt::util::vec::polyMerge(PolyMergeCallback* callback, PolyCtx* pc) {

    std::vector<GEOSGeometry*> polys;
    GEOSGeometry* geom = nullptr;
    int id;

    while(pc->merging && !pc->geomBuf.empty()) {

        // Get the ID.
        id = pc->geomBuf.front().first;
        polys.swap(pc->geomBuf.front().second);
        pc->geomBuf.pop_front();

        if(polys.empty())
            continue;

        GEOSGeometry* multi = GEOSGeom_createCollection_r(pc->gctx, GEOS_GEOMETRYCOLLECTION, polys.data(), polys.size());
        geom = GEOSUnaryUnion_r(pc->gctx, multi);
        GEOSGeom_destroy_r(pc->gctx, multi);
        polys.clear();

        // If we're removing dangles, throw away all but the
        // largest single polygon. If it was originally a polygon, there are no dangles.
        int numGeoms;
        if(pc->removeDangles
                && (numGeoms = GEOSGetNumGeometries_r(pc->gctx, geom)) > 1) {
            int idx = 0;
            double a, area = 0;
            for(int i = 0; i < numGeoms; ++i) {
                const GEOSGeometry* p = GEOSGetGeometryN_r(pc->gctx, geom, i);
                GEOSArea_r(pc->gctx, p, &a);
                if(a > area) {
                    area = a;
                    idx = i;
                }
            }
            GEOSGeometry *g = GEOSGeom_clone_r(pc->gctx, GEOSGetGeometryN_r(pc->gctx, geom, idx)); // Force copy.
            GEOSGeom_destroy_r(pc->gctx, geom);
            geom = g;
        }

        // If we're removing holes, extract the exterior rings of all constituent polygons.
        if(pc->removeHoles) {
            std::vector<GEOSGeometry*> geoms0;
            for(int i = 0; i < GEOSGetNumGeometries_r(pc->gctx, geom); ++i) {
                const GEOSGeometry* p = GEOSGetGeometryN_r(pc->gctx, geom, i);
                const GEOSGeometry* l = GEOSGetExteriorRing_r(pc->gctx, p);
                const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(pc->gctx, l);
                GEOSGeometry* r = GEOSGeom_createLinearRing_r(pc->gctx, GEOSCoordSeq_clone_r(pc->gctx, seq));
                GEOSGeometry* npoly = GEOSGeom_createPolygon_r(pc->gctx, r, 0, 0);
                geoms0.push_back(npoly);
            }
            GEOSGeometry* g = GEOSGeom_createCollection_r(pc->gctx, GEOSGeomTypes::GEOS_MULTIPOLYGON, geoms0.data(), geoms0.size()); // Do not copy -- take ownership.
            GEOSGeom_destroy_r(pc->gctx, geom);
            geom = g;
        }

        // If the result is not a multi, make it one.
        if(GEOSGeomTypeId_r(pc->gctx, geom) != GEOSGeomTypes::GEOS_MULTIPOLYGON) {
            std::vector<GEOSGeometry*> geoms0;
            geoms0.push_back(geom);
            // Collection keeps ownership of the geom.
            GEOSGeometry* g = GEOSGeom_createCollection_r(pc->gctx, GEOSGeomTypes::GEOS_MULTIPOLYGON, geoms0.data(), 1);
            geom = g;
        }

        if(!geom) {
            std::cerr << "Warning: Null geometry." << std::endl;
        } else {
            (*callback)(id, geom, pc);
        }
    }
}


void tt::util::vec::PolyMergeCallback::operator()(int& id, GEOSGeometry* geom, PolyCtx* pc) {
	pc->geoms.emplace_back(id, geom);
}

GEOSGeometry* tt::util::vec::polyMakeGeom(GEOSContextHandle_t gctx, double x0, double y0, double x1, double y1, double eps, int dims) {

	// Build the geometry.
	GEOSCoordSequence* seq = GEOSCoordSeq_create_r(gctx, 5, dims);
	GEOSCoordSeq_setX_r(gctx, seq, 0, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 0, y0);
	GEOSCoordSeq_setX_r(gctx, seq, 1, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 1, y1);
	GEOSCoordSeq_setX_r(gctx, seq, 2, x1);
	GEOSCoordSeq_setY_r(gctx, seq, 2, y1);
	GEOSCoordSeq_setX_r(gctx, seq, 3, x1);
	GEOSCoordSeq_setY_r(gctx, seq, 3, y0);
	GEOSCoordSeq_setX_r(gctx, seq, 4, x0);
	GEOSCoordSeq_setY_r(gctx, seq, 4, y0);
	GEOSGeometry* ring = GEOSGeom_createLinearRing_r(gctx, seq);
	GEOSGeometry* poly = GEOSGeom_createPolygon_r(gctx, ring, 0, 0);
	GEOSGeometry* prec = GEOSGeom_setPrecision_r(gctx, poly, eps, 0);
	GEOSGeom_destroy_r(gctx, poly);
	return prec;
}


GEOSGeometry* tt::util::vec::setGeomZ(GEOSContextHandle_t gctx, const GEOSGeometry* geom, double z) {
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
        for(int i = 0; i < size; ++i) {
            GEOSCoordSeq_getX_r(gctx, cs, i, &x);
            GEOSCoordSeq_getY_r(gctx, cs, i, &y);
            GEOSCoordSeq_setX_r(gctx, cs0, i, x);
            GEOSCoordSeq_setY_r(gctx, cs0, i, y);
            GEOSCoordSeq_setZ_r(gctx, cs0, i, z);
        }
        return GEOSGeom_createLinearRing_r(gctx, cs0);
    } else{
        _runerr("Unexpected geom type: " << type);
    }
	return nullptr;
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
void tt::util::vec::makeCrownDataset(GDALDataset*& ds, OGRLayer*& layer,
		const std::string& filename, 
		const std::string& driver,
		const std::string& layerName,
		const std::string& projection, 
		OGRwkbGeometryType type,
		const std::string& idField, 
		const std::vector<std::pair<std::string, OGRFieldType>>& fields) {

	ds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, nullptr, nullptr, nullptr);

	std::string drvl = lowercase(driver);

	if(!ds) {
		// Get the vector driver.
		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		if(!drv)
			_runerr("Failed to find driver for " << driver << ".");

		// Create an output dataset for the polygons.
		char** dopts = NULL;
		if(drvl == "sqlite")
			dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

		ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
		CSLDestroy(dopts);
		if(!ds)
			_runerr("Failed to create dataset " << filename << ".");

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
		_runerr("Failed to create layer " << layerName << ".");

	// Create the field set for the database.
	for(const std::pair<std::string, OGRFieldType>& field : fields) {
		OGRFieldDefn dfn(field.first.c_str(), field.second);
		layer->CreateField(&dfn, field.first == idField); // Flag as the key if it matches idField.
	}

	// Dispose of the spatial reference object.
	if(sr)
		sr->Release();
}

