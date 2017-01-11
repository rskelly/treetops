#include <cstdio>
#include <vector>
#include <memory>
#include <string>
#include <cassert>

#include "lasreader.hpp"

//#define NDEBUG

using namespace geotools::las;
using namespace geotools::util;
using namespace geotools::raster;

LASFilter::LASFilter() :
		m_radius(0),
		m_resX(0), m_resY(0),
		m_originX(0), m_originY(0) {
}

void LASFilter::setClasses(const std::set<uint8_t> &classes) {
	if(classes.empty()) {
		for(int i = 0; i < 255; ++i)
			m_classes[i] = 1;
	} else {
		for(int i = 0; i < 255; ++i)
			m_classes[i] = 0;
		for(const uint8_t cls : classes)
			m_classes[cls] = 1;
	}
}

void LASFilter::setGrid(double resX, double resY, double originX, double originY) {
	m_resX = resX;
	m_resY = resY;
	m_originX = originX;
	m_originY = originY;
}

void LASFilter::setRadius(double radius) {
	m_radius = g_sq(radius);
}

bool LASFilter::keep(const LASPoint &pt) const{
	if(!m_classes[pt.cls]) return false;
	if(m_radius) {
		assert(m_resX != 0 && m_resY != 0);
		int col = (int) ((pt.x - m_originX) / m_resX);
		int row = (int) ((pt.y - m_originY) / m_resY);
		double gx = col * m_resX + m_resX * 0.5;
		double gy = row * m_resY + m_resY * 0.5;
		double dist = g_sq(gx - pt.x) + g_sq(gy - pt.y);
		if(dist > m_radius) return false;
	}
	return true;
}

LASReaderCallback::~LASReaderCallback() {}

void _read(void *buf, size_t l1, size_t l2, std::FILE *f) {
	if (!std::fread(buf, l1, l2, f))
		g_runerr("Failed to read " << l1 << "x" << l2 << " from LAS file.");
}

void LASReader::init() {
	//g_debug("LASReader: open " << m_file);
	m_f = std::fopen(m_file.c_str(), "rb");
	if (!m_f)
		g_runerr("Failed to open " << m_file << "; errno: " << errno);

	std::fseek(m_f, 4, SEEK_SET);
	_read((void *) &m_sourceId, sizeof (uint16_t), 1, m_f);

	std::fseek(m_f, 24, SEEK_SET);
	char maj, min;
	_read((void *) &maj, sizeof (char), 1, m_f);
	_read((void *) &min, sizeof (char), 1, m_f);
	//m_version = "" + maj + "." + min;

	std::fseek(m_f, 94, SEEK_SET);
	_read((void *) &m_headerSize, sizeof (uint16_t), 1, m_f);

	std::fseek(m_f, 96, SEEK_SET);
	_read((void *) &m_offset, sizeof (uint32_t), 1, m_f);

	std::fseek(m_f, 104, SEEK_SET);
	_read((void *) &m_pointFormat, sizeof (uint8_t), 1, m_f);
	_read((void *) &m_pointLength, sizeof (uint16_t), 1, m_f);

	uint32_t pc, pcbr[5];
	_read((void *) &pc, sizeof (uint32_t), 1, m_f);
	_read((void *) &pcbr, 5 * sizeof (uint32_t), 1, m_f);
	m_pointCount = pc;
	for (int i = 0; i < 5; ++i)
		m_pointCountByReturn[i] = pcbr[i];

	_read((void *) &m_xScale, sizeof (double), 1, m_f);
	_read((void *) &m_yScale, sizeof (double), 1, m_f);
	_read((void *) &m_zScale, sizeof (double), 1, m_f);

	_read((void *) &m_xOffset, sizeof (double), 1, m_f);
	_read((void *) &m_yOffset, sizeof (double), 1, m_f);
	_read((void *) &m_zOffset, sizeof (double), 1, m_f);

	_read((void *) &m_xMax, sizeof (double), 1, m_f);
	_read((void *) &m_xMin, sizeof (double), 1, m_f);
	_read((void *) &m_yMax, sizeof (double), 1, m_f);
	_read((void *) &m_yMin, sizeof (double), 1, m_f);
	_read((void *) &m_zMax, sizeof (double), 1, m_f);
	_read((void *) &m_zMin, sizeof (double), 1, m_f);

	// TODO: Extended point count.

	LASPoint::setScale(m_xScale, m_yScale, m_zScale);
	reset();
}

LASReader::LASReader(const std::string &file) :
	m_f(nullptr),
	m_file(file),
	m_filter(nullptr) {
	init();
}

LASReader::~LASReader() {
	if (m_f)
		std::fclose(m_f);
	if(m_filter)
		delete m_filter;
	if (m_buf.get())
		delete m_buf.release();
}

void LASReader::setFilter(const LASFilter *filter) {
	m_filter = new LASFilter(*filter);
}

std::string LASReader::filename() const {
	return m_file;
}

void LASReader::reset() {
	m_curPoint = 0;
	std::fseek(m_f, m_offset, SEEK_SET);
}

bool LASReader::loadBatch() {
	if (!m_f || m_curPoint >= m_pointCount)
		return false;
	m_batchSize = g_min(BATCH_SIZE, m_pointCount - m_curPoint);
	if (m_buf.get())
		delete m_buf.release();
	m_buf.reset(new Buffer(m_batchSize * m_pointLength));
	//g_debug(" -- loading " << m_batchSize << "; " << m_pointCount << "; " << m_pointLength);
	try {
		_read((void *) m_buf->buf, m_batchSize * m_pointLength, 1, m_f);
	} catch(const std::exception &ex) {
		g_debug(" -- " << ex.what());
		g_runerr("Failed to read batch. LAS file may be shorter than indicated by header: " << ex.what());
	}
	return true;
}

bool LASReader::isValid() const {
	Bounds b = bounds();
	if(!(b.width() > 0 && b.height() > 0)) return false;
	if(!m_pointCount) return false;
	return true;
}

bool LASReader::next(LASPoint &pt) {
	while(true) {
		if (m_curPoint >= m_pointCount)
			return false;
		if (m_curPoint % BATCH_SIZE == 0 && !loadBatch())
			return false;
		char *buf = ((char *) m_buf->buf) + (m_curPoint % BATCH_SIZE) * m_pointLength;
		pt.readLAS(buf, m_pointFormat);
		++m_curPoint;
		if(m_filter && !m_filter->keep(pt))
			continue;
		return true;
	}
}

Bounds LASReader::bounds() const {
	return Bounds(m_xMin, m_yMin, m_xMax, m_yMax, m_zMin, m_zMax);
}

uint64_t LASReader::pointCount() const {
	return m_pointCount;
}
   

bool __lr__cancel = false;

void LASMultiReader::init(const std::vector<std::string> &files) {
	m_bounds.collapse();
	for(const std::string &file : files) {
		if(*m_cancel) return;
		LASReader r(file);
		if(!r.isValid()) {
			g_warn("Invalid LAS file: " << file);
			continue;
		}
		m_files.push_back(file);
		m_bounds.extend(r.bounds());
		m_blockBounds.push_back(r.bounds());
		m_pointCount += r.pointCount();
	}
	m_cols = bounds().maxCol(m_resolutionX) + 1;
}

LASMultiReader::LASMultiReader(const std::vector<std::string> &files,
		double resolutionX, double resolutionY, bool *cancel) :
	m_cancel(cancel),
	m_idx(0),
	m_cols(0),
	m_cellCount(0),
	m_pointCount(0),
	m_resolutionX(resolutionX), m_resolutionY(resolutionY),
	m_reader(nullptr),
	m_filter(nullptr) {

	if(m_cancel == nullptr)
		m_cancel = &__lr__cancel;

	init(files);
}
        
LASMultiReader::~LASMultiReader() {
	if(m_reader)
		delete m_reader;
	if(m_filter)
		delete m_filter;
}

void LASMultiReader::setFilter(const LASFilter *filter) {
	m_filter = new LASFilter(*filter);
}

void LASMultiReader::buildFinalizer(const LASReaderCallback *callback) {
	try {
		int cols = m_cols = m_bounds.maxCol(m_resolutionX) + 1;
		int rows =          m_bounds.maxRow(m_resolutionY) + 1;
		bool fileChanged;
		uint32_t file = 0;
		LASPoint pt;

		// Resize and zero the finalizer.
		m_finalizer.resize((size_t) cols * rows);
		std::fill(m_finalizer.begin(), m_finalizer.end(), 0);

		// Count the points in each cell.
		while(next(pt, nullptr, nullptr, &fileChanged)) {
			if(*m_cancel) return;
			int col = m_bounds.toCol(pt.x, m_resolutionX);
			int row = m_bounds.toRow(pt.y, m_resolutionY);
			m_finalizer[(size_t) row * cols + col]++;
			if(fileChanged && callback)
				callback->status((float) ++file / m_files.size());
		}

		// Count the valid cells.
		m_cellCount = 0;
		for(const uint32_t &c : m_finalizer)
			if(c) ++m_cellCount;

	} catch(const std::exception &ex) {
		m_finalizer.resize(0);
		g_runerr("Failed to initialize finalizer. LAS bounds may be incorrect in header: " << ex.what());
	}
	reset();
}

void LASMultiReader::reset() {
	if(m_reader)
		delete m_reader;
	m_reader = nullptr;
	m_idx = 0;
}

uint64_t LASMultiReader::cellCount() const {
	return m_cellCount;
}
    
uint64_t LASMultiReader::pointCount() const {
	return m_pointCount;
}


void LASMultiReader::setBounds(const Bounds &bounds) {
	m_bounds.collapse();
	m_bounds.extend(bounds);
}
    
Bounds LASMultiReader::bounds() const {
	return m_bounds;
}
    
std::vector<Bounds> LASMultiReader::blockBounds() const {
	return m_blockBounds;
}
    
bool LASMultiReader::next(LASPoint &pt, bool *final, uint64_t *finalIdx, bool *fileChanged) {
	if(!m_reader || !m_reader->next(pt)) {
		if(m_idx >= m_files.size())
			return false;
		if(m_reader)
			delete m_reader;
		m_reader = new LASReader(m_files[m_idx++]);
		m_reader->setFilter(m_filter);
		if(!m_reader->next(pt))
			return false;
		// If the file has been switched, set the indicator to true
		if(fileChanged != nullptr)
			*fileChanged = true;
	} else if(fileChanged != nullptr) {
		// If the file is not switched, set the indicator to false.
		*fileChanged = false;
	}
	// If a pointer is given for finalIdx, set it.
	if(finalIdx != nullptr) {
		int col = m_bounds.toCol(pt.x, m_resolutionX);
		int row = m_bounds.toRow(pt.y, m_resolutionY);
		size_t idx = (size_t) row * m_cols + col;
		uint32_t count = m_finalizer[idx];
		if(count == 0)
			g_runerr("Finalizer reached zero. This is impossible!");
		m_finalizer[idx]--;
		if(count == 1) {
			*finalIdx = idx;
			if(final != nullptr)
				*final = true;
		} else {
			if(final != nullptr)
				*final = false;
		}
	}
	return true;
}
