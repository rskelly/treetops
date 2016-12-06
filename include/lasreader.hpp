#ifndef __LASREADER_HPP__
#define __LASREADER_HPP__

#include <cstdio>
#include <vector>
#include <memory>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "util.hpp"
#include "laspoint.hpp"
#include "raster.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_with_holes_2<K>                       Polygon_with_holes_2;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef K::Point_2                                          Point_2;

using namespace geotools::las;
using namespace geotools::util;
using namespace geotools::raster;

class LASReader {
private:
    const static uint64_t BATCH_SIZE = 1000000;
    std::FILE *m_f;
    std::string m_file;
    uint16_t m_sourceId;
    std::string m_version;
    uint16_t m_headerSize;
    uint32_t m_offset;
    uint8_t m_pointFormat;
    uint16_t m_pointLength;
    uint64_t m_pointCount;
    uint64_t m_pointCountByReturn[5];
    double m_xScale;
    double m_yScale;
    double m_zScale;
    double m_xOffset;
    double m_yOffset;
    double m_zOffset;
    double m_xMin;
    double m_yMin;
    double m_zMin;
    double m_xMax;
    double m_yMax;
    double m_zMax;

    uint64_t m_curPoint;
    uint64_t m_batchSize;

    std::unique_ptr<Buffer> m_buf;
    std::queue<LASPoint> m_pts;

    void _read(void *buf, size_t l1, size_t l2, std::FILE *f) {
        if (!std::fread(buf, l1, l2, f))
            g_runerr("Failed to read " << l1 << "x" << l2 << " from LAS file.");
    }

    void load() {
        m_f = std::fopen(m_file.c_str(), "rb");
        if (!m_f)
            g_runerr("Failed to open " << m_file);

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

        g_debug(" -- " << m_file << " has " << m_pointCount << " points");
        LASPoint::setScale(m_xScale, m_yScale, m_zScale);
        reset();
    }

public:

    LASReader(const std::string &file) :
        m_f(nullptr),
        m_file(file) {
        load();
    }

    ~LASReader() {
        if (m_f)
            std::fclose(m_f);
        if (m_buf.get())
            delete m_buf.release();
    }

    void reset() {
        m_curPoint = 0;
        std::fseek(m_f, m_offset, SEEK_SET);
    }

    bool loadBatch() {
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
            g_runerr("Failed to read batch. LAS file may be shorter than indicated by header.");
        }
        return true;
    }

    bool next(LASPoint &pt) {
        if (m_curPoint >= m_pointCount)
            return false;
        if (m_curPoint % BATCH_SIZE == 0 && !loadBatch())
            return false;
        char *buf = ((char *) m_buf->buf) + (m_curPoint % BATCH_SIZE) * m_pointLength;
        pt.readLAS(buf, m_pointFormat);
        ++m_curPoint;
        return true;
    }

    Bounds bounds() {
        return Bounds(m_xMin, m_yMin, m_xMax, m_yMax, m_zMin, m_zMax);
    }

    uint64_t pointCount() {
        return m_pointCount;
    }
   

};

class LASMultiReader {
private:
    std::vector<std::string> m_files;
    std::vector<std::unique_ptr<LASReader> > m_readers;
    LASReader *m_reader;
    uint32_t m_idx;
    uint32_t m_cols;
    std::unique_ptr<MemRaster<uint32_t> > m_finalizer;
    double m_resolutionX;
    double m_resolutionY;
    Bounds m_bounds;
    std::vector<Bounds> m_blockBounds;
    bool m_cancel;
    uint64_t m_cells;
    
    void load(const std::vector<std::string> &files) {
        m_bounds.collapse();
        for(const std::string &file : files) {
            if(m_cancel) return;
            std::unique_ptr<LASReader> r(new LASReader(file));
            m_files.push_back(file);
            m_bounds.extend(r->bounds());
            m_blockBounds.push_back(r->bounds());
            m_readers.push_back(std::move(r));
        }
        m_cols = bounds().maxCol(m_resolutionX) + 1;
    }

    void buildFinalizer() {
        if(m_finalizer.get())
            delete m_finalizer.release();
        int cols = m_bounds.maxCol(m_resolutionX) + 1;
        int rows = m_bounds.maxRow(m_resolutionY) + 1;
        m_cols = cols;
        g_debug(" -- finalizer: " << cols << ", " << rows);
        m_finalizer.reset(new MemRaster<uint32_t>(cols, rows, false));
        m_finalizer->fill(0);
        LASPoint pt;
        std::set<uint64_t> cells;
        try {
            while(next(pt, nullptr, nullptr)) {
                if(m_cancel) return;
                int col = m_bounds.toCol(pt.x, m_resolutionX);
                int row = m_bounds.toRow(pt.y, m_resolutionY);
                m_finalizer->set(col, row, m_finalizer->get(col, row) + 1);
                cells.insert((uint64_t) row * m_cols + col);
            }
        } catch(const std::exception &ex) {
            g_runerr("Failed to initialize finalizer. LAS bounds may be incorrect in header.");
        }
        g_debug(" -- finalizer: " << m_finalizer->cols() << ", " << m_finalizer->rows() << "; " << m_bounds.print());
        m_cells = cells.size();
        reset();
    }
    
public:
    LASMultiReader(const std::vector<std::string> &files, double resolutionX, double resolutionY) :
        m_reader(nullptr),
        m_idx(0),
        m_resolutionX(resolutionX), m_resolutionY(resolutionY),
        m_cancel(false),
        m_cells(0) {
        load(files);
    }
        
    void init() {
        buildFinalizer();
    }

    void reset() {
        m_idx = 0;
        m_reader = nullptr;
        m_cancel = false;
        for(const auto &r : m_readers)
            r->reset();
    }
    
    uint64_t cellCount() {
        return m_cells;
    }
    
    void cancel() {
        m_cancel = true;
    }
    
    void setBounds(const Bounds &bounds) {
        m_bounds.collapse();
        m_bounds.extend(bounds);
    }
    
    Bounds bounds() const {
        return m_bounds;
    }
    
    std::vector<Bounds> blockBounds() const {
        return m_blockBounds;
    }
    
    bool next(LASPoint &pt, bool *final, uint64_t *finalIdx) {
        if(!m_reader || !m_reader->next(pt)) {
            if(m_idx >= m_readers.size())
                return false;
            if(m_reader)
                m_reader->reset();
            m_reader = m_readers[m_idx++].get();
            if(!m_reader->next(pt))
                return false;
        }
        if(finalIdx != nullptr) {
            int col = m_bounds.toCol(pt.x, m_resolutionX);
            int row = m_bounds.toRow(pt.y, m_resolutionY);
            uint32_t count = m_finalizer->get(col, row);
            //g_debug(" -- pt " << pt.x << ", " << pt.y << "; " << col << ", " << row << "; " << count << "; " << m_cols);
            if(count == 0)
                g_runerr("Finalizer reached zero. This is impossible!");
            m_finalizer->set(col, row, count - 1);
            if(count == 1) {
                *finalIdx = (uint64_t) row * m_cols + col;
                *final = true;
                //g_debug(" -- final " << *finalIdx);
            } else {
                *final = false;
            }
        }
        return true;
    }
};
#endif
