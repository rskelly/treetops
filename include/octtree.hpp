/*
 * octree.hpp
 *
 *  Created on: Jan 10, 2017
 *      Author: rob
 */

#ifndef __OCTREE_HPP__
#define __OCTREE_HPP__

#include "util.hpp"
#include "laspoint.hpp"

using namespace geotools::las;
using namespace geotools::util;

#define FILE_SIZE 1000000

class OctMemManager {
private:
	std::map<uint64_t, uint64_t> m_offsets;
	uint64_t m_id;
	std::queue<uint64_t> m_released;
	std::map<uint64_t, MappedFile*> m_mappedFiles;

public:
	OctMemManager() :
		m_id(0) {
	}

	~OctMemManager() {
		for(auto &it : m_mappedFiles)
			delete it.second;
	}

	MappedFile *getMappedFile(uint64_t offset) {
		uint64_t idx = offset / FILE_SIZE;
		if(m_mappedFiles.find(idx) == m_mappedFiles.end())
			m_mappedFiles[idx] = Util::mapFile(Util::tmpFile(), FILE_SIZE, true);
		return m_mappedFiles[idx];
	}

	void writePoint(uint64_t id, const LASPoint &pt) {
		uint64_t offset = id * 1000 * LASPoint::dataSize() + m_offsets[id] * LASPoint::dataSize();
		MappedFile *mf = getMappedFile(offset);
		pt.write(mf->data() + offset);
	}

	void readPoint(uint64_t id, LASPoint &pt) {
		uint64_t offset = id * 1000 * LASPoint::dataSize() + m_offsets[id] * LASPoint::dataSize();
		MappedFile *mf = getMappedFile(offset);
		pt.read(mf->data() + offset);
	}

	void readPoints(uint64_t id, std::vector<LASPoint*> &pts) {
		uint64_t start = id * 1000 * LASPoint::dataSize();
		uint64_t end = id * 1000 * LASPoint::dataSize() + m_offsets[id] * LASPoint::dataSize();
		MappedFile *mf = getMappedFile(start);
		for(uint64_t offset = start; offset < end; ++offset) {
			LASPoint *pt = new LASPoint();
			pt->read(mf->data() + offset);
			pts.push_back(pt);
		}
	}

	uint64_t nextId() {
		uint64_t id;
		if(!m_released.empty()) {
			id = m_released.front(); m_released.pop();
		} else {
			id = ++m_id;
			m_offsets[id] = 0;
		}
		return id;
	}

	void releaseId(uint64_t id) {
		m_released.push(id);
		m_offsets[id] = 0;
	}
};

class OctNode {
private:
	double m_midx, m_midy;
	Bounds m_bounds;
	OctNode *m_nodes[4];
	uint64_t m_count;
	uint64_t m_mid;
	OctMemManager *m_manager;
	bool m_split;

	// Get a node for the quadrant.
	OctNode *getNode(uint8_t idx) {
		if(!m_nodes[idx]) {
			Bounds b;
			b.collapse();
			switch(idx) {
			case 0:
				b.extend(m_bounds.minx(), m_bounds.miny());
				b.extend(m_bounds.minx() + m_bounds.width() / 2.0, m_bounds.miny() + m_bounds.height() / 2.0);
				break;
			case 1:
				b.extend(m_bounds.minx(), m_bounds.miny() + m_bounds.height() / 2.0);
				b.extend(m_bounds.minx() + m_bounds.width() / 2.0, m_bounds.maxy());
				break;
			case 2:
				b.extend(m_bounds.minx() + m_bounds.width() / 2.0, m_bounds.miny());
				b.extend(m_bounds.maxx(), m_bounds.miny() + m_bounds.height() / 2.0);
				break;
			case 3:
				b.extend(m_bounds.minx() + m_bounds.width() / 2.0, m_bounds.miny() + m_bounds.height() / 2.0);
				b.extend(m_bounds.maxx(), m_bounds.maxy());
				break;
			}
			m_nodes[idx] = new OctNode(b);
		}
		return m_nodes[idx];
	}

	void split() {
		m_split = true;
		std::vector<LASPoint*> pts;
		m_manager->readPoints(m_mid, pts);
		m_manager->releaseId(m_mid);
		for(const LASPoint *pt : pts) {
			addPoint(*pt);
			delete pt;
		}
	}

public:
	OctNode(Bounds &bounds, OctMemManager *manager) :
		m_count(0),
		m_mid(0),
		m_manager(manager),
		m_split(false) {

		m_nodes = {nullptr, nullptr, nullptr, nullptr};
		m_bounds.collapse();
		m_bounds.extend(bounds);
		m_midx = m_bounds.minx() + m_bounds.width() / 2.0;
		m_midy = m_bounds.miny() + m_bounds.height() / 2.0;
		m_mid = m_manager->nextId();
	}

	void addPoint(const LASPoint &pt) {
		if(m_count < 1000) {
			m_manager->writePoint(m_mid, pt);
		} else {
			if(!m_split)
				split();
			if(pt.x < m_midx) {
				if(pt.y < m_midy) {
					getNode(0)->addPoint(pt);
				} else {
					getNode(1)->addPoint(pt);
				}
			} else {
				if(pt.y < m_midy) {
					getNode(2)->addPoint(pt);
				} else {
					getNode(3)->addPoint(pt);
				}
			}
		}
	}
};

class OctTree : public OctNode {
	OctNode *m_nodes[4]; // Clockwise, zero northwest.

public:
	OctTree(Bounds &bounds) :
		m_manager(new OctMemManager()),
		OctNode(bounds, m_manager) {
	}

	~OctTree() {
		delete m_manager;
	}
};


#endif /* __OCTREE_HPP__ */
