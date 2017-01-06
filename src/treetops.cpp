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
 *    maximum-value kernel to locate local maxima. Locates the tree top height
 *    from the original CHM.
 * 3) Delineate crowns. Uses the tree tops as seeds in a region-growing 
 *    algorithm that creates tree crown boundaries with configurable limits.
 *
 *  Created on: May 3, 2016
 *      Author: Rob Skelly 
 *       Email: rob@dijital.ca
 */

#include <queue>
#include <iostream>
#include <atomic>

#include <omp.h>

#include <boost/algorithm/string/join.hpp>

#include "geotools.hpp"
#include "sqlite.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "treetops.hpp"

using namespace geotools::raster;
using namespace geotools::util;
using namespace geotools::db;

using namespace geotools::treetops::config;
using namespace geotools::treetops::util;
using namespace geotools::treetops;

// Dummy cancel variable.
bool __tt_cancel = false;

namespace geotools {

	namespace treetops {

		namespace util {

			// Represents a grid cell, and maintains some properties of the
			// seed that originated it.
			class Node {
			public:
				uint64_t id;
				int c, r, tc, tr;
				double z, tz;

				Node(uint64_t id, int c, int r, double z, int tc, int tr, double tz) :
						id(id), c(c), r(r), tc(tc), tr(tr), z(z), tz(tz) {
				}

				Node(Point *p) :
					id(0), c(0), r(0), tc(0), tr(0), z(0), tz(0) {
					id = atoi(p->fields["id"].c_str());
					c = tc = atoi(p->fields["sc"].c_str());
					r = tr = atoi(p->fields["sr"].c_str());
					z = tz = atof(p->fields["sz"].c_str());
				}

			};

			// Returns true if the pixel at the center of the given circular window is
			// the maximum value in the window.
			bool isMaxCenter(Grid<float> &raster, int col, int row, int window, double *max) {
				*max = 0;
				int mc = 0, mr = 0;
				float v;
				float d = g_max(std::sqrt(2.0), g_sq((float) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && (v = raster.get(c, r)) > *max) {
							*max = v;
							mc = c;
							mr = r;
						}
					}
				}
				return mc == col && mr == row;
			}

			// Returns the max value of pixels in the kernel.
			void getKernelMax(Grid<float> &raster, int col, int row, int window,
				double *max, uint16_t *mc, uint16_t *mr) {
				*max = G_DBL_MAX_NEG;
				float v;
				float d = g_max(std::sqrt(2.0), g_sq((float) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && (v = raster.get(c, r)) > *max) {
							*max = v;
							*mc = c;
							*mr = r;
						}
					}
				}
			}

			// Set all cells in the circular kernel to zero except the center.
			void zeroKernel(Grid<uint8_t> &raster, int col, int row, int window) {
				float d = g_max(std::sqrt(2.0), g_sq((float) window / 2));
				for (int r = g_max(0, row - window / 2); r < g_min(raster.rows(), row + window / 2 + 1); ++r) {
					for (int c = g_max(0, col - window / 2); c < g_min(raster.cols(), col + window / 2 + 1); ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && c != col && r != row)
							raster.set(c, r, 0);
					}
				}
			}

			void savePoints(std::list<Top*> &tops, SQLite &db) {
				std::vector<Point*> points(tops.size());
				int i = 0;
				for (Top *t : tops) {
					std::map<std::string, std::string> fields;
					fields["id"] = std::to_string(t->id);
					fields["ox"] = std::to_string(t->ox);
					fields["oy"] = std::to_string(t->oy);
					fields["oz"] = std::to_string(t->oz);
					fields["sx"] = std::to_string(t->sx);
					fields["sy"] = std::to_string(t->sy);
					fields["sz"] = std::to_string(t->sz);
					fields["sc"] = std::to_string(t->sc);
					fields["sr"] = std::to_string(t->sr);
					points[i++] = new Point(t->sx, t->sy, t->oz, fields);
					delete t;
				}
				tops.clear();
				db.begin();
				db.addPoints(points);
				db.commit();
			}

		} // util

	} // trees

} // geotools


// TreetopsConfig implementation

TreetopsConfig::TreetopsConfig() :
	srid(0),
	buildIndex(false),
	tableCacheSize(1024 * 1024),
	rowCacheSize(24 * 1024 * 1024),
	threads(1),
	doSmoothing(false),
	smoothWindowSize(3),
	smoothSigma(0.8),
	doTops(false),
	doCrowns(false),
	crownsRadius(10.0),
	crownsHeightFraction(0.65),
	crownsMinHeight(4.0) {
}

void TreetopsConfig::checkSmoothing() const {
	if (!doSmoothing)
		g_argerr("Not configured to perform smoothing.");
	if (smoothOriginalCHM.empty())
		g_argerr("Smoothing: CHM filename must not be empty.");
	if (smoothSmoothedCHM.empty())
		g_argerr("Smoothing: Output filename must not be empty.");
	if (smoothSigma <= 0 || smoothSigma > 1)
		g_argerr("Smoothing: Std. deviation must be 0 < n <= 1. " << smoothSigma << " given.");
	if (smoothWindowSize % 2 == 0 || smoothWindowSize < 3)
		g_argerr("Smoothing: The window must be odd and >=3.");
}

void TreetopsConfig::checkTops() const {
	if (!doTops)
		g_argerr("Not configured to find treetops.");
	if (topsOriginalCHM.empty())
		g_argerr("Tops: Unsmoothed CHM filename must not be empty.");
	if (topsSmoothedCHM.empty())
		g_argerr("Tops: Smoothed CHM filename must not be empty.");
	if (topsTreetopsDatabase.empty())
		g_argerr("Tops: Treetops database filename must not be empty.");
	if (topsThresholds.empty()) {
		g_argerr("Tops: At least one threshold must be configured.");
	} else {
		double lastHeight;
		uint8_t lastWindow = 0;
		for(const auto &it : topsThresholds) {
			if(it.first < 0.0)
				g_argerr("Threshold heights below zero are not allowed.");
			if(it.second % 2 == 0 || it.second < 3)
				g_argerr("Window size must be odd and >= 3.");
			if(lastWindow) {
				if(lastWindow >= it.second)
					g_argerr("Each window must be larger than the previous one.");
				if(lastHeight >= it.first)
					g_argerr("Each height must be larger than the previous one.");
			}
			lastWindow = it.second;
			lastHeight = it.first;
		}
	}
}

void TreetopsConfig::checkCrowns() const {
	if (!doCrowns)
		g_argerr("Not configured to find crowns.");
	if (crownsRadius <= 0.0)
		g_argerr("Crowns: The maximum crown radius must be > 0. " << crownsRadius << " given.");
	if (crownsHeightFraction <= 0.0 || crownsHeightFraction > 1.0)
		g_argerr("Crowns: The crown height fraction must be between 0 and 1. " << crownsHeightFraction << " given.");
	if (crownsCrownsRaster.empty())
	g_argerr("Crowns: Output raster filename must not be empty.")
	if (crownsTreetopsDatabase.empty())
		g_argerr("Crowns: Treetops database filename must not be empty.");
	if (crownsSmoothedCHM.empty())
		g_argerr("Crowns: Smoothed CHM filename must not be empty.");
}

void TreetopsConfig::checkMerge() const {
	g_runerr("Not implemented.");
}

void TreetopsConfig::check() const {
	if (doSmoothing)
		checkSmoothing();
	if (doTops)
		checkTops();
	if (doCrowns)
		checkCrowns();
}

bool TreetopsConfig::canRun() const {
	try {
		check();
		return true;
	} catch (const std::exception &ex) {
		// Do nothing.
	}
	return false;
}

std::string TreetopsConfig::thresholds() const {
	std::vector<std::string> p(topsThresholds.size() * 2);
	int i = 0;
	for(const auto &it : topsThresholds) {
		p[i++] = std::to_string(it.first);
		p[i++] = std::to_string(it.second);
	}
	return boost::algorithm::join(p, ",");
}

void TreetopsConfig::parseThresholds(const std::string &str) {
	std::stringstream ss(str);
	std::string item1, item2;
	topsThresholds.clear();
	while(std::getline(ss, item1, ',')) {
		if(!std::getline(ss, item2, ','))
			break;
		if(item1.empty() || item2.empty()) continue;
		float height = atof(item1.c_str());
		uint8_t window = atoi(item2.c_str());
		if(window == 0) continue;
		topsThresholds[height] = window;
	}
}

// Top implementation

Top::Top(uint64_t id, double ox, double oy, double oz, double sx, double sy, double sz, int32_t sc, int32_t sr) :
	id(id), ox(ox), oy(oy), oz(oz), sx(sx), sy(sy), sz(sz), sc(sc), sr(sr) {
}

Top::Top() :
	id(0), ox(0), oy(0), oz(0), sx(0), sy(0), sz(0), sc(0), sr(0) {
}


// Treetops implementation

void Treetops::setCallbacks(Callbacks *callbacks) {
	m_callbacks = callbacks;
}

void Treetops::smooth(const TreetopsConfig &config, bool *cancel) {
	config.checkSmoothing();

	if (!cancel)
		cancel = &__tt_cancel;

	Raster<float> in(config.smoothOriginalCHM);
	Raster<float> out(config.smoothSmoothedCHM, 1, in);
	out.setNodata(in.nodata());
	in.smooth(out, config.smoothSigma, config.smoothWindowSize, m_callbacks, cancel);
}

void Treetops::treetops(const TreetopsConfig &config, bool *cancel) {

	config.checkTops();

	cancel = cancel == nullptr ? &__tt_cancel : cancel;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.01f);
		m_callbacks->statusCallback("Starting treetops...");
	}

	if (*cancel)
		return;

	if(m_callbacks)
		m_callbacks->statusCallback("Loading rasters...");

	// Initialize input rasters.
	Raster<float> original(config.topsOriginalCHM);
	Raster<float> smoothed(config.topsSmoothedCHM);

	if (*cancel)
		return;

	// Prepare database.
	if(m_callbacks)
		m_callbacks->statusCallback("Preparing database...");

	std::map<std::string, int> fields;
	fields["id"] = 1;
	fields["sx"] = 2;
	fields["sy"] = 2;
	fields["sz"] = 2;
	fields["sc"] = 1;
	fields["sr"] = 1;
	fields["ox"] = 2;
	fields["oy"] = 2;
	fields["oz"] = 2;
	SQLite db(config.topsTreetopsDatabase, SQLite::POINT, config.srid, fields, true);
	db.makeFast();
	db.dropGeomIndex();
	db.setCacheSize(config.tableCacheSize);

	if (*cancel)
		return;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Processing...");
	}

	MemRaster<uint8_t> topsGrid(original.cols(), original.rows(), true);
	topsGrid.fill(0);

	uint64_t total = config.topsThresholds.size() * topsGrid.rows();
	std::atomic<uint64_t> status(0);

	auto it = config.topsThresholds.begin();

	for(size_t i = 0; i < config.topsThresholds.size(); ++i) {

		uint8_t window = it->second;
		double threshold = it->first;

		it++;

		MemRaster<float> smooth(smoothed.cols(), 1024 + window, true);
		smooth.setNodata(smoothed.nodata());

		#pragma omp parallel for
		for(uint16_t j = 0; j < smoothed.rows() / 1024 + 1; ++j) {

			int32_t b = j * 1024;

			if(m_callbacks)
				m_callbacks->statusCallback("Reading...");

			uint16_t readOffset = b > 0 ? b - window / 2 : b;
			uint16_t writeOffset = b > 0 ? 0 : window / 2;
			#pragma omp critical(__tops_read)
			smoothed.readBlock(0, readOffset, smooth, 0, writeOffset);

			std::list<std::tuple<int32_t, int32_t, uint8_t> > tops;

			if(m_callbacks)
				m_callbacks->statusCallback("Processing...");

			for (int32_t row = window / 2; row < g_min(1024, smoothed.rows() - b) - window / 2; ++row) {
				if (*cancel) break;
				for (int32_t col = window / 2; col < smooth.cols() - window / 2; ++col) {

					double v = smooth.get(col, row);
					if(v < threshold)
						continue;

					double max;

					if (isMaxCenter(smooth, col, row, window, &max))
						tops.push_back(std::make_tuple(b + col, b + row - window / 2, window));
				}
				if (m_callbacks)
					m_callbacks->stepCallback(0.02f + (float) ++status / total * 0.48f);
			}

			if(m_callbacks)
				m_callbacks->statusCallback("Writing...");

			#pragma omp critical(__tops_write)
			{
				for (const auto &top : tops) {
					topsGrid.set(std::get<0>(top), std::get<1>(top), std::get<2>(top));
					zeroKernel(topsGrid, std::get<0>(top), std::get<1>(top), std::get<2>(top));
				}
			}
		}
	}

	Raster<uint8_t> tmp("/tmp/tmp.tif", 1, original);
	tmp.writeBlock(topsGrid);

	total = topsGrid.rows();
	status = 0;

	std::atomic<uint64_t> topId(0);

	#pragma omp parallel
	{

		std::list<Top*> tops;

		#pragma omp for
		for(int32_t row = 0; row < topsGrid.rows(); ++row) {
			if(*cancel) continue;
			if (m_callbacks)
				m_callbacks->statusCallback("Processing...");
			for(int32_t col = 0; col < topsGrid.cols(); ++col) {

				uint8_t window = topsGrid.get(col, row);
				if(!window) continue;

				// TODO: Can't use peak within window, because it may be on the slope of another crown
				// Get the original height from the unsmoothed raster.
				//double smax, omax;
				//uint16_t smc, smr, omc, omr;
				//getKernelMax(original, col, row, window, &omax, &omc, &omr);
				//getKernelMax(smoothed, col, row, window, &smax, &smc, &smr);

				Top * pt = new Top(++topId,
					original.toCentroidX(col), // omc; center of pixel
					original.toCentroidY(row), // omr
					original.get(col, row), // omax,
					original.toCentroidX(col), // smc; center of pixel
					original.toCentroidY(row), // smr
					smoothed.get(col, row), // smax,
					col, row // smc, smr
				);
				tops.push_back(pt);
			}

			if (m_callbacks)
				m_callbacks->stepCallback(0.5f + (float) ++status / total * 0.48f);

			if (tops.size() >= 5000) {
				if (m_callbacks)
					m_callbacks->statusCallback("Inserting points...");
				#pragma omp critical
				savePoints(tops, db);
			}
		}
		if (m_callbacks)
			m_callbacks->statusCallback("Inserting points...");
		#pragma omp critical
		savePoints(tops, db);
	}

	if (m_callbacks)
		m_callbacks->stepCallback(0.99f);

	if (*cancel)
		return;
	if (config.buildIndex) {
		if(m_callbacks)
			m_callbacks->statusCallback("Building index...");
		db.createGeomIndex();
	}
	db.makeSlow();

	if (m_callbacks) {
		m_callbacks->stepCallback(1.0);
		m_callbacks->statusCallback("Done.");
	}

}

void Treetops::treecrowns(const TreetopsConfig &config, bool *cancel) {
	config.checkCrowns();

	if (!cancel)
		cancel = &__tt_cancel;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.01f);
		m_callbacks->statusCallback("Starting crowns...");
	}

	if (*cancel)
		return;

	// Initialize the rasters.
	if(m_callbacks)
		m_callbacks->statusCallback("Preparing rasters...");

	Raster<float> inrast(config.crownsSmoothedCHM);
	Raster<uint32_t> outrast(config.crownsCrownsRaster, 1, inrast);
	outrast.setNodata(0);
	outrast.fill(0);

	MemRaster<float> smooth(inrast.cols(), inrast.rows(), 1);
	smooth.setNodata(inrast.nodata());
	smooth.writeBlock(inrast);

	MemRaster<uint32_t> blk(inrast.cols(), inrast.rows(), 1);
	blk.fill(0);

	double nodata = inrast.nodata();

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Preparing database...");
	}

	if (*cancel)
		return;

	// Initialize the database, get the treetop count and estimate the buffer size.

	SQLite db(config.crownsTreetopsDatabase);

	if (m_callbacks) {
		m_callbacks->stepCallback(0.03f);
		m_callbacks->statusCallback("Processing crowns...");
	}

	// Build the list of offsets for D8 or D4 search.
	size_t offsetCount = 8;
	int offsets[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
	//size_t offsetCount = 4;
	//int offsets[][2] = {{0, -1}, {-1, 0}, {1, 0}, {0, 1}};

	// To keep track of visited cells.
	std::vector<bool> visited((uint64_t) inrast.size());

	uint64_t geomCount = db.getGeomCount();
	int count = 100000;
	int offset = 0;

	if (m_callbacks)
		m_callbacks->statusCallback("Processing crowns...");

	while(true) {

		std::vector<Point*> tops;
		db.getPoints(tops, count, offset);
		offset += count;
		if(tops.empty())
			break;

		// Convert the Tops to Nodes.
		std::queue<Node*> q;
		for (Point *p : tops) {
			if (*cancel)
				break;
			q.push(new Node(p));
			delete p;
		}

		if (m_callbacks)
			m_callbacks->stepCallback(0.03f + ((float) (offset + count) / geomCount) * 0.97f);

		// Run through the queue.
		while (!*cancel && q.size()) {
			Node *n = q.front();
			q.pop();

			blk.set(n->c, n->r, (uint32_t) n->id);

			for(size_t i = 0; i < offsetCount; ++i) {
				int c = n->c + offsets[i][0];
				int r = n->r + offsets[i][1];

				if (r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols())
					continue;

				uint32_t idx = (uint64_t) r * inrast.cols() + c;
				if (visited[idx])
					continue;

				double v = smooth.get(idx);
				if (v != nodata           												// is not nodata
					&& v < n->z		      												// is less than the neighbouring pixel
					&& v >= config.crownsMinHeight 										// is greater than the min height
					&& (n->tz - v) / n->tz <= config.crownsHeightFraction 				// is greater than the threshold height
					&& g_sq(n->tc - c) + g_sq(n->tr - r) <= g_sq(config.crownsRadius) 	// is within the radius
				) {
					q.push(new Node(n->id, c, r, v, n->tc, n->tr, n->tz));
					blk.set(c, r, n->id);
					visited[idx] = true;
				}
			}
			delete n;
		}

		if (*cancel)
			continue;
	}

	// Filter out isolated pixels if D8 is used. See DTOOLS-22
	if(offsetCount == 9) {

		if(m_callbacks)
			m_callbacks->statusCallback("Filtering degenerate polygons...");

		for(int32_t row = 1; row < blk.rows() - 1; ++row) {
			for(int32_t col = 1; col < blk.cols() - 1; ++col) {
				uint32_t v = blk.get(col, row);
				if(blk.get(col, row - 1) != v
						&& blk.get(col - 1, row) != v
						&& blk.get(col + 1, row) != v
						&& blk.get(col, row + 1))
					blk.set(col, row, 0);
			}
		}
	}

	if(m_callbacks)
		m_callbacks->statusCallback("Writing output...");

	outrast.writeBlock(blk);

	if(!config.crownsCrownsDatabase.empty()) {
		if(m_callbacks)
			m_callbacks->statusCallback("Polygonizing...");
		outrast.polygonize(config.crownsCrownsDatabase, 1, m_callbacks, cancel);
	}

	if (m_callbacks) {
		m_callbacks->stepCallback(1.0);
		m_callbacks->statusCallback("Done.");
	}

}

void Treetops:: merge(const TreetopsConfig &config, bool *cancel) {
	config.checkMerge();

	// 1) Load tree tops into 2d tree (sqlite might be fine).
	// 2) Find pairs of tops within specified 3d distance
	// 3) Identify crowns to investigate using contained tops to locate
	// 4) Reload pixels from raster if necessary.

	// Things to investigate
	// - 3d proximity of tops
	// - circularity of combinations of crowns (i.e, if crowns are roughly
	//   circular, non-circular crowns may need to be merged)
	// -- if 2 or more tops are near, find the centroid between them
	//    and merge their crowns. if the crown meets some criterion for roundness
	//    like RMS, merge them.
	// - etc.


}


