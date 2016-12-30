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

				Node(const Top &top) :
						id(top.id), c(top.col), r(top.row), tc(top.col), tr(top.row), z(
								top.z), tz(top.z) {
				}
			};

			// Returns true if the pixel at the center of the given raster is
			// the maximum value in the raster.
			bool isMaxCenter(MemRaster<float> &raster, int col, int row, int window,
					double *max) {
				int cc = col + window / 2;
				int cr = row + window / 2;
				float nd = raster.nodata();
				if (raster.get(cc, cr) == nd)
					return false;
				*max = 0;
				int mc = 0, mr = 0;
				for (int r = row; r < row + window; ++r) {
					for (int c = col; c < col + window; ++c) {
						float v = raster.get(c, r);
						if (v != nd && v > *max) {
							*max = v;
							mc = c;
							mr = r;
						}
					}
				}
				return mc == cc && mr == cr;
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
	topsMinHeight(4.0),
	topsWindowSize(7),
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
	if (topsWindowSize % 2 == 0 || topsWindowSize < 3)
		g_argerr("Tops: Treetops window size must be an odd number >= 3. " << topsWindowSize << " given.");
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


// Top implementation

Top::Top(uint64_t id, double x, double y, double z, double uz, int col, int row) :
		id(id), x(x), y(y), z(z), uz(uz), col(col), row(row) {
}

Top::Top() :
		id(0), x(0), y(0), z(0), uz(0), col(0), row(0) {
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

	// Initialize input rasters.
	g_debug(" -- tops: opening rasters");
	Raster<float> original(config.topsOriginalCHM);
	Raster<float> smoothed(config.topsSmoothedCHM);

	if (*cancel)
		return;

	// Prepare database.
	if(m_callbacks)
		m_callbacks->statusCallback("Preparing database...");

	std::map<std::string, int> fields;
	fields["id"] = 1;
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

	int32_t bufSize = 256; // The number of rows each thread works on at a time.
	std::atomic<uint64_t> topCount(0);
	std::atomic<int32_t> curRow(0);     // Used for status indicator.
	std::atomic<uint64_t> topId(0);

	#pragma omp parallel
	{
		MemRaster<float> blk(original.cols(), bufSize + config.topsWindowSize);
		blk.setNodata(original.nodata());
		blk.fill(blk.nodata());

		std::map<uint64_t, std::unique_ptr<Top> > tops0;
		int32_t bufSize0;

		#pragma omp for
		for (int32_t brow = 0; brow < original.rows(); brow += bufSize) {
			if (*cancel)
				continue;

			bufSize0 = g_min(bufSize,
					original.rows() - brow - config.topsWindowSize);
			if (bufSize0 < bufSize)
				blk.init(original.cols(), bufSize0 + config.topsWindowSize);

			if (*cancel)
				continue;

			#pragma omp critical(__c)
			smoothed.readBlock(0, brow, blk);

			for (int32_t row = 0; row < bufSize0; ++row) {
				if (*cancel)
					break;
				for (int32_t col = 0; col < original.cols() - config.topsWindowSize; ++col) {
					if (*cancel)
						break;

					int32_t r = row + config.topsWindowSize / 2;
					int32_t c = col + config.topsWindowSize / 2;
					double max;

					if (blk.get(c, r) >= config.topsMinHeight
							&& isMaxCenter(blk, col, row, config.topsWindowSize, &max)) {

						// Compute the id based on the cell.
						uint64_t id = ((uint64_t) c << 32) | r;
						// Get the original height from the unsmoothed raster.
						double umax = original.get(c, r + brow, 1);

						std::unique_ptr<Top> pt(
							new Top(++topId, original.toCentroidX(c), // center of pixel
							original.toCentroidY(r + brow), max, umax, c, r + brow)
						);
						tops0[id] = std::move(pt);
					}
				}
			}

			if (m_callbacks) {
				curRow += bufSize0;
				m_callbacks->stepCallback(0.02f + (float) curRow / original.rows() * 0.48f);
			}
		}

		uint64_t topCount0 = tops0.size();
		uint64_t b = 0;
		uint64_t batch = db.maxAddPointCount();

		topCount += topCount0;

		std::vector<std::unique_ptr<Point> > points;
		for (const auto &it : tops0) {
			if (*cancel)
				break;
			const std::unique_ptr<Top> &t = it.second;
			std::map<std::string, std::string> fields;
			fields["id"] = std::to_string(t->id);
			std::unique_ptr<Point> pt(new Point(t->x, t->y, t->uz, fields));
			points.push_back(std::move(pt));

			if (++b % batch == 0 || b >= topCount0) {
				if (*cancel)
					break;
				if(m_callbacks)
					m_callbacks->statusCallback("Inserting points...");
				#pragma omp critical(__b)
				{
					db.begin();
					db.addPoints(points);
					db.commit();
				}
				points.clear();
				if (m_callbacks)
					m_callbacks->stepCallback((float) b / topCount0 * 0.3f + 0.48f);
			}
			if(m_callbacks)
				m_callbacks->statusCallback("Processing...");
		}
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

	uint32_t batchSize = 10000;

	// Initialize the rasters.
	if(m_callbacks)
		m_callbacks->statusCallback("Preparing rasters...");

	Raster<float> inrast(config.crownsSmoothedCHM);
	Raster<uint32_t> outrast(config.crownsCrownsRaster, 1, inrast);
	outrast.setNodata(0, 1);
	outrast.fill(0, 1);

	double nodata = inrast.nodata();

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Preparing database...");
	}

	if (*cancel)
		return;

	// Initialize the database, get the treetop count and estimate the buffer size.
	uint64_t geomCount;
	SQLite db(config.crownsTreetopsDatabase);
	db.getGeomCount(&geomCount);

	if (m_callbacks) {
		m_callbacks->stepCallback(0.03f);
		m_callbacks->statusCallback("Processing crowns...");
	}

	// Build the list of offsets for D8 or D4 search.
	int offsets[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

	// The number of extra rows above and below the buffer.
	int32_t bufRows = (int) std::ceil(g_abs(config.crownsRadius / inrast.resolutionY()));
	// The height of the row, not including disposable buffer. Use bufRows as lower bound to
	// avoid read error later (keeps min row index to >=0)
	int32_t rowStep = g_max(bufRows,
			(int) g_abs(std::ceil((double) batchSize / geomCount * inrast.rows()) / inrast.resolutionY()));
	// The total height of the buffer
	int32_t rowHeight = rowStep + bufRows * 2;
	std::atomic<uint32_t> curRow(0);

	#pragma omp parallel
	{
		// To keep track of visited cells.
		std::vector<bool> visited((uint64_t) inrast.cols() * rowHeight);
		MemRaster<float> buf(inrast.cols(), rowHeight);
		MemRaster<uint32_t> blk(inrast.cols(), rowHeight);
		buf.fill(inrast.nodata());
		blk.fill(0);

		#pragma omp for
		for (int32_t row = 0; row < inrast.rows(); row += rowStep) {
			if (*cancel)
				continue;
			curRow += rowStep;

			// Load the tree tops for the strip.
			Bounds bounds(inrast.toX(0), inrast.toY(row - bufRows),
					inrast.toX(inrast.cols()),
					inrast.toY(row + rowStep + bufRows));
			std::vector<std::unique_ptr<Point> > tops;

			#pragma omp critical(__crowns_getpoints)
			db.getPoints(tops, bounds);

			#pragma omp critical(__crowns_readbuf)
			inrast.readBlock(0, row == 0 ? row : row - bufRows, buf, 0, row == 0 ? bufRows : 0);

			// Convert the Tops to Nodes.
			std::queue<std::unique_ptr<Node> > q;
			for (const std::unique_ptr<Point> &t : tops) {
				if (*cancel)
					break;
				int col = inrast.toCol(t->x);
				int row = inrast.toRow(t->y);
				int id = atoi(t->fields["id"].c_str());
				std::unique_ptr<Node> nd(new Node(id, col, row, t->z, col, row, t->z));
				q.push(std::move(nd));
			}

			if (m_callbacks) {
				m_callbacks->stepCallback(0.03f + ((float) curRow / inrast.rows()) * 0.97f);
				m_callbacks->statusCallback("Processing crowns...");
			}

			// Run through the queue.
			while (!*cancel && q.size()) {
				std::unique_ptr<Node> n = std::move(q.front());
				q.pop();

				blk.set(n->c, n->r - row + bufRows, (uint32_t) n->id);

				for(size_t i = 0; i < 8; ++i) {
					int c = n->c + offsets[i][0];
					int r = n->r + offsets[i][1];

					if (r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols())
						continue;
					if (r - row + bufRows < 0 || r - row + bufRows >= buf.rows())
						continue;

					uint32_t idx = (uint64_t) (r - row + bufRows) * inrast.cols() + c;
					if (visited[idx])
						continue;

					double v = buf.get(idx);
					if (v != nodata           												// is not nodata
						&& v < n->z		      												// is less than the neighbouring pixel
						&& v >= config.crownsMinHeight 										// is greater than the min height
						&& (v / n->tz) >= config.crownsHeightFraction 						// is greater than the threshold height
						&& g_sq(n->tc - c) + g_sq(n->tr - r) <= g_sq(config.crownsRadius) 	// is within the radius
					) {
						std::unique_ptr<Node> nd(new Node(n->id, c, r, v, n->tc, n->tr, n->tz));
						q.push(std::move(nd));
						visited[idx] = true;
					}
				}
			}

			if (row > 0 && (row + bufRows) >= inrast.rows())
				continue;

			if (*cancel)
				continue;

			if(m_callbacks)
				m_callbacks->statusCallback("Writing output...");

			#pragma omp critical(__b)
			outrast.writeBlock(0, row, blk, 0, bufRows);

		}
	}

	if(!config.crownsCrownsDatabase.empty()) {
		if(m_callbacks)
			m_callbacks->statusCallback("Polygonizing...");
		outrast.polygonize(config.crownsCrownsDatabase, 1);
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

