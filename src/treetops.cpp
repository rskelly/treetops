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

				Node(const Top &top) :
						id(top.id), c(top.col), r(top.row), tc(top.col), tr(top.row), z(
								top.z), tz(top.z) {
				}
			};

			// Returns true if the pixel at the center of the given circular window is
			// the maximum value in the window.
			bool isMaxCenter(Grid<float> &raster, int col, int row, int window, double *max) {
				*max = 0;
				int mc = 0, mr = 0;
				float v;
				float d = g_sq((float) window / 2);
				for (int r = row - window / 2; r < row + window / 2; ++r) {
					for (int c = col - window / 2; c < col + window / 2; ++c) {
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
				float d = g_sq((float) window / 2);
				for (int r = row - window / 2; r < row + window / 2; ++r) {
					for (int c = col - window / 2; c < col + window / 2; ++c) {
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
				float d = g_sq((float) window / 2);
				for (int r = row - window / 2; r < row + window / 2; ++r) {
					for (int c = col - window / 2; c < col + window / 2; ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && c != col && r != row)
							raster.set(c, r, 0);
					}
				}
			}

			void savePoints(std::list<Top*> &tops, SQLite &db) {
				std::vector<Point*> points(tops.size());
				for (const Top *t : tops) {
					std::map<std::string, std::string> fields;
					fields["id"] = std::to_string(t->id);
					points.push_back(new Point(t->x, t->y, t->uz, fields));
				}
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
	topsMinHeight(0.0),
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
	if(topsMinHeight >= topsThresholds.begin()->first)
		g_argerr("The minimum height must be lower than the lowest threshold.");
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
	std::vector<std::string> p;
	for(const auto &it : topsThresholds) {
		p.push_back(std::to_string(it.first));
		p.push_back(std::to_string(it.second));
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

	MemRaster<uint8_t> topsGrid(original.cols(), original.rows());
	topsGrid.fill(0);

	for(const auto &it : config.topsThresholds) {

		uint8_t window = it.second;
		double threshold = it.first;

		for (int32_t row = window / 2; row < smoothed.rows() - window / 2; ++row) {
			if (*cancel)
				break;
			for (int32_t col = window / 2; col < smoothed.cols() - window / 2; ++col) {

				double v = smoothed.get(col, row);
				if(v < config.topsMinHeight || v > threshold)
					continue;

				double max;

				if (isMaxCenter(smoothed, col, row, window, &max)) {

					topsGrid.set(col, row, window);
					zeroKernel(topsGrid, col, row, window);

				}
			}
		}
	}

	std::list<Top*> tops;
	uint64_t topId = 0;

	for(int32_t row = 0; row < topsGrid.rows(); ++row) {
		for(int32_t col = 0; col < topsGrid.cols(); ++col) {

			uint8_t window = topsGrid.get(col, row);
			if(!window)
				continue;

			// Get the original height from the unsmoothed raster.
			double smax, omax;
			uint16_t mc, mr;
			getKernelMax(smoothed, col, row, window, &smax, &mc, &mr);
			getKernelMax(original, col, row, window, &omax, &mc, &mr);

			Top * pt = new Top(++topId,
				original.toCentroidX(mc), // center of pixel
				original.toCentroidY(mr),
				omax, smax, mc, mr
			);
			tops.push_back(pt);
		}

		if (tops.size() >= 1000) {

			if (m_callbacks)
				m_callbacks->statusCallback("Inserting points...");

			savePoints(tops, db);
			tops.clear();

			if (m_callbacks)
				m_callbacks->statusCallback("Processing...");
		}

	}


	if (m_callbacks)
		m_callbacks->statusCallback("Inserting points...");

	savePoints(tops, db);
	tops.clear();

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
		for (int32_t row0 = 0; row0 < inrast.rows(); row0 += rowStep) {
			if (*cancel)
				continue;
			curRow += rowStep;

			// Load the tree tops for the strip.
			Bounds bounds(inrast.toX(0), inrast.toY(row0 - bufRows),
					inrast.toX(inrast.cols()),
					inrast.toY(row0 + rowStep + bufRows));
			std::vector<std::unique_ptr<Point> > tops;

			#pragma omp critical(__crowns_getpoints)
			db.getPoints(tops, bounds);

			#pragma omp critical(__crowns_readbuf)
			inrast.readBlock(0, row0 == 0 ? row0 : row0 - bufRows, buf, 0, row0 == 0 ? bufRows : 0);

			// Convert the Tops to Nodes.
			std::queue<std::unique_ptr<Node> > q;
			for (const std::unique_ptr<Point> &t : tops) {
				if (*cancel)
					break;
				int col = inrast.toCol(t->x);
				int row = inrast.toRow(t->y);
				uint64_t id = atoi(t->fields["id"].c_str());
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

				blk.set(n->c, n->r - row0 + bufRows, (uint32_t) n->id);

				for(size_t i = 0; i < 8; ++i) {
					int c = n->c + offsets[i][0];
					int r = n->r + offsets[i][1];

					if (r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols())
						continue;
					if (r - row0 + bufRows < 0 || r - row0 + bufRows >= buf.rows())
						continue;

					uint32_t idx = (uint64_t) (r - row0 + bufRows) * inrast.cols() + c;
					if (visited[idx])
						continue;

					double v = buf.get(idx);
					if (v != nodata           												// is not nodata
						&& v < n->z		      												// is less than the neighbouring pixel
						&& v >= config.crownsMinHeight 										// is greater than the min height
						&& (n->tz - v) / n->tz <= config.crownsHeightFraction 				// is greater than the threshold height
						&& g_sq(n->tc - c) + g_sq(n->tr - r) <= g_sq(config.crownsRadius) 	// is within the radius
					) {
						std::unique_ptr<Node> nd(new Node(n->id, c, r, v, n->tc, n->tr, n->tz));
						q.push(std::move(nd));
						visited[idx] = true;
					}
				}
			}

			if (row0 > 0 && (row0 + bufRows) >= inrast.rows())
				continue;

			if (*cancel)
				continue;

			if(m_callbacks)
				m_callbacks->statusCallback("Writing output...");

			#pragma omp critical(__b)
			outrast.writeBlock(0, row0, blk, 0, bufRows);

		}
	}

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


