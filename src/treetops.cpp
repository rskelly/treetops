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
#include "util.hpp"
#include "raster.hpp"
#include "treetops.hpp"

using namespace geotools::raster;
using namespace geotools::util;

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
				int c, r;
				double z;
				int tc, tr;
				double tz;

				Node(uint64_t id, int c, int r, double z, int tc, int tr, double tz) :
						id(id),
						c(c), r(r), z(z),
						tc(tc), tr(tr), tz(tz) {
				}

				Node(Top *p) :
					id(p->id), 
					c(p->sc), r(p->sr), z(p->sz),
					tc(p->sc), tr(p->sr), tz(p->sz) {
				}
			};

			// Returns true if the pixel at the center of the given circular window is
			// the maximum value in the window.
			bool isMaxCenter(Grid &raster, int col, int row, int window, double *max, double *nulls, int band) {
				*max = 0;
				int mc = 0, mr = 0;
				int n = 0, ntotal = 0;
				double nodata = raster.props().nodata();
				double d = g_max(std::sqrt(2.0), g_sq((double) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						double d0 = g_sq((double) col - c) + g_sq((double) row - r);
						if (d0 <= d) {
							++ntotal;
							double v = raster.getFloat(c, r, band);
							if(v == nodata) {
								++n;
							} else if(v > *max) {
								*max = v;
								mc = c;
								mr = r;
							}
						}
					}
				}
				*nulls = (double) n / ntotal;
				return mc == col && mr == row;
			}

			// Returns the max value of pixels in the kernel.
			void getKernelMax(Grid &raster, int col, int row, int window,
				double *max, uint16_t *mc, uint16_t *mr, int band) {
				*max = G_DBL_MAX_NEG;
				float v;
				float d = g_max(std::sqrt(2.0), g_sq((float) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && (v = raster.getFloat(c, r, band)) > *max) {
							*max = v;
							*mc = c;
							*mr = r;
						}
					}
				}
			}

			// Set all cells in the circular kernel to zero except the center.
			void zeroKernel(Grid &raster, int col, int row, int window, int band) {
				int cols = raster.props().cols();
				int rows = raster.props().rows();
				float d = g_max(std::sqrt(2.0), g_sq((float) window / 2));
				for (int r = g_max(0, row - window / 2); r < g_min(rows, row + window / 2 + 1); ++r) {
					for (int c = g_max(0, col - window / 2); c < g_min(cols, col + window / 2 + 1); ++c) {
						float d0 = g_sq((float) col - c) + g_sq((float) row - r);
						if (d0 <= d && c != col && r != row)
							raster.setInt(c, r, 0, band);
					}
				}
			}

			void findCrownMax(Raster &chm, Raster &crowns,
					std::unordered_map<uint32_t, std::tuple<double, double, double> > &heights, int chmBand, int crownsBand) {
				const GridProps &chmProps = chm.props();
				for(long i = 0; i < crowns.props().size(); ++i) {
					uint32_t id = crowns.getInt(i, crownsBand);
					double v = chm.getFloat(i);
					if(heights.find(id) == heights.end() || v > std::get<2>(heights[id]))
						heights[id] = std::make_tuple(chmProps.toCentroidX(i % chmProps.cols()), chmProps.toCentroidY(i / chmProps.rows()), v);
				}
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
	smoothOriginalCHMBand(1),
	doTops(false),
	topsMaxNulls(0.2),
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
	if (smoothSmoothedCHMDriver.empty())
		g_argerr("Smoothing: Output driver must not be empty.");
	if (smoothSigma <= 0 || smoothSigma > 1)
		g_argerr("Smoothing: Std. deviation must be 0 < n <= 1. " << smoothSigma << " given.");
	if (smoothWindowSize % 2 == 0 || smoothWindowSize < 3)
		g_argerr("Smoothing: The window must be odd and >=3.");
}

void TreetopsConfig::checkTops() const {
	if (!doTops)
		g_argerr("Not configured to find treetops.");
	if (topsSmoothedCHM.empty())
		g_argerr("Tops: Smoothed CHM filename must not be empty.");
	if (topsTreetopsDatabase.empty())
		g_argerr("Tops: Database filename must not be empty.");
	if (topsTreetopsDatabaseDriver.empty())
		g_argerr("Tops: Database driver must not be empty.");
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
	if (crownsCrownsRasterDriver.empty())
		g_argerr("Crowns: Output raster driver must not be empty.")
	if(!crownsCrownsDatabase.empty() && crownsCrownsDatabaseDriver.empty())
		g_argerr("Crowns: If database file is given, driver must also be given.");
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
		topsThresholds.push_back(std::make_pair(height, window));
	}
}

// Top implementation

Top::Top(uint64_t id, uint64_t parentID, double ox, double oy, double oz, double sx, double sy, double sz, int32_t sc, int32_t sr) :
	id(id), parentID(parentID), ox(ox), oy(oy), oz(oz), sx(sx), sy(sy), sz(sz), sc(sc), sr(sr) {
}

Top::Top() :
	id(0), parentID(0), ox(0), oy(0), oz(0), sx(0), sy(0), sz(0), sc(0), sr(0) {
}


// Treetops implementation

void Treetops::setCallbacks(Callbacks *callbacks) {
	m_callbacks = callbacks;
}

void Treetops::smooth(const TreetopsConfig &config, bool *cancel) {
	config.checkSmoothing();

	if (!cancel)
		cancel = &__tt_cancel;

	Raster in(config.smoothOriginalCHM);
	GridProps pr = GridProps(in.props());
	pr.setWritable(true);
	pr.setDriver(config.smoothSmoothedCHMDriver);
	Raster out(config.smoothSmoothedCHM, pr);
	in.smooth(out, config.smoothSigma, config.smoothWindowSize,
			config.smoothOriginalCHMBand, m_callbacks, cancel);
	out.flush();
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
	Raster smoothed(config.topsSmoothedCHM);

	if (*cancel)
		return;

	// Prepare database.
	if(m_callbacks)
		m_callbacks->statusCallback("Preparing database...");

	std::unordered_map<std::string, FieldType> fields;
	fields["id"] = FieldType::FTInt;
	fields["parentID"] = FieldType::FTInt;
	fields["originalX"] = FieldType::FTDouble;
	fields["originalY"] = FieldType::FTDouble;
	fields["originalZ"] = FieldType::FTDouble;
	fields["smoothedX"] = FieldType::FTDouble;
	fields["smoothedY"] = FieldType::FTDouble;
	fields["smoothedZ"] = FieldType::FTDouble;
	fields["smoothedCol"] = FieldType::FTInt;
	fields["smoothedRow"] = FieldType::FTInt;
	
	TTDB db(config.topsTreetopsDatabase, "data", fields, GeomType::GTPoint, config.srid, true);
	db.dropGeomIndex();
	db.setCacheSize(config.tableCacheSize);

	if (*cancel)
		return;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Processing...");
	}

	GridProps pr = GridProps(smoothed.props());
	pr.setDataType(DataType::Byte);
	MemRaster topsWindowGrid(pr, true);
	topsWindowGrid.fillInt(0);

	pr.setDataType(DataType::UInt32);
	MemRaster topsIDGrid(pr, true);
	topsIDGrid.fillInt(0);

	uint64_t topId = 0;

	uint64_t total = config.topsThresholds.size() * topsWindowGrid.props().rows();
	std::atomic<uint64_t> status(0);

	uint16_t bufSize = 256;

	omp_set_num_threads(1);

	// Iterate over the thresholds from lowest height to highest.
	for(const auto &it : config.topsThresholds) {

		uint8_t window = it.second;
		double threshold = it.first;

		#pragma omp parallel
		{

			GridProps pr = GridProps(smoothed.props());
			pr.setDataType(DataType::Float64);
			MemRaster smooth(pr, true);

			#pragma omp for
			for(uint16_t j = 0; j < smoothed.props().rows() / bufSize + 1; ++j) {
				if (*cancel) continue;

				int32_t b = j * bufSize;

				if(m_callbacks)
					m_callbacks->statusCallback("Reading...");

				uint16_t readOffset = b > 0 ? b - window / 2 : 0;  // If this is the first row, read from zero, otherwise -(size / 2)
				uint16_t writeOffset = b > 0 ? 0 : window / 2;     // If this is the first row, write to (size / 2), otherwise 0.
				smooth.fillFloat(smooth.props().nodata());
				#pragma omp critical(__tops_read)
				smoothed.writeToBlock(smooth, smooth.props().cols(), smooth.props().rows(),
						0, readOffset, 0, writeOffset);

				std::list<std::tuple<int32_t, int32_t, uint8_t> > tops;

				if(m_callbacks)
					m_callbacks->statusCallback("Processing...");

				// Iterate over the raster, applying the kernel to find the maxium at the centre.
				for (int32_t row = window / 2; row < g_min(bufSize, smoothed.props().rows() - b) - window / 2; ++row) {
					for (int32_t col = window / 2; col < smooth.props().cols() - window / 2; ++col) {

						double v = smooth.getFloat(col, row);
						if(v < threshold) continue;

						double max;
						double nulls; // The proportion of null pixels.
						bool isMax = isMaxCenter(smooth, col, row, window, &max, &nulls, 1);
						if (isMax && nulls <= config.topsMaxNulls)
							tops.push_back(std::make_tuple(col, b + row - window / 2, window));
					}
					if (m_callbacks)
						m_callbacks->stepCallback(0.02f + (float) ++status / total * 0.48f);
				}

				if(m_callbacks)
					m_callbacks->statusCallback("Writing...");

				// Write the tops to the tops raster.
				#pragma omp critical(__tops_write)
				{
					for (const auto &top : tops) {
						topsWindowGrid.setInt(std::get<0>(top), std::get<1>(top), std::get<2>(top));
						topsIDGrid.setInt(std::get<0>(top), std::get<1>(top), ++topId);
					}
				}
			}
		}
	}

	pr = GridProps(smoothed.props());
	MemRaster topsParentGrid(pr, true);
	topsParentGrid.fillInt(0);

	uint8_t minWindow = config.topsThresholds.begin()->second;

	for(uint16_t row = 0; row < pr.rows(); ++row) {
		for(uint16_t col = 0; col < pr.cols(); ++col) {

			uint8_t window = topsWindowGrid.getInt(col, row);
			if(window <= minWindow) continue;

			for(int32_t r = g_max(0, row - window / 2); r < g_min(pr.rows(), row + window / 2 + 1); ++r) {
				for(int32_t c = g_max(0, col - window / 2); c < g_min(pr.cols(), col + window / 2 + 1); ++c) {
					if(c == col || r == row || (g_sq(c - col) + g_sq(r - row)) > g_sq(window)) continue;

					uint8_t window0 = topsWindowGrid.getInt(c, r);
					if(window0 < window && !topsParentGrid.getInt(c, r))
						topsParentGrid.setInt(c, r, topsIDGrid.getInt(col, row)); //  Save the ID of the parent.
				}
			}
		}
	}

	total = topsIDGrid.props().rows();
	status = 0;

	#pragma omp parallel
	{

		std::list<std::unique_ptr<Top> > tops;

		#pragma omp for
		for(int32_t row = 0; row < topsIDGrid.props().rows(); ++row) {
			if(*cancel) continue;
			if (m_callbacks)
				m_callbacks->statusCallback("Processing...");
			for(int32_t col = 0; col < topsIDGrid.props().cols(); ++col) {

				uint8_t window = topsWindowGrid.getInt(col, row);
				if(!window) continue;

				std::unique_ptr<Top> t(new Top(topsIDGrid.getInt(col, row),
					topsParentGrid.getInt(col, row),
					0, // original.toCentroidX(col),
					0, // original.toCentroidY(row),
					0, // original.get(col, row),
					pr.toCentroidX(col),
					pr.toCentroidY(row),
					smoothed.getInt(col, row),
					col, row
				));
				tops.push_back(std::move(t));
			}

			if (m_callbacks)
				m_callbacks->stepCallback(0.5f + (float) ++status / total * 0.48f);

			if (tops.size() >= 5000) {
				if (m_callbacks)
					m_callbacks->statusCallback("Inserting points...");
				#pragma omp critical(__save_points)
				{
					db.begin();
					db.addTops(tops);
					db.commit();
					tops.clear();
				}

			}
		}
		if (m_callbacks)
			m_callbacks->statusCallback("Inserting points...");
		#pragma omp critical(__save_points)
		{
			db.begin();
			db.addTops(tops);
			db.commit();
			tops.clear();
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
		m_callbacks->statusCallback("Crowns: Loading raster...");
	}

	// Load the smoothed CHM.
	Raster inrast(config.crownsSmoothedCHM);
	const GridProps &iprops = inrast.props();
	double nodata = iprops.nodata();

	// Initialize the output raster for writing.
	GridProps pr = GridProps(iprops);
	pr.setBands(1);
	pr.setDataType(DataType::UInt32);
	pr.setNoData(0);
	Raster outrast(config.crownsCrownsRaster, pr);

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Crowns: Preparing database...");
	}

	// Initialize the database, get the treetop count and estimate the buffer size.
	TTDB db(config.crownsTreetopsDatabase, "data");

	if (m_callbacks) {
		m_callbacks->stepCallback(0.03f);
		m_callbacks->statusCallback("Crowns: Processing...");
	}

	int bufSize = 256;
	int16_t radius = (int16_t) std::ceil(config.crownsRadius / g_abs(inrast.props().resolutionX()));

	{
		// Build the list of offsets for D8 search.
		size_t offsetCount = 8;
		int offsets[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

		uint64_t geomCount = db.getGeomCount();
		std::atomic<uint64_t> status(0);

		omp_set_num_threads(1);

		#pragma omp parallel
		{

			// A buffer for reading the CHM.
			GridProps pr = GridProps(inrast.props());
			pr.setDataType(DataType::Float64);
			pr.setSize(pr.cols(), bufSize + radius * 2 + 1);
			MemRaster buf(pr);

			// A block for writing the crowns.
			pr.setDataType(DataType::UInt32);
			MemRaster blk(pr);
			blk.fillInt(0);

			#pragma omp for
			for(int i = 0; i < inrast.props().rows() / bufSize + 1; ++i) {
				if (*cancel) continue;

				int b = i * bufSize;

				std::vector<std::unique_ptr<Top> > tops;

				Bounds bounds(iprops.toX(0), iprops.toY(b - radius),
						iprops.toX(iprops.cols()), iprops.toY(b + bufSize + radius));

				#pragma omp critical(__crowns_db)
				db.getTops(tops, bounds);
				if(tops.empty())
					continue;

				status += tops.size();

				// Convert the Tops to Nodes.
				std::queue<std::unique_ptr<Node> > q;
				for (std::unique_ptr<Top> &t : tops) {
					std::unique_ptr<Node> n(new Node(t.get()));
					q.push(std::move(n));
				}
				tops.clear();

				if (*cancel)
					continue;

				if (m_callbacks)
					m_callbacks->statusCallback("Loading raster...");

				std::vector<bool> visited((size_t) iprops.cols() * (bufSize + radius * 2 + 1));

				uint16_t readOffset = b > 0 ? b - radius : 0;
				uint16_t writeOffset = b > 0 ? 0 : radius;
				buf.fillFloat(iprops.nodata());
				blk.fillInt(0);

				#pragma omp critical(__crowns_in)
				inrast.writeToBlock(buf, buf.props().cols(), buf.props().rows(),
						0, readOffset, 0, writeOffset);

				if (m_callbacks)
					m_callbacks->statusCallback("Delineating crowns...");

				// Run through the queue.
				while (!*cancel && q.size()) {
					std::unique_ptr<Node> n = std::move(q.front());
					q.pop();

					int c = n->c;
					int r = n->r - b + radius;

					blk.setInt(c, r, (uint32_t) n->id);
					for(size_t i = 0; i < offsetCount; ++i) {
						c = n->c + offsets[i][0];
						r = n->r + offsets[i][1];

						if ((r - b + radius) < 0 || c < 0 ||
								(r - b + radius) >= buf.props().rows() || c >= buf.props().cols())
							continue;

						size_t idx = (size_t) (r - b + radius) * iprops.cols() + c;
						if (visited[idx])
							continue;

						double v = buf.getFloat(idx);
						if (v != nodata           												// is not nodata
							&& v < n->z		      												// is less than the neighbouring pixel
							&& v >= config.crownsMinHeight 										// is greater than the min height
							&& (n->tz - v) / n->tz <= config.crownsHeightFraction 				// is greater than the threshold height
							&& g_sq(n->tc - c) + g_sq(n->tr - r) <= g_sq(config.crownsRadius) 	// is within the radius
						) {
							blk.setInt(idx, n->id);
							visited[idx] = true;
							std::unique_ptr<Node> n0(new Node(n->id, c, r, v, n->tc, n->tr, n->tz));
							q.push(std::move(n0));
						}
					}
				}
				if(m_callbacks) {
					m_callbacks->statusCallback("Writing output...");
					m_callbacks->stepCallback(0.03f + (float) status / geomCount * 0.95f);
				}

				#pragma omp critical(__crowns_out)
				blk.writeToBlock(outrast, iprops.cols(), bufSize, 0, radius, 0, b);
			}
		}
	}

	{
		if(m_callbacks)
			m_callbacks->statusCallback("Finding tops from original CHM...");

		std::unordered_map<uint32_t, std::tuple<double, double, double> > heights;
		Raster chm(config.crownsOriginalCHM);
		findCrownMax(chm, outrast, heights, 1, 1);

		for(int i = 0; i < iprops.rows() / bufSize + 1; ++i) {
			if (*cancel) continue;

			int b = i * bufSize;
			Bounds bounds(iprops.toX(0), iprops.toY(b - radius),
					iprops.toX(iprops.cols()), iprops.toY(b + bufSize + radius));

			if (m_callbacks)
				m_callbacks->statusCallback("Loading tops...");

			std::vector<std::unique_ptr<Top> > tops;
			db.getTops(tops, bounds);
			if(tops.empty())
				continue;

			if (m_callbacks)
				m_callbacks->statusCallback("Updating tops...");

			auto end = heights.end();
			for(const std::unique_ptr<Top> &t : tops) {
				if(heights.find(t->id) != end) {
					auto tup = heights[t->id];
					t->ox = std::get<0>(tup);
					t->oy = std::get<1>(tup);
					t->oz = std::get<2>(tup);
				}
			}

			if(m_callbacks)
				m_callbacks->statusCallback("Saving tops...");

			db.begin();
			db.updateTops(tops);
			db.commit();
		}
	}

	if(m_callbacks)
		m_callbacks->stepCallback(0.99f);

	if(!config.crownsCrownsDatabase.empty()) {
		if(m_callbacks)
			m_callbacks->statusCallback("Polygonizing...");
		outrast.polygonize(config.crownsCrownsDatabase, "crowns", db.srid(), 1, m_callbacks, cancel);
		if(m_callbacks)
			m_callbacks->statusCallback("Deleting invalid polygons...");
		TTDB db(config.crownsCrownsDatabase, "crowns");
		db.begin();
		db.deleteFeature("id", 0);
		db.commit();
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


