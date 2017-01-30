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
			// the maximum value in the window. Assigns the max pixel value to max,
			// and the proportion of nulls [0-1] to nulls.
			bool isMaxCenter(Grid &raster, int col, int row, int window, double *max, double *nulls, int band) {
				*max = 0;
				int mc = 0, mr = 0;
				int n = 0, ntotal = 0;
				double nodata = raster.props().nodata();
				double distance = g_max(std::sqrt(2.0), g_sq((double) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						double dist0 = g_sq((double) col - c) + g_sq((double) row - r);
						if (dist0 <= distance) {
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

			// Returns the max value of pixels in the kernel. Sets the
			// col and row of the pixel to mr and mc, and the max to max.
			void getKernelMax(Grid &raster, int col, int row, int window,
				double *max, int *mc, int *mr, int band) {
				*max = G_DBL_MAX_NEG;
				double distance = g_max(std::sqrt(2.0), g_sq((double) window / 2));
				for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
					for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
						double dist0 = g_sq((double) col - c) + g_sq((double) row - r);
						double v = raster.getFloat(c, r, band);
						if (dist0 <= distance && v > *max) {
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
				double d = g_max(std::sqrt(2.0), g_sq((double) window / 2));
				for (int r = g_max(0, row - window / 2); r < g_min(rows, row + window / 2 + 1); ++r) {
					for (int c = g_max(0, col - window / 2); c < g_min(cols, col + window / 2 + 1); ++c) {
						double d0 = g_sq((double) col - c) + g_sq((double) row - r);
						if (d0 <= d && c != col && r != row)
							raster.setInt(c, r, 0, band);
					}
				}
			}

			// For each ID represented in the crowns raster, finds the highest pixel value
			// in chm within the pixels corresponding to that ID. Produces
			// a map relating the ID to a tuple containing the 3D coordinate of the
			// highest pixel.
			void findCrownMax(Raster &chm, Raster &crowns,
					std::unordered_map<uint32_t, std::tuple<double, double, double> > &heights,
					int chmBand, int crownsBand, Status &status, bool *cancel) {

				const GridProps &chmProps = chm.props();
				const GridProps &crownProps = crowns.props();

				int rows = crownProps.rows();
				for(long i = 0; i < crownProps.size(); ++i) {
					uint32_t id = crowns.getInt(i, crownsBand);
					double v = chm.getFloat(i, chmBand);
					if(heights.find(id) == heights.end() || v > std::get<2>(heights[id])) {
						double x = chmProps.toCentroidX(i % chmProps.cols());
						double y = chmProps.toCentroidY(i / chmProps.rows());
						heights[id] = std::make_tuple(x, y, v);
					}
					if(i % rows == 0) {
						if(*cancel)
							break;
						status.update((float) i / crownProps.size());
					}
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
		int lastWindow = 0;
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
	return Util::join(p, ",");
}

void TreetopsConfig::parseThresholds(const std::string &str) {
	std::stringstream ss(str);
	std::string item1, item2;
	topsThresholds.clear();
	while(std::getline(ss, item1, ',')) {
		if(!std::getline(ss, item2, ','))
			break;
		if(item1.empty() || item2.empty()) continue;
		double height = atof(item1.c_str());
		int window = atoi(item2.c_str());
		if(window == 0) continue;
		topsThresholds.push_back(std::make_pair(height, window));
	}
}

// Top implementation

Top::Top(uint64_t id, uint64_t parentID, double ox, double oy, double oz, double sx, double sy, double sz, int sc, int sr) :
	id(id), parentID(parentID), ox(ox), oy(oy), oz(oz), sx(sx), sy(sy), sz(sz), sc(sc), sr(sr) {
}

Top::Top() :
	id(0), parentID(0), ox(0), oy(0), oz(0), sx(0), sy(0), sz(0), sc(0), sr(0) {
}


// TTDB implementation

TTDB::TTDB(const std::string &file, const std::string &layer, const std::string &driver,
    const std::unordered_map<std::string, FieldType> &fields, GeomType type,
	int srid, bool replace) :
    DB(file, layer, driver, fields, type, srid, replace) {}

TTDB::TTDB(const std::string &file, const std::string &layer) :
    DB(file, layer) {}

void TTDB::addTop(const std::unique_ptr<Top> &top) {
    if(m_type != GeomType::GTPoint)
        g_runerr("This dataset is not a point dataset.");
    OGRFeature* feat = OGRFeature::CreateFeature(m_fdef);
    feat->SetField("id", (GIntBig) top->id);
    feat->SetField("parentId", (GIntBig) top->parentID);
    feat->SetField("originalX", top->ox);
    feat->SetField("originalY", top->oy);
    feat->SetField("originalZ", top->oz);
    feat->SetField("smoothedX", top->sx);
    feat->SetField("smoothedY", top->sy);
    feat->SetField("smoothedZ", top->sz);
    feat->SetField("smoothedCol", top->sc);
    feat->SetField("smoothedRow", top->sr);
    OGRPoint geom(top->sx, top->sy, top->sz);
    feat->SetGeometry(&geom);
    OGRErr e = m_layer->CreateFeature(feat);
    if(CPLE_None != e)
        g_runerr("Failed to add feature to " << m_file << ".");
    OGRFeature::DestroyFeature(feat);
}

void TTDB::addTops(const std::list<std::unique_ptr<Top> > &tops) {
    for(const std::unique_ptr<Top> &t : tops)
        addTop(t);
}

void TTDB::getTops(std::vector<std::unique_ptr<Top> > &tops, const geotools::util::Bounds &bounds) {
    if(m_type != GeomType::GTPoint)
        g_runerr("This dataset is not a point dataset.");
    m_layer->SetSpatialFilterRect(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy());
    OGRFeature *feat;
    while((feat = m_layer->GetNextFeature())) {
        std::unique_ptr<Top> t(new Top(
            feat->GetFieldAsInteger64("id"),
            feat->GetFieldAsInteger64("parentID"),
            feat->GetFieldAsDouble("originalX"),
            feat->GetFieldAsDouble("originalY"),
            feat->GetFieldAsDouble("originalZ"),
            feat->GetFieldAsDouble("smoothedX"),
            feat->GetFieldAsDouble("smoothedY"),
            feat->GetFieldAsDouble("smoothedZ"),
            feat->GetFieldAsInteger("smoothedCol"),
            feat->GetFieldAsInteger("smoothedRow")
        ));
        tops.push_back(std::move(t));
        OGRFeature::DestroyFeature(feat);
    }
    m_layer->SetSpatialFilter(NULL);
}

void TTDB::updateTop(const std::unique_ptr<Top> &top) {
    std::stringstream ss;
    ss << "\"id\"=" << top->id;
    m_layer->SetAttributeFilter(ss.str().c_str());
    OGRFeature *feat = m_layer->GetNextFeature();
    if(!feat)
        g_runerr("Failed to find feature with ID: id=" << top->id);
    feat->SetField("parentID", (GIntBig) top->parentID);
    feat->SetField("originalX", top->ox);
    feat->SetField("originalY", top->oy);
    feat->SetField("originalZ", top->oy);
    feat->SetField("smoothedX", top->sx);
    feat->SetField("smoothedY", top->sy);
    feat->SetField("smoothedZ", top->sz);
    feat->SetField("smoothedCol", top->sc);
    feat->SetField("smoothedRow", top->sr);
    if(CPLE_None != m_layer->SetFeature(feat))
        g_runerr("Failed to save feature: " << top->id << ".");
    OGRFeature::DestroyFeature(feat);
}

void TTDB::updateTops(std::vector<std::unique_ptr<Top> > &tops) {
    for(const std::unique_ptr<Top> &top : tops)
        updateTop(top);
}



// Treetops implementation

void Treetops::setCallbacks(Callbacks *callbacks) {
	m_callbacks = callbacks;
}

void Treetops::smooth(const TreetopsConfig &config, bool *cancel) {
	config.checkSmoothing();

	if (!cancel)
		cancel = &__tt_cancel;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.01f);
		m_callbacks->statusCallback("Smoothing...");
	}

	Raster in(config.smoothOriginalCHM);

	GridProps pr = GridProps(in.props());
	pr.setWritable(true);
	pr.setDriver(config.smoothSmoothedCHMDriver);

	Raster out(config.smoothSmoothedCHM, pr);

	in.smooth(out, config.smoothSigma, config.smoothWindowSize,
			config.smoothOriginalCHMBand, m_callbacks, cancel);
	out.flush();
}

// TODO: Raster band selection.
void Treetops::treetops(const TreetopsConfig &config, bool *cancel) {
	config.checkTops();

	if(!cancel)
		cancel = &__tt_cancel;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.01f);
		m_callbacks->statusCallback("Treetops: Preparing...");
	}

	// Initialize input rasters.
	Raster smoothed(config.topsSmoothedCHM);

	if (*cancel)
		return;

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
	
	TTDB db(config.topsTreetopsDatabase, "data", config.topsTreetopsDatabaseDriver,
			fields, GeomType::GTPoint, config.srid, true);
	//db.dropGeomIndex();
	//db.setCacheSize(config.tableCacheSize);

	if (*cancel)
		return;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.02f);
		m_callbacks->statusCallback("Treetops: Starting...");
	}

	GridProps pr = GridProps(smoothed.props());

	// This raster tracks the size of the window
	// that was used to identify a top.
	pr.setDataType(DataType::Byte);
	MemRaster topsWindowGrid(pr, true);
	topsWindowGrid.fillInt(0);

	// This raster tracks the top IDs.
	pr.setDataType(DataType::UInt32);
	MemRaster topsIDGrid(pr, true);
	topsIDGrid.fillInt(0);

	uint64_t topId = 0;
	int bufSize = 256; // TODO: Configurable or dynamic buffer size.

	// The total number of actions to perform, and a status
	// variable to keep track of them. (For status updates).
	uint64_t total = config.topsThresholds.size() * smoothed.props().rows();
	std::atomic<uint64_t> status(0);

	// Iterate over the thresholds from lowest height to highest.
	for(const auto &it : config.topsThresholds) {

		int window = it.second;
		double threshold = it.first;

		#pragma omp parallel
		{

			// Prepare a buffer for the smoothing.
			GridProps pr = GridProps(smoothed.props());
			pr.setDataType(DataType::Float64);
			MemRaster smooth(pr, true);

			#pragma omp for
			for(int j = 0; j < smoothed.props().rows() / bufSize + 1; ++j) {

				if (*cancel)
					continue;

				// The offset of the current buffer, in rows.
				int b = j * bufSize;

				if(m_callbacks)
					m_callbacks->statusCallback("Treetops: Reading buffer...");

				int readOffset = b > 0 ? b - window / 2 : 0;  // If this is the first row, read from zero, otherwise -(size / 2)
				int writeOffset = b > 0 ? 0 : window / 2;     // If this is the first row, write to (size / 2), otherwise 0.
				smooth.fillFloat(smooth.props().nodata());
				#pragma omp critical(__tops_read)
				smoothed.write(smooth, smooth.props().cols(), smooth.props().rows(),
						0, readOffset, 0, writeOffset);

				std::list<std::tuple<int, int, int> > tops;

				if(m_callbacks)
					m_callbacks->statusCallback("Treetops: Finding tops...");

				// Iterate over the raster, applying the kernel to find the maxium at the centre.
				for (int row = window / 2; row < g_min(bufSize, smoothed.props().rows() - b) - window / 2; ++row) {
					for (int col = window / 2; col < smooth.props().cols() - window / 2; ++col) {

						double v = smooth.getFloat(col, row);
						if(v < threshold) continue;

						double max;
						double nulls; // The proportion of null pixels.
						bool isMax = isMaxCenter(smooth, col, row, window, &max, &nulls, 1);
						if (isMax && nulls <= config.topsMaxNulls)
							tops.push_back(std::make_tuple(col, b + row - window / 2, window));
					}
				}

				if(m_callbacks)
					m_callbacks->statusCallback("Treetops: Writing...");

				// Write the tops to the tops raster.
				#pragma omp critical(__tops_write)
				{
					for (const auto &top : tops) {
						topsWindowGrid.setInt(std::get<0>(top), std::get<1>(top), std::get<2>(top));
						topsIDGrid.setInt(std::get<0>(top), std::get<1>(top), ++topId);
					}
				}

				if (m_callbacks)
					m_callbacks->stepCallback(0.02f + (float) ++status / total * 0.31f);
			}
		}
	}


	if(m_callbacks)
		m_callbacks->statusCallback("Treetops: Finding parent tops...");

	// A new buffer to keep track of parent top IDs.
	pr = GridProps(smoothed.props());
	MemRaster topsParentGrid(pr, true);
	topsParentGrid.fillInt(0);

	// The smallest window used in this job.
	int minWindow = config.topsThresholds.begin()->second;

	for(int row = 0; row < pr.rows(); ++row) {

		if (*cancel)
			return;

		for(int col = 0; col < pr.cols(); ++col) {

			int window = topsWindowGrid.getInt(col, row);
			if(window <= minWindow) continue; 			// If the top was found by the smallest window, skip it.

			// Find the tops which are inside of the window of another top's window.
			for(int r = g_max(0, row - window / 2); r < g_min(pr.rows(), row + window / 2 + 1); ++r) {
				for(int c = g_max(0, col - window / 2); c < g_min(pr.cols(), col + window / 2 + 1); ++c) {
					if(c == col || r == row || (g_sq(c - col) + g_sq(r - row)) > g_sq(window)) continue;

					int window0 = topsWindowGrid.getInt(c, r);
					if(window0 < window && !topsParentGrid.getInt(c, r))
						topsParentGrid.setInt(c, r, topsIDGrid.getInt(col, row)); //  Save the ID of the parent.
				}
			}
		}
		if (m_callbacks)
			m_callbacks->stepCallback(0.33f + (float) row / pr.rows() * 0.33f);
	}

	// Finally, scrape up all the tops and put them in the DB.
	if(m_callbacks)
		m_callbacks->statusCallback("Treetops: Saving tops...");

	total = topsIDGrid.props().rows();
	status = 0;

	#pragma omp parallel
	{

		std::list<std::unique_ptr<Top> > tops;

		#pragma omp for
		for(int row = 0; row < topsIDGrid.props().rows(); ++row) {

			if(*cancel)
				continue;

			for(int col = 0; col < topsIDGrid.props().cols(); ++col) {

				int window = topsWindowGrid.getInt(col, row);
				if(!window) continue;

				double z;
				#pragma omp critical(__smoothed_read)
				z = smoothed.getFloat(col, row);

				std::unique_ptr<Top> t(new Top(topsIDGrid.getInt(col, row),
					topsParentGrid.getInt(col, row),
					0, // original.toCentroidX(col),
					0, // original.toCentroidY(row),
					0, // original.get(col, row),
					pr.toCentroidX(col),
					pr.toCentroidY(row),
					z,
					col, row
				));
				tops.push_back(std::move(t));
			}

			++status;

			if (tops.size() >= 5000) {
				#pragma omp critical(__save_points)
				{
					db.begin();
					db.addTops(tops);
					db.commit();
					tops.clear();
				}
				if (m_callbacks)
					m_callbacks->stepCallback(0.66f + (float) status / total * 0.32f);
			}
		}

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
			m_callbacks->statusCallback("Treetops: Building index...");
		db.createGeomIndex();
	}

	if (m_callbacks) {
		m_callbacks->statusCallback("Treetops: Done.");
		m_callbacks->stepCallback(1.0);
	}

}

void Treetops::treecrowns(const TreetopsConfig &config, bool *cancel) {
	config.checkCrowns();

	if (!cancel)
		cancel = &__tt_cancel;

	if (m_callbacks) {
		m_callbacks->stepCallback(0.01f);
		m_callbacks->statusCallback("Crowns: Preparing...");
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

	if (m_callbacks)
		m_callbacks->stepCallback(0.02f);

	// Initialize the database, get the treetop count and estimate the buffer size.
	TTDB db(config.crownsTreetopsDatabase, "data");

	if (m_callbacks) {
		m_callbacks->stepCallback(0.03f);
		m_callbacks->statusCallback("Crowns: Processing...");
	}

	// Buffer for reading/processing
	// TODO: Configurable/dynamic buffer size.
	int bufSize = 256;
	int16_t radius = (int16_t) std::ceil(config.crownsRadius / g_abs(inrast.props().resolutionX()));

	{
		// Build the list of offsets for D8 search.
		size_t offsetCount = 8;
		int offsets[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

		// Values for status indicator.
		uint64_t geomCount = db.getGeomCount();
		std::atomic<uint64_t> status(0);

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

				if (*cancel)
					continue;

				int b = i * bufSize;

				std::vector<std::unique_ptr<Top> > tops;

				Bounds bounds(iprops.toX(0), iprops.toY(b - radius),
						iprops.toX(iprops.cols()), iprops.toY(b + bufSize + radius));

				if (m_callbacks)
					m_callbacks->statusCallback("Crowns: Loading tops...");

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
					m_callbacks->statusCallback("Crowns: Loading raster...");

				std::vector<bool> visited((size_t) iprops.cols() * (bufSize + radius * 2 + 1));

				int readOffset = b > 0 ? b - radius : 0;
				int writeOffset = b > 0 ? 0 : radius;
				buf.fillFloat(iprops.nodata());
				blk.fillInt(0);

				#pragma omp critical(__crowns_in)
				inrast.write(buf, buf.props().cols(), buf.props().rows(),
						0, readOffset, 0, writeOffset);

				if (m_callbacks)
					m_callbacks->statusCallback("Crowns: Delineating crowns...");

				// Run through the queue.
				while (!*cancel && q.size()) {
					std::unique_ptr<Node> n = std::move(q.front());
					q.pop();

					int c = n->c;
					int r = n->r - b + radius;

					blk.setInt(c, r, (uint) n->id);
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
				if(m_callbacks)
					m_callbacks->statusCallback("Crowns: Writing output...");

				#pragma omp critical(__crowns_out)
				blk.write(outrast, iprops.cols(), bufSize, 0, radius, 0, b);

				if(m_callbacks)
					m_callbacks->stepCallback(0.03f + (float) status / geomCount * 0.30f);
			}
		}
	}

	{
		if(m_callbacks)
			m_callbacks->statusCallback("Crowns: Finding original top heights...");

		if (*cancel)
			return;

		std::unordered_map<uint, std::tuple<double, double, double> > heights;
		Raster chm(config.crownsOriginalCHM);
		Status status(m_callbacks, 0.33f, 0.50f);
		findCrownMax(chm, outrast, heights, 1, 1, status, cancel);

		int steps = iprops.rows() / bufSize + 1;
		for(int i = 0; i < steps; ++i) {

			if (*cancel)
				continue;

			int b = i * bufSize;
			Bounds bounds(iprops.toX(0), iprops.toY(b - radius),
					iprops.toX(iprops.cols()), iprops.toY(b + bufSize + radius));

			std::vector<std::unique_ptr<Top> > tops;
			db.getTops(tops, bounds);
			if(tops.empty())
				continue;

			auto end = heights.end();
			for(const std::unique_ptr<Top> &t : tops) {
				if(heights.find(t->id) != end) {
					auto tup = heights[t->id];
					t->ox = std::get<0>(tup);
					t->oy = std::get<1>(tup);
					t->oz = std::get<2>(tup);
				}
			}

			db.begin();
			db.updateTops(tops);
			db.commit();

			if(m_callbacks)
				m_callbacks->stepCallback(0.50f + (float) i / steps * 0.25f);
		}
	}

	if(!config.crownsCrownsDatabase.empty()) {
		if(m_callbacks)
			m_callbacks->statusCallback("Polygonizing...");
		Status status(m_callbacks, 0.75f, 0.99f);
		outrast.polygonize(config.crownsCrownsDatabase, "crowns", db.srid(), 1, &status, cancel);
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


