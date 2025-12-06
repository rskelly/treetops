/*
 *  Created on: May 3, 2016
 *      Author: Rob Skelly
 *       Email: rob@dijital.ca
 */

#include "geo.hpp"

#include <queue>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <cstdint>

#include <ogr_feature.h>
#include <ogrsf_frmts.h>

#include "ds/simple_interval_tree.hpp"
#include "ds/mqtree.hpp"
#include "crowns.hpp"
#include "grid.hpp"
#include "util.hpp"

using namespace dijital;
using namespace dijital::grid;
using namespace dijital::util;
using namespace dijital::ds;
using namespace dijital::crowns;
using namespace dijital::crowns::util;
using namespace dijital::crowns::config;


dijital::crowns::util::Node::Node(int id, int c, int r, float z, int tc, int tr, float tz) :
		id(id),
		c(c), r(r), z(z),
		tc(tc), tr(tr), tz(tz) {
}

dijital::crowns::util::Node::Node(Top& p) :
	id(p.id),
	c(p.sc), r(p.sr), z(p.sz),
	tc(p.sc), tr(p.sr), tz(p.sz) {
}

dijital::crowns::util::Top::Top() :
		Top(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) {
}

dijital::crowns::util::Top::Top(size_t id, size_t parentId,
	float ox, float oy, float oz,
	float sx, float sy, float sz,
	float groundZ,
	int sc, int sr) :
		pos(0),
		id(id), parentID(parentId),
		ox(ox), oy(oy), oz(oz),
		sx(sx), sy(sy), sz(sz),
		groundZ(groundZ),
		sc(sc), sr(sr) {
}

float dijital::crowns::util::Top::x() const {
	return sx;
}

float dijital::crowns::util::Top::y() const {
	return sy;
}

void dijital::crowns::util::Top::update(size_t id, size_t parentId,
		float ox, float oy, float oz,
		float sx, float sy, float sz,
		int sc, int sr) {
	this->id = id;
	this->parentID = parentId;
	this->ox = ox;
	this->oy = oy;
	this->oz = oz;
	this->sx = sx;
	this->sy = sy;
	this->sz = sz;
	this->sc = sc;
	this->sr = sr;
}

bool dijital::crowns::util::isMaxCenter(std::vector<float>& raster,
		int col, int row, int cols, int rows, int window, float nodata,
		float& maxOut, float& nullsOut) {

	float v, max = dijital::minvalue<float>();
	int mc = -1, mr = -1, n = 0;
	for(int r = -window / 2; r < window / 2 + 1; ++r) {
		for(int c = -window / 2; c < window / 2 + 1; ++c) {
			int cc = col + c;
			int rr = row + r;
			if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows)
					&& (v = raster[rr * cols + cc]) != nodata) {
				if(v > max) {
					max = v;
					mc = cc;
					mr = rr;
				}
			} else {
				++n;
			}
		}
	}

	// Compute the proportion of nulls; +1 for the centre cell.
	nullsOut = (float) n / (window * window);
	maxOut = max;
	return mc == col && mr == row;
}

// Returns the max value of pixels in the kernel. Sets the
// col and row of the pixel to mr and mc, and the max to max.
void dijital::crowns::util::getKernelMax(Band<float>& raster, int col, int row, int window,
	float& max, int& mc, int& mr) {
	max = dijital::minvalue<float>();
	float distance = dijital::max(std::sqrt(2.0f), dijital::sq((float) window / 2));
	for (int r = row - window / 2; r < row + window / 2 + 1; ++r) {
		for (int c = col - window / 2; c < col + window / 2 + 1; ++c) {
			float dist0 = dijital::sq((float) col - c) + dijital::sq((float) row - r);
			float v = raster.get(c, r);
			if (dist0 <= distance && v > max) {
				max = v;
				mc = c;
				mr = r;
			}
		}
	}
}

// Set all cells in the circular kernel to zero except the center.
void dijital::crowns::util::zeroKernel(Band<int>& raster, int col, int row, int window) {
	int cols = raster.props().cols();
	int rows = raster.props().rows();
	float d = dijital::max(std::sqrt(2.0f), dijital::sq((float) window / 2));
	for (int r = dijital::max(0, row - window / 2); r < dijital::min(rows, row + window / 2 + 1); ++r) {
		for (int c = dijital::max(0, col - window / 2); c < dijital::min(cols, col + window / 2 + 1); ++c) {
			float d0 = dijital::sq((float) col - c) + dijital::sq((float) row - r);
			if (d0 <= d && c != col && r != row)
				raster.set(c, r, 0);
		}
	}
}

/**
 * For each ID represented in the crowns raster, finds the highest pixel value
 * in CHM within the pixels corresponding to that ID. Produces
 * a map relating the ID to a tuple containing the 3D coordinate of the
 * highest pixel.
 */
void dijital::crowns::util::findCrownMax(Band<float>& chm, Band<uint32_t>& crowns,
		std::unordered_map<size_t, CrownConfig>& heights) {

	const GridProps& chmProps = chm.props();
	const GridProps& crownProps = crowns.props();
	int cols = crownProps.cols();
	int rows = crownProps.rows();

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			size_t id = crowns.get(col, row);
			if(id > 0) {
				float z = chm.get(col, row);
				CrownConfig& c = heights[id];
				if(z > c.z) {
					c.x = chmProps.toX(col);
					c.y = chmProps.toY(row);
					c.z = z;
				}
			}
		}
		if(Monitor::get().canceled())
			break;
		Monitor::get().status((float) row / rows, "");
	}
}

// Returns true if a pixel, represented by c, r, z is a valid crown
// pixel, given the thresholds and the location of the treetop, represented
// by n.
bool dijital::crowns::util::verifyCrownPixel(int col, int row, float z, float nodata, float res,
		const std::unique_ptr<Node>& n, CrownsAppConfig& config,
		const dijital::ds::SimpleIntervalTree<float, size_t>& st) {

	if(z == nodata || z >= n->z) // is not nodata and is not less than the neighbouring pixel
		return false;

	size_t idx;
	if(st.find(z, &idx)) {
		const CrownThreshold& ct = config.crownsThresholds()[idx];
		float radius = dijital::sq(ct.radius/ res);										// The squared radius in pixels.
		if(z >= ct.threshold 														// is greater than the min height
			&& (n->tz - z) / n->tz <= ct.fraction									// is greater than the threshold height
			&& dijital::sq((float) n->c - col) + dijital::sq((float) n->r - row) <= radius) {	// is within the radius
			return true;
		}
	}
	return false;
}

// Returns the largest radius threshold.
float dijital::crowns::util::maxCrownRadius(CrownsAppConfig& config) {

	float m = 0, m0;
	for(const CrownThreshold& ct : config.crownsThresholds()) {
		if((m0 = ct.radius) > m)
			m = m0;
	}
	return m;
}

// Locates the treetop height from the original CHM. Only searches
// within the delineated crown for the highest pixel.
void dijital::crowns::util::updateOriginalCHMHeights(CrownsAppConfig& config) {

	mqtree<Top>& qt = config.topsTree();
	CrownDB& cdb = config.crownDB();

	Monitor::get().status(0.0f, "Crowns: Finding original top heights...");

	if (Monitor::get().canceled())
		return;

	// Find the crown maximums.
	std::unordered_map<size_t, CrownConfig> heights;
	{
		// Load the original CHM.
		Band<float>& chm = config.chm();
		// Initialize the crowns raster for reading.
		Band<uint32_t>& crowns = config.crowns();
		findCrownMax(chm, crowns, heights);
	}

	Monitor::get().status(0.5f, "Crowns: Updating tops...");

	qt.reset();
	Top top;
	int stat = 0;
	int count = cdb.count();
	cdb.reset();
	while(cdb.next(top)) {
		Monitor::get().status(0.5f + (float) ++stat / count * 0.49f);
		if(heights.find(top.id) != heights.end()) {
			const CrownConfig& c = heights[top.id];
			top.ox = c.x;
			top.oy = c.y;
			top.oz = c.z;
			cdb.update(top);
		}
	}

	Monitor::get().status(1.0f, "");
}


void dijital::crowns::util::delineateCrowns(CrownsAppConfig& config) {

	Monitor::get().status(0.0f, "Crowns: Preparing...");

	// Load the smoothed CHM.
	Band<float>& smoothed = config.smoothed();
	Band<uint32_t>& crowns = config.crowns();
	mqtree<Top>& qt = config.topsTree();

	crowns.fill(0);

	// The interval tree keeps track of ranges of completed rows
	dijital::ds::SimpleIntervalTree<float, size_t> st;
	const std::vector<CrownThreshold>& thresh = config.crownsThresholds();
	float maxRadius = 0;
	for(size_t i = 0; i < thresh.size(); ++i) {
		st.add(thresh[i].fraction, i);
		if(thresh[i].radius > maxRadius)
			maxRadius = thresh[i].radius;
	}

	const GridProps& iprops = smoothed.props();
	int cols = iprops.cols();
	int rows = iprops.rows();
	float nodata = iprops.nodata();
	float resX = iprops.resX();

	Monitor::get().status(0.02f, "Crowns: Calculating tiles...");

	// Build the list of offsets for D8 search.
	size_t offsetCount = 8;
	int offsets[][2] = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};

	// Set up list of tiles for piecewise handling.
	// Tiles will have a buffer added, then removed before writing back.
	int tileSize = 512;
	std::list<std::pair<int, int>> tiles;
	for(int tr = 0; tr < rows; tr += tileSize) {
		for(int tc = 0; tc < cols; tc += tileSize)
			tiles.emplace_back(tc, tr);
	}

	// A list to track visited pixels.
	std::vector<bool> visited((size_t) cols * (size_t) rows);
	std::fill(visited.begin(), visited.end(), false);

	// Collect unique top IDs to count status.
	size_t topCount = 0;
	std::unordered_set<size_t> topSet;
	{
		Top t;
		qt.reset();
		while(qt.next(t)) {
			topSet.insert(t.id);
			++topCount;
		}
	}

	Monitor::get().status(0.03f, "Crowns: Delineating...");

	// Convert the Tops to Nodes, add to work queue.
	std::queue<Node> q;
	Bounds<float> bounds;
	Top query;
	for(const auto& tile : tiles) {
		// Create a bounding box to search for tops; this is the
		// buffer size, plus a fringe equal to max radius.
		bounds.set(
				iprops.toX(tile.first), iprops.toY(tile.second),
				iprops.toX(tile.first + tileSize), iprops.toY(tile.second + tileSize)
		);
		bounds.extendX(bounds.minx() - maxRadius);
		bounds.extendY(bounds.miny() - maxRadius);
		bounds.extendX(bounds.maxx() + maxRadius);
		bounds.extendY(bounds.maxy() + maxRadius);

		// Search radius is the diagonal of the box plus maxRadius.
		float rad = std::sqrt(std::pow(bounds.width() / 2, 2) + std::pow(bounds.height() / 2, 2)) + maxRadius;

		// Search starts at the box's center.
		query.update(0, 0, 0, 0, 0, bounds.midx(), bounds.midy(), 0, 0, 0);

		// Search for the tops.
		std::list<Top> tops;
		if(!qt.search(query, rad, std::back_inserter(tops)))
			continue;

		// Enqueue the tops for processing. TODO: Would be nice to search and enqueue in one step.
		for(Top& top : tops)
			q.emplace(top);

		// Run through the queue.
		while (!Monitor::get().canceled() && !q.empty()) {

			// Get the next node.
			Node n = std::move(q.front());
			q.pop();

			// If this is an original top, remove from the status set and update monitor.
			if(n.c == n.tc && n.r == n.tr) {
				topSet.erase(n.id);
				Monitor::get().status(1.0f - (float) topSet.size() / topCount);
			}

			// Calculate pixel index; set the ID at the current pixel in the crowns raster.
			size_t i = (size_t) n.r * (size_t) cols + (size_t) n.c;
			if(visited[i])
				continue;
			crowns.set(n.c, n.r, n.id);
			visited[i] = true;

			// Execute the kernel. For each neighbour, we're going to check the pixel height
			// to see if it's part of the crown that owns the current node.
			for(size_t j = 0; j < offsetCount; ++j) {

				// Add offsets.
				int cc = n.c + offsets[j][0];
				int rr = n.r + offsets[j][1];

				// Is out of range, continue.
				if(cc < 0 || rr < 0	|| cc >= cols || rr >= rows)
					continue;

				// Calculate the offset pixel index.
				size_t ii = (size_t) rr * (size_t) cols + (size_t) cc;
				if(visited[ii])
					continue;

				// Get the height from the smoothed raster.
				float z = smoothed.get(cc, rr);

				// If is not nodata and is not less than the neighbouring pixel, is valid.
				if(z == nodata || z >= n.z || std::isnan(z))
					continue;

				// Find the index of the threshold corresponding to the given top height.
				size_t idx;
				if(st.find(z, &idx)) {
					const CrownThreshold& ct = thresh[idx];			// Get the crown threshold object.
					float pradius = std::pow(ct.radius / resX, 2);	// The squared radius in pixels.

					float radius = std::pow((float) cc - n.tc, 2) + std::pow((float) rr - n.tr, 2);
					float frac = (n.tz - z) / n.tz;

					// Check that the top meets the threshold:
					if(z < ct.threshold 				// is greater than the min height;
							|| frac > ct.fraction		// is within the height fraction;
							|| radius > pradius)		// is within the radius.
						continue;

					// The pixel is a member of the current crown; add it to the queue.
					q.emplace(n.id, cc, rr, z, n.tc, n.tr, n.tz);
				}
			}
		}
	}

	config.crowns().flush();

	Monitor::get().status(1.0f, "");

}

void crownWriter(CrownDB* db, PolygonContext* pc) {

	int id = 0;
	Top top;
	GEOSGeometry* geom0 = nullptr;
	GEOSGeometry* geom = nullptr;

	while (!Monitor::get().canceled() && (pc->writeRunning || !pc->geoms.empty())) {

		// Get an ID and the list of polys from the queue.
		{
			// Wait for if the queue is empty.
			std::unique_lock<std::mutex> lk(pc->gmtx);
			while (!Monitor::get().canceled() && pc->writeRunning && pc->geoms.empty())
				pc->gcv.wait(lk);
			// If the wakeup is spurious, skip.
			if (pc->geoms.empty())
				continue;
			// Get the ID and geom.
			id = pc->geoms.front().first;
			geom = pc->geoms.front().second;
			pc->geoms.pop_front();
		}

		if (Monitor::get().canceled()) {
			GEOSGeom_destroy_r(pc->gctx, geom);
			continue;
		}

		CrownGeom cgeom0(geom0, nullptr);
		if(db->find(id, top, cgeom0, false)) {
			CrownGeom cgeom(geom, nullptr);
			db->update(top, cgeom);
			if(geom)
				GEOSGeom_destroy_r(pc->gctx, geom);
		} else {
			g_warn("Top " << id << " not found.");
		}

	}
}

void dijital::crowns::util::polygonizeCrowns(CrownsAppConfig& config) {

	if(config.crownsDatabase().empty())
		g_runerr("No crowns database filename given.")

	Monitor::get().status(0.00f, "Crowns: Polygonizing...");

	Band<uint32_t>& crowns = config.crowns();
	CrownDB& cdb = config.crownDB();

	// The polygonization context will be passed into the poly threads.
	PolygonContext pc;
	pc.removeDangles = config.crownsRemoveDangles();
	pc.removeHoles = config.crownsRemoveHoles();
	pc.writeRunning = true;
	pc.dimensions = 3; // Need to set the tree height later.

	// Start the writer thread.
	std::thread th(crownWriter, &cdb, &pc);

	// Create the functor that will accept polygon objects.
	struct fn {
		void operator()(int id, GEOSGeometry* geom, PolygonContext* pc) {
			std::lock_guard<std::mutex> lk(pc->gmtx);
			pc->geoms.emplace_back(id, geom);
			pc->gcv.notify_one();
		}
	};

	// Run the polygonization.
	struct fn f;
	crowns.polygonize(f, &pc);

	// Cleanup any waiting to be written.
	pc.writeRunning = false;
	pc.gcv.notify_all();

	// Wait for the writer to exit and join.
	if (th.joinable())
		th.join();

}

