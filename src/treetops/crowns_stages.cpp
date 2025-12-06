/*
* crowns_stages.cpp
*
* Contains high-level processing stages for the top and crown delineation program. 
*
* Created on: Apr 1, 2024
* Author: Rob Skelly <rob@dijital.ca>
*/

#include "geo.hpp"

#include <atomic>
#include <fstream>

#include "ds/mqtree.hpp"
#include "crowns.hpp"
#include "util.hpp"
#include "grid.hpp"

using namespace dijital::ds;
using namespace dijital::util;
using namespace dijital::grid;
using namespace dijital::crowns;
using namespace dijital::crowns::config;
using namespace dijital::crowns::util;

namespace {

	/**
	 * Creates a list of Px (pixel) objects wherein the column and row coordinates are offsets
	 * with respect to the center of the circle at 0, 0. The diameter of the circle
	 * is given by size.
	 * \param size The diameter of the circular kernel.
	 * \returns A list of Px (pixel) objects.
	*/
	std::vector<Px> circularKernel(int size, int, bool) {
		std::vector<Px> k;
		int rad = std::pow(size / 2, 2);
		for(int r = 0; r < size; ++r) {
			for(int c = 0; c < size; ++c) {
				Px& px = k.emplace_back();
				if(std::pow(r - size / 2, 2) + std::pow(c - size / 2, 2) <= rad) {
					px.col = c - size / 2;
					px.row = r - size / 2;
				}
			}
		}
		return k;
	}

} // anon

/**
 * Find the maxima, create and store tops.
 * \param config The configuration object.
 */
void dijital::crowns::util::stage1(CrownsAppConfig& config) {

	Band<float>& smoothed = config.smoothed();
	Band<uint32_t>& topsWindowGrid = config.topsWindow();
	Band<uint32_t>& topsIDGrid = config.topsID();

	if(!config.smoothExisted()) {
		Band<float>& chm = config.chm();
		chm.smooth(smoothed, config.smoothSigma(), config.smoothWindowSize());
		smoothed.save(config.smoothedCHM(), config.smoothedCHMDriver());
		config.setSmoothExisted(true);
	}

	// Counter to create unique treetop IDs.
	int topId = 0;
	int cols = smoothed.props().cols();
	int rows = smoothed.props().rows();
	float nodata = smoothed.props().nodata();

	// Build a list of circular window offsets for the window sizes. Find the largest window size.
	int maxWindow = 0;
	std::unordered_map<int, std::vector<Px> > circWindows;
	for(const TopThreshold& t : config.topsThresholds()) {
		circWindows[t.window] = circularKernel(t.window, 0, false);
		if(t.window > maxWindow)
			maxWindow = t.window;
	}

	// The tile buffer size is 1/2 of the largest window; added to each edge of the tile.
	int bufSize = (int) std::ceil((float) maxWindow / 2);
	int tileSize = std::min(4096, std::max(cols, rows));
	int tileBufSize = tileSize + bufSize * 2;
	std::vector<float> tile(std::pow(tileBufSize, 2));

	std::mutex mtx;
	size_t status = 0;
	size_t total = std::max(1, (cols / tileSize * rows / tileSize)) * config.topsThresholds().size();

	// Iterate over the thresholds from lowest height to highest.
	for(const TopThreshold& t : config.topsThresholds()) {

		// Iterate over the tiles.
		for(int tr = 0; tr < rows; tr += tileSize) {
			for(int tc = 0; tc < cols; tc += tileSize) {

				smoothed.getTile(tile.data(), tc, tr, tileSize, tileSize, bufSize, bufSize);

				Monitor::get().status((float) ++status / total, "Finding tops.");

				if (Monitor::get().canceled())
					continue;

				// To store treetops with col/row and window size.
				std::list<TopConfig> tops;

				// Iterate over the raster, applying the kernel to find the maxium at the centre.
				float v;
				float max;   // The maximum elevation (i.e., the top height)
				float nulls; // The proportion of null pixels.
				bool isMax;
				for (int row = 0; row < tileSize; ++row) {
					for (int col = 0; col < tileSize; ++col) {
						int cc = tc + col;
						int rr = tr + row;
						if(cc < 0 || rr < 0 || cc >= cols || rr >= rows)
							continue;
						int c = col + bufSize;
						int r = row + bufSize;
						if((v = tile[r * tileBufSize + c]) >= t.threshold) {
							isMax = isMaxCenter(tile, c, r, tileBufSize, tileBufSize, t.window, nodata, max, nulls);
							if (isMax && nulls <= config.topsMaxNulls())
								tops.emplace_back(cc, rr, t.window);
						}
					}
				}

				{
					std::lock_guard<std::mutex> lk(mtx);
					for (const TopConfig& top : tops) {
						++topId;
						topsWindowGrid.set(top.col, top.row, top.window);
						topsIDGrid.set(top.col, top.row, topId);
					}
				}
			}
		}
	}
}

/**
 * Locates the parent treetops and assigns the IDs.
 * \param config The configuration object.
 */
void dijital::crowns::util::stage2(CrownsAppConfig& config) {

	Band<uint32_t>& topsWindowGrid = config.topsWindow();
	Band<uint32_t>& topsParentGrid = config.topsParent();
	Band<uint32_t>& topsIDGrid = config.topsID();

	const GridProps& props = topsWindowGrid.props();
	int rows = props.rows();
	int cols = props.cols();

	// To keep track of progress.
	int status = 0;

	// The smallest window used in this job.
	int minWindow = config.topsThresholds()[0].window;
	// Build a list of circular window offsets for the window sizes.
	std::unordered_map<int, std::vector<Px> > circWindows;
	for(const TopThreshold& t : config.topsThresholds())
		circWindows[t.window] = circularKernel(t.window / 2, 0, false);

	for(int row = 0; row < rows; ++row) {

		Monitor::get().status((float) ++status / rows, "Finding parent tree tops.");

		if (Monitor::get().canceled())
			continue;

		// Collect the locations of tops with a > minimum window size.
		std::list<TopConfig> windows;
		{
			int window;
			for(int col = 0; col < cols; ++col) {
				if((window = topsWindowGrid.get(col, row)) > minWindow)
					windows.emplace_back(col, row, window);
			}
		}

		// Iterate over the windows, searching each for other tops.
		for(const TopConfig& top : windows) {
			// Iterate over the col/row offsets for a window of this size.
			for(const Px& px : circWindows[top.window]) {
				int c = top.col + px.col;
				int r = top.row + px.row;
				if(c >= 0 && r >= 0 && c < cols && r < rows) {
					int window0 = topsWindowGrid.get(c, r);
					if(window0 > 0 && window0 < top.window && !topsParentGrid.get(c, r))
						topsParentGrid.set(c, r, topsIDGrid.get(top.col, top.row)); // TODO: This was cc, r; Save the ID of the parent.
				}
			}
		}
	}
}

//  Creates the treetop objects using the rasters as inputs.
void dijital::crowns::util::stage3(CrownsAppConfig& config) {

	mqtree<dijital::crowns::util::Top>& qt = config.topsTree();
	Band<float>& smoothed = config.smoothed();
	Band<uint32_t>& topsWindowGrid = config.topsWindow();
	Band<uint32_t>& topsParentGrid = config.topsParent();
	Band<uint32_t>& topsIDGrid = config.topsID();

	const GridProps& props = smoothed.props();
	int rows = props.rows();
	int cols = props.cols();
	int count = 0;
	std::atomic<int> status(0);

	for(int row = 0; row < rows; ++row) {

		Monitor::get().status((float) ++status / rows, "Creating tree top objects.");

		if(Monitor::get().canceled())
			continue;

		for(int col = 0; col < cols; ++col) {
			int window;
			if((window = topsWindowGrid.get(col, row))) {
				double z = smoothed.get(col, row);
				uint32_t topID = topsIDGrid.get(col, row);
				uint32_t parentID = topsParentGrid.get(col, row);
				double x = props.toX(col);
				double y = props.toY(row);
				qt.add(dijital::crowns::util::Top(
					topID, parentID,
					0, 0, 0, 	// Can't fill these until crowns are available
					x, y, z,
					0,			// TODO: Ground z.
					col, row
				));
				++count;
			}
		}

	}

}
