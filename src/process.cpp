#include <vector>
#include <unordered_set>

#include "process.hpp"
#include "grid.hpp"
#include "settings.hpp"
#include "treetops.hpp"
#include "ds/interval_tree.hpp"
#include "ds/mqtree.hpp"

using namespace tt::proc;
using namespace tt::grid;
using namespace tt::config;
using namespace tt::data;
using namespace tt::ds;

namespace {

	constexpr int intmax = std::numeric_limits<int>::max();
    constexpr int intmin = std::numeric_limits<int>::lowest();

    class px {
    public:
        int col;
        int row;
        int id;
        px(int id, int col, int row) :
            id(id),
            col(col),
            row(row) {}

        px() : px(0, 0, 0) {}
    };

    void fixSmoothParams(int& windowSize, float& sigma) {
        if(windowSize < 3)
            windowSize = 3;
        if(windowSize % 2 == 0)
            ++windowSize;
        if(sigma < 0 || std::isnan(sigma))
            sigma = 1;
    }

    /**
	 * Creates a list of Px (pixel) objects wherein the column and row coordinates are offsets
	 * with respect to the center of the circle at 0, 0. The diameter of the circle
	 * is given by size.
	 * \param size The diameter of the circular kernel.
	 * \returns A list of Px (pixel) objects.
	*/
	std::vector<px> circularKernel(int size, int, bool) {
		std::vector<px> k;
		int rad = std::pow(size / 2, 2);
		for(int r = 0; r < size; ++r) {
			for(int c = 0; c < size; ++c) {
				px& px = k.emplace_back();
				if(std::pow(r - size / 2, 2) + std::pow(c - size / 2, 2) <= rad) {
					px.col = c - size / 2;
					px.row = r - size / 2;
				}
			}
		}
		return k;
	}

    /**
     * Returns true if the pixel at the center of the given circular window is
     * the maximum value in the window. Assigns the max pixel value to max,
     * and the proportion of nulls [0-1] to nulls.
     *
     * \param raster The source raster as vector of floats.
     * \param col The column of interest.
     * \param row The row of interest.
     * \param cols The number of columns in the raster.
     * \param rows The number of rows in the raster.
     * \param window The size of the window.
     * \param maxOut The maximum pixel value in the window.
     * \param nullsOut The proportion of pixels that are null.
     * \return True if the center pixel is the maximum.
     */
    bool isMaxCenter(Grid<float>& grid,
		int col, int row, int window, float nodata,
		float& maxOut, float& nullsOut) {

        int cols = grid.cols();
        int rows = grid.rows();

        float v, max = std::numeric_limits<float>::lowest();
        int mc = -1, mr = -1, n = 0;
        for(int r = -window / 2; r < window / 2 + 1; ++r) {
            for(int c = -window / 2; c < window / 2 + 1; ++c) {
                int cc = col + c;
                int rr = row + r;
                if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows)
                        && (v = grid.get(cc, rr)) != nodata) {
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

    template <class T>
    class Bounds {
    public:
        T lx, ty, rx, by;
        Bounds() : Bounds(intmax, intmin, intmin, intmax) {}
        Bounds(T lx, T ty, T rx, T by) : lx(lx), ty(ty), rx(rx), by(by) {}
        void extend(T x, T y) {
            if(x < lx) lx = x;
            if(x > rx) rx = x;
            if(y > ty) ty = y;
            if(y < by) by = y;
        }
        void set(T lx, T ty, T rx, T by) {
            this->lx = lx; this->ty = ty; this->rx = rx; this->by = by;
        }
        int contains(int x, int y) const {
            return x <= rx && x >= lx && y <= ty && y >= by;
        }
        int width() const {
            return rx - lx;
        }
        int height() const {
            return ty - by;
        }
        int midx() {
            return lx + width() / 2.0;
        }
        int midy() {
            return by + height() / 2.0;
        }
    };    
} // anon

Processor::Processor(Settings* settings) :
        m_settings(settings) {
}

void Processor::run() {

	std::unique_ptr<Grid<float>> grid = std::make_unique<Grid<float>>();
	grid->load(m_settings->get("originalCHM", ""), m_settings->get("originalCHMBand", 1));

	if(m_settings->get("doSmoothing", true)) {
		std::unique_ptr<Grid<float>> smoothed = std::make_unique<Grid<float>>();
		smoothGrid(*grid, *smoothed);
		grid.reset(smoothed.release());
	}

	std::vector<Treetop> tops;
    Grid<int> topsWindowGrid;
    Grid<int> topsIDGrid;
	findTops(*grid, tops, topsWindowGrid, topsIDGrid);

    topsWindowGrid.save(join(tt::util::parent(m_settings->get("originalCHM", "")), "tops_windows.tif"));
    topsIDGrid.save(join(tt::util::parent(m_settings->get("originalCHM", "")), "tops_ids.tif"));
}

/**
 * Step 1: smooth the grid.
 */
void Processor::smoothGrid(Grid<float>& grid, Grid<float>& smoothed) {
    int windowSize = m_settings->get("smoothWindowSize", 0);
    float sigma = m_settings->get("smoothSigma", 1.0f);
    // Copy the properties into the smoothed grid.
    smoothed.copyOther(grid);
    // Apply smoothing.
    grid.smooth(smoothed, windowSize, sigma);
    // Save the output to the smoothed file.
    smoothed.save(m_settings->get("smoothedCHM", ""));
}

/**
 * Step 2: find the tops.
 * This creates:
 * 1) A list of treetop objects.
 * 2) A window grid which stores the size of the window used to locate the top, at the position of each top.
 * 3) An ID grid, which stores the ID of a top at its position.
 */
void Processor::findTops(Grid<float>& grid, std::vector<Treetop>& tops, 
        Grid<int>& topsWindowGrid, Grid<int>& topsIDGrid) {

	// Counter to create unique treetop IDs.
	int topId = 0;
	int cols = grid.cols();
	int rows = grid.rows();
	float nodata = grid.nodata();
    
    topsWindowGrid.copyOther(grid);
    topsIDGrid.copyOther(grid);

    // The maximum proportion of nulls allowed in a window.
    float topsMaxNulls = 0.5; //m_settings->get("topsMaxNulls", 0.5);

	// Build a list of circular window offsets for the window sizes. Find the largest window size.
	int maxWindow = 0;
	std::unordered_map<int, std::vector<px> > circWindows;
	for(const TopThreshold& t : m_settings->topThresholds()) {
		circWindows[t.window] = circularKernel(t.window, 0, false);
		if(t.window > maxWindow)
			maxWindow = t.window;
	}

	// Iterate over the thresholds from lowest height to highest.
	for(const TopThreshold& t : m_settings->topThresholds()) {

        // To store treetops with col/row and window size.
        std::list<Treetop> tops;

        // Iterate over the raster, applying the kernel to find the maxium at the centre.
        float v;
        float max;   // The maximum elevation (i.e., the top height)
        float nulls; // The proportion of null pixels.
        bool isMax;
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                if((v = grid.get(col, row)) >= t.threshold) {
                    isMax = isMaxCenter(grid, col, row, t.window, nodata, max, nulls);
                    if (isMax && nulls <= topsMaxNulls)
                        tops.emplace_back(++topId, col, row, t.window);
                }
            }
        }

        for (const Treetop& top : tops) {
            topsWindowGrid.set(top.col, top.row, top.window);
            topsIDGrid.set(top.col, top.row, top.id);
        }
    }

}

/**
 * Step 3: delineate crowns.
 */
void Processor::delineateCrowns(Grid<float>& grid, Grid<int>& crowns, Grid<int>& ids, Grid<int>& windows) {

	crowns.fill(0);
	float res = grid.xRes();

	// The interval tree keeps track of ranges of completed rows
	int maxRadius = 0;
	MQTree<Treetop> qt;
	IntervalTree<float, size_t> st;
	const std::vector<CrownThreshold>& thresh = m_settings->crownThresholds();
    
	for(size_t i = 0; i < thresh.size(); ++i) {
		st.add(thresh[i].fraction, i);
		if(thresh[i].radius > (maxRadius * res))
			maxRadius = std::ceil(thresh[i].radius / res);
	}

	int cols = grid.cols();
	int rows = grid.rows();
	float nodata = grid.nodata();

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
	std::vector<bool> visited(cols * rows);
	std::fill(visited.begin(), visited.end(), false);

	// Collect unique top IDs to count status.
	size_t topCount = 0;
	std::unordered_set<size_t> topSet;
	{
		Treetop t;
		qt.reset();
		while(qt.next(t)) {
			topSet.insert(t.id);
			++topCount;
		}
	}

	/*
	// Convert the Tops to Nodes, add to work queue.
	std::queue<Treetop> q;
	Bounds<int> bounds;
	Treetop query;
	for(const auto& tile : tiles) {
		// Create a bounding box to search for tops; this is the
		// buffer size, plus a fringe equal to max radius.
		bounds.set(tile.first - maxRadius, tile.second - maxRadius, tile.first + tileSize + maxRadius, tile.second + tileSize + maxRadius);

		// Search radius is the diagonal of the box plus maxRadius.
		float rad = std::sqrt(std::pow(bounds.width() / 2, 2) + std::pow(bounds.height() / 2, 2)) + maxRadius;

		// Search starts at the box's center.
		query.update(0, 0, 0, 0, 0, bounds.midx(), bounds.midy(), 0, 0, 0);

		// Search for the tops.
		std::list<Treetop> tops;
		if(!qt.search(query, rad, std::back_inserter(tops)))
			continue;

		// Enqueue the tops for processing. TODO: Would be nice to search and enqueue in one step.
		for(Treetop& top : tops)
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
	*/
}

/**
 * Step 4: merge crowns.
 */
void Processor::mergeCrowns() {}

/**
 * Step 5: polygonize crowns (and clean up, if configured).
 */
void Processor::polygonizeCrowns() {}

/**
 * Step 6: save outputs.
 */
void Processor::saveOutputs() {

}

