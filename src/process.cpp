#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "process.hpp"
#include "grid.hpp"
#include "vector.hpp"
#include "settings.hpp"
#include "treetops.hpp"
#include "ds/interval_tree.hpp"
#include "ds/mqtree.hpp"
#include "ds.hpp"

using namespace tt::proc;
using namespace tt::grid;
using namespace tt::config;
using namespace tt::data;
using namespace tt::ds;
using namespace tt::vec;

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

Processor::Processor(Settings& settings) :
        m_settings(&settings) {
}

void Processor::run() {

	std::unique_ptr<Grid<float>> grid;			// Source grid.
	std::unique_ptr<Grid<float>> smoothed;		// Smoothed grid.
	Grid<float>* working;						// Points to the smoothed or original grid, depending on settings.

	grid = std::make_unique<Grid<float>>();
	grid->load(m_settings->get("originalCHM", ""), m_settings->get("originalCHMBand", 1));

	if(m_settings->get("doSmoothing", true)) {
		smoothed = std::make_unique<Grid<float>>();
		smoothGrid(*grid, *smoothed);
		working = smoothed.get();
	} else {
		working = grid.get();
	}

	std::vector<Treetop> tops;
    Grid<int> topsWindowGrid;
    Grid<int> topsIDGrid;
	findTops(*working, tops, topsWindowGrid, topsIDGrid);

    topsWindowGrid.save(join(tt::util::parent(m_settings->get("originalCHM", "")), "tops_windows.tif"));
    topsIDGrid.save(join(tt::util::parent(m_settings->get("originalCHM", "")), "tops_ids.tif"));

	Grid<int> crowns(topsIDGrid);
	delineateCrowns(tops, *working, crowns, topsIDGrid, topsWindowGrid);
	crowns.save(join(tt::util::parent(m_settings->get("originalCHM", "")), "crowns.tif"));

	if(true /* settings -> update tops */)
		updateTops(tops, *grid, crowns);

	CrownDB db;
	polygonizeCrowns(tops, crowns, db);
	db.saveCrowns(join(tt::util::parent(m_settings->get("originalCHM", "")), "crowns.sqlite"), "Spatialite", "crowns", crowns.crs());
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
                        tops.emplace_back(++topId, col, row, t.window, max);
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
void Processor::delineateCrowns(std::vector<Treetop>& tops, Grid<float>& grid, Grid<int>& crowns, Grid<int>& ids, Grid<int>& windows) {

	crowns.fill(0);
	float res = grid.xRes();

	// The interval tree keeps track of ranges of completed rows
	int maxRadius = 0;
	MQTree<Treetop> qt;
	IntervalTree<float, int> st;
	const std::vector<CrownThreshold>& thresh = m_settings->crownThresholds();
    
	for(int i = 0; i < thresh.size(); ++i) {
		st.add(thresh[i].fraction, i);
		if(thresh[i].radius > (maxRadius * res))
			maxRadius = std::ceil(thresh[i].radius / res);
	}

	int cols = grid.cols();
	int rows = grid.rows();
	float nodata = grid.nodata();

	// Build the list of offsets for D8 search.
	int offsetCount = 8;
	int offsets[][2] = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};

	// A list to track visited pixels.
	std::vector<bool> visited(cols * rows);
	std::fill(visited.begin(), visited.end(), false);

	// Convert the Tops to Nodes, add to work queue.
	std::queue<Node> q;
	Treetop query;

	// Enqueue the tops for processing. TODO: Would be nice to search and enqueue in one step.
	for(Treetop& top : tops)
		q.emplace(top);

	// Run through the queue.
	while(!q.empty()) {

		// Get the next node.
		Node n = std::move(q.front());
		q.pop();

		// Calculate pixel index; set the ID at the current pixel in the crowns raster.
		int i = n.r * cols + n.c;
		if(visited[i])
			continue;
		crowns.set(n.c, n.r, n.id);
		visited[i] = true;

		// Execute the kernel. For each neighbour, we're going to check the pixel height
		// to see if it's part of the crown that owns the current node.
		for(int j = 0; j < offsetCount; ++j) {

			// Add offsets.
			int cc = n.c + offsets[j][0];
			int rr = n.r + offsets[j][1];

			// Is out of range, continue.
			if(cc < 0 || rr < 0	|| cc >= cols || rr >= rows)
				continue;

			// Calculate the offset pixel index.
			int ii = rr * cols + cc;
			if(visited[ii])
				continue;

			// Get the height from the smoothed raster.
			float z = grid.get(cc, rr);

			// If is not nodata and is not less than the neighbouring pixel, is valid.
			if(z == nodata || z >= n.z || std::isnan(z))
				continue;

			// Find the index of the threshold corresponding to the given top height.
			int idx;
			if(st.find(z, &idx)) {
				const CrownThreshold& ct = thresh[idx];			// Get the crown threshold object.
				float pradius = std::pow(ct.radius / res, 2);	// The squared radius in pixels.

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

/**
 * Step 4: update tops with the max height from the unsmoothed raster within the delineated crown.
 */
void Processor::updateTops(const std::vector<Treetop>& tops, Grid<float>& grid, Grid<int>& crowns) {

	std::unordered_map<int, Treetop*> topMap;
	for(const Treetop& top : tops)
		topMap.emplace(top.id, &top);

	for(int r = 0; r < crowns.rows(); ++r) {
		for(int c = 0; c < crowns.cols(); ++c) {
			int id = crowns.get(c, r);
			if(topMap.contains(id)) {
				Treetop* top = topMap.at(id);
				float v = grid.get(c, r);
				if(v != grid.nodata() && v > top->oz) {
					top->oz = v;
					top->ox = grid.toX(c);
					top->oy = grid.toY(r);
				}
			}
		}
	}

}

/**
 * Step 5: polygonize crowns (and clean up, if configured).
 */
void Processor::polygonizeCrowns(const std::vector<Treetop>& tops, Grid<int>& crowns, CrownDB& db) {

	// The polygonization context will be passed into the poly threads.
	PolyCtx pc;
	pc.removeDangles = m_settings->get("crownsRemoveDangles", false);
	pc.removeHoles = m_settings->get("crownsRemoveHoles", false);
	pc.running = true;
	pc.dimensions = 3; // Need to set the tree height later.

	// Create the functor that will accept polygon objects.
    // TODO: This will be important for threading.
	PolyMergeCallback callback;
	crowns.polygonize(&callback, &pc);

	std::unordered_map<int, const Treetop*> topMap;
	for(const Treetop& top : tops)
		topMap[top.id] = &top;

	for(std::pair<int, GEOSGeometry*>& geom : pc.geoms) {
		const Treetop* top = topMap[geom.first];
		db.insert(*top, geom.second, pc.gctx);
	}

	// Cleanup any waiting to be written.
	pc.running = false;

}


/**
 * Step 5: save outputs.
 */
void Processor::saveOutputs() {

}

