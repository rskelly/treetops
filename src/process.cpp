#include <vector>

#include "process.hpp"
#include "grid.hpp"
#include "settings.hpp"


using namespace tt::proc;
using namespace tt::grid;
using namespace tt::config;
using namespace tt::data;

namespace {

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
        for(int r = -window / 2; r < window / 2; ++r) {
            for(int c = -window / 2; c < window / 2; ++c) {
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
                        tops.emplace_back(col, row, t.window);
                }
            }
        }

        for (const Treetop& top : tops) {
            ++topId;
            topsWindowGrid.set(top.col, top.row, top.window);
            topsIDGrid.set(top.col, top.row, topId);
        }
    }

}

/**
 * Step 3: delineate crowns.
 */
void Processor::delineateCrowns() {}

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

