#include <vector>

#include "process.hpp"
#include "grid.hpp"
#include "settings.hpp"

using namespace tt::proc;
using namespace tt::grid;
using namespace tt::config;

void fixSmoothParams(int& windowSize, float& sigma) {
    if(windowSize < 3)
        windowSize = 3;
    if(windowSize % 2 == 0)
        ++windowSize;
    if(sigma < 0 || std::isnan(sigma))
        sigma = 1;
}

Processor::Processor(Settings* settings) :
        m_settings(settings) {
}

/**
 * Step 1: smooth the grid.
 */
void Processor::smoothGrid() {
    int windowSize = m_settings->get("smoothWindowSize", 0);
    float sigma = m_settings->get("smoothSigma", 1.0f);
    Grid<float> grid(m_settings->get("originalCHM", ""), m_settings->get("originalCHMBand", 1));
    Grid<float> smoothed(grid);
    grid.smooth(smoothed, windowSize, sigma);
    smoothed.save(m_settings->get("smoothedCHM", ""));
}

/**
 * Step 2: find the tops.
 */
void Processor::findTops() {}

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

