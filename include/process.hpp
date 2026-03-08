#ifndef __PROCESS_HPP__
#define __PROCESS_HPP__

#include "settings.hpp"
#include "treetops.hpp"
#include "grid.hpp"
#include "vector.hpp"

using namespace tt::data;
using namespace tt::grid;
using namespace tt::vec;

namespace tt {
namespace proc {

    class Processor {
    private:
        tt::config::Settings* m_settings;
        std::vector<Treetop> m_tops;

    public:

        Processor(tt::config::Settings&);

        void run();

        /**
         * Step 1: smooth the grid.
         */
        void smoothGrid(Grid<float>&, Grid<float>&);

        /**
         * Step 2: find the tops.
         */
        void findTops(Grid<float>&, std::vector<Treetop>&, Grid<int>&, Grid<int>&);

        /**
         * Step 3: delineate crowns.
         * Source grid, crown grid, ID grid, window grid.
         */
        void delineateCrowns(std::vector<Treetop>&, Grid<float>&, Grid<int>&, Grid<int>&, Grid<int>&);

        /**
         * Step 4: polygonize crowns (and clean up, if configured).
         */
        void polygonizeCrowns(const std::vector<Treetop>&, tt::grid::Grid<int>&, tt::vec::CrownDB&);

        /**
         * Update the top object with the max height in the crown from the unsmoothed raster.
         */
        void updateTops(std::vector<Treetop>&, tt::grid::Grid<float>&, tt::grid::Grid<int>&);

        /**
         * Step 5: save outputs.
         */
        void saveOutputs();

    };

} // process
} // tt

#endif // __PROCESS_HPP__