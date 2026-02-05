#include "settings.hpp"
#include "treetops.hpp"
#include "grid.hpp"

using namespace tt::data;
using namespace tt::grid;

namespace tt {
namespace proc {

    class Processor {
    private:
        tt::config::Settings* m_settings;
        std::vector<Treetop> m_tops;

    public:

        Processor(tt::config::Settings*);

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
        void delineateCrowns(Grid<float>&, Grid<int>&, Grid<int>&, Grid<int>&);

        /**
         * Step 4: merge crowns.
         */
        void mergeCrowns();

        /**
         * Step 5: polygonize crowns (and clean up, if configured).
         */
        void polygonizeCrowns();

        /**
         * Step 6: save outputs.
         */
        void saveOutputs();

    };

} // process
} // tt