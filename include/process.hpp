#include "settings.hpp"
#include "treetops.hpp"
#include "grid.hpp"


namespace tt {
namespace proc {

    class Processor {
    private:
        tt::config::Settings* m_settings;
        std::vector<tt::data::Treetop> m_tops;

    public:

        Processor(tt::config::Settings*);

        void run();

        /**
         * Step 1: smooth the grid.
         */
        void smoothGrid(tt::grid::Grid<float>&, tt::grid::Grid<float>&);

        /**
         * Step 2: find the tops.
         */
        void findTops(tt::grid::Grid<float>&, std::vector<tt::data::Treetop>&, tt::grid::Grid<int>&, tt::grid::Grid<int>&);

        /**
         * Step 3: delineate crowns.
         */
        void delineateCrowns();

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