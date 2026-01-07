#include "settings.hpp"

namespace tt {
namespace proc {

    class Processor {
    private:
        tt::config::Settings* m_settings;

    public:

        Processor(tt::config::Settings*);

        /**
         * Step 1: smooth the grid.
         */
        void smoothGrid();

        /**
         * Step 2: find the tops.
         */
        void findTops();

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