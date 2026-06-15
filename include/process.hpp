#ifndef __PROCESS_HPP__
#define __PROCESS_HPP__

#include "config.hpp"
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
        const tt::config::Config* m_settings;
        std::vector<Treetop> m_tops;

    public:

        explicit Processor(const tt::config::Config&);

        void run();

        void smoothGrid(Grid<float>&, Grid<float>&);

        void findTops(Grid<float>&, std::vector<Treetop>&, Grid<int>&, Grid<int>&);

        void delineateCrowns(std::vector<Treetop>&, Grid<float>&, Grid<int>&, Grid<int>&, Grid<int>&);

        void polygonizeCrowns(const std::vector<Treetop>&, tt::grid::Grid<int>&, tt::vec::CrownDB&);

        void updateTops(std::vector<Treetop>&, tt::grid::Grid<float>&, tt::grid::Grid<int>&);

        void saveOutputs(
            Grid<int>& topsWindowGrid,
            Grid<int>& topsIDGrid,
            Grid<int>& crowns,
            CrownDB& db,
            const std::string& projection);

    };

} // process
} // tt

#endif // __PROCESS_HPP__
