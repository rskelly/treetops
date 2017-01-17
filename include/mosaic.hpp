#ifndef __MOSAIC_HPP__
#define __MOSAIC_HPP__

#include <vector>
#include <string>

#include "geotools.hpp"
#include "raster.hpp"
#include "util.hpp"

namespace geotools {

    namespace raster {

        class G_DLL_EXPORT Mosaic {
        private:
            geotools::util::Callbacks *m_callbacks;

        public:
            Mosaic();
            void setCallbacks(geotools::util::Callbacks *callbacks);
            void mosaic(const std::vector<std::string> &files, const std::string &outfile, float distance, int tileSize, int threads = 1);
        };

    } // raster

} // geotools

bool readInput(MemRaster &buf, Raster &input) const {
	int col = iCol - buffer;
	int row = iRow - buffer;
	int cols = tileSize + buffer * 2;
	int rows = tileSize + buffer * 2;
	int cOff = 0;
	int rOff = 0;
	if (col < 0) {
		cOff = -col;
		cols -= col;
		col = 0;
	}
	if (row < 0) {
		rOff = -row;
		rows -= row;
		row = 0;
	}
	if (cols <= 0 || rows <= 0)
		return false;

	input.readBlock(col, row, buf, cOff, rOff, cols, rows);
	return true;
}

#endif
