/**
 * "Feathers" the edges of a data region (i.e. not nodata) to the specified
 * distance in map units, using the specified curve.
 */
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <atomic>

#include <omp.h>

#include "geotools.hpp"
#include "mosaic.hpp"
#include "util.hpp"
#include "raster.hpp"

using namespace geotools::raster;

namespace geotools {

	namespace raster {

		namespace util {

			/**
			 * Returns a value between 0 and 1 following the tan curve.
			 * The range of step is expected to be 0 -> steps and is clamped.
			 */
			float tanCurve(float step, float steps) {
				step = g_min(steps, g_max(0.0, step));
				return tanh(((step - steps / 2.0) / (steps / 2.0)) * G_PI) * 0.5 + 0.5;
			}

			/**
			 * Returns true if the given pixel is next to a nodata pixel (but not if
			 * it is one), or the edge of the grid.
			 */
			bool isEdgePixel(Grid &fillGrid, int col, int row, int cols, int rows) {
				if (fillGrid.getInt(col, row) == 0)
					return false;
				for (int r = row - 1; r < row + 2; ++r) {
					for (int c = col - 1; c < col + 2; ++c) {
						if (c <= 0 || r <= 0 || c >= cols - 1 || r >= rows - 1
								|| fillGrid.getInt(c, r) == 0)
							return true;
					}
				}
				return false;
			}

			/**
			 * Feathers the edges of data regions in a raster grid by setting the alpha
			 * value for a pixel in proportion to its distance from the nearest null or edge.
			 */
			bool feather(Grid &srcGrid, Grid &dstGrid, float distance,
					float nodata, float resolution) {
				// Fill grid is used to keep track of where the edges are as they're "snowed in"
				// Starts out as a mask of non-nodata pixels from the source.
				int cols = srcGrid.props().cols();
				int rows = srcGrid.props().rows();
				MemRaster fillGrid(srcGrid.props());
				int valid = 0;
				for (size_t i = 0; i < (size_t) rows * cols; ++i) {
					if (srcGrid.getFloat(i) == nodata) {
						fillGrid.setFloat(i, 0);
					} else {
						fillGrid.setFloat(i, 1);
						++valid;
					}
				}

				if (valid == 0)
					return false;

				// The number of steps is just the number of pixels needed to
				// cover the distance of the fade.
				float step = 0.0;
				float steps = g_max(1.0, distance / resolution);
				bool found = false;
				// "Snow in" the alpha mask. The loop will exit when no more edge pixels can be found.
				do {
					found = false;
					for (int row = 0; row < rows; ++row) {
						for (int col = 0; col < cols; ++col) {
							if (isEdgePixel(fillGrid, col, row, cols, rows)) {
								fillGrid.setInt(col, row, 2); // Set edge to dirty.
								dstGrid.setFloat(col, row, tanCurve(step, steps)); // TODO: Configurable curves.
								found = true;
							}
						}
					}
					// Reset dirty edges to 0
					for (size_t i = 0; i < (size_t) rows * cols; ++i)
						if (fillGrid.getInt(i) == 2)
							fillGrid.setInt(i, 0);
					step += 1.0;
				} while (found && step <= steps);

				return true;
			}

			/**
			 * Blends two rasters together using the alpha grid for blending.
			 */
			void blend(Grid &imGrid, Grid &bgGrid, Grid &alpha,
					float imNodata, float bgNodata, int buffer) {
				int cols = imGrid.props().cols();
				int rows = imGrid.props().rows();
				for (int r = buffer; r < rows - buffer; ++r) {
					for (int c = buffer; c < cols - buffer; ++c) {
						float bv = bgGrid.getFloat(c, r);
						float iv = imGrid.getFloat(c, r);
						if (!(bv == bgNodata || iv == imNodata)) {
							float av = alpha.getFloat(c, r);
							bgGrid.setFloat(c, r, bv * (1.0 - av) + iv * av);
						}
					}
				}
			}

			int __tile_id = 0;

			class Tile {
			public:
				int id;
				int tileSize;
				int buffer;
				int iCol;
				int iRow;
				int oCol;
				int oRow;

				Tile(int tileSize, int buffer, int iCol, int iRow, int oCol, int oRow) :
						id(++__tile_id), tileSize(tileSize), buffer(buffer), iCol(iCol), iRow(
								iRow), oCol(oCol), oRow(oRow) {
				}

				bool readInput(geotools::raster::MemRaster &buf, geotools::raster::Raster &input) const {
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

					input.writeToBlock(buf, cols, rows, col, row, cOff, rOff);
					return true;
				}

				bool readOutput(MemRaster &buf, Raster &output) const {
					int col = oCol - buffer;
					int row = oRow - buffer;
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
					buf.writeToBlock(output, cols, rows, col, row, cOff, rOff);
					return true;
				}

				void writeOutput(MemRaster &buf, Raster &output) const {
					if (!(oCol >= output.props().cols() || oRow >= output.props().rows()))
						buf.writeToBlock(output, tileSize, tileSize, oCol, oRow, buffer, buffer);
				}

				void print() const {
					std::cerr << "[Tile: in: " << iCol << "," << iRow << "; out: " << oCol
							<< "," << oRow << "]" << std::endl;
				}
			};

		} // util

		Mosaic::Mosaic() :
				m_callbacks(nullptr) {
		}

		void Mosaic::setCallbacks(geotools::util::Callbacks *callbacks) {
			m_callbacks = callbacks;
		}

		/**
		 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
		 * mask for the others. The spatial reference system and resolution of all layers must match.
		 * The distance argument determines the distance over which the edges are feathered. The
		 * overviews argument, if true, forces the construction of new overviews. If this is not
		 * used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
		 */
		// TODO: Background must be larger than other files. Remedy this.
		void Mosaic::mosaic(const std::vector<std::string> &files,
				const std::string &outfile, float distance, int tileSize, int threads) {

			g_debug("mosaic");

			using namespace geotools::raster::util;

			if (m_callbacks) {
				m_callbacks->stepCallback(0.0f);
				m_callbacks->overallCallback(0.0f);
			}

			if (distance <= 0.0)
				g_argerr("Invalid distance: " << distance);
			if (outfile.size() == 0)
				g_argerr("No output file given.");
			if (files.size() < 2)
				g_argerr("Less than 2 files. Nothing to do.");
			if (tileSize <= 0)
				g_argerr("Tile size must be greater than zero.");

			if (threads > 0) {
				g_debug("running with " << threads << " threads");
				omp_set_dynamic(1);
				omp_set_num_threads(threads);
			} else {
				g_argerr("Run with >=1 thread.");
			}

			// Open the BG file for reading only.
			g_debug(" -- opening base file.");
			Raster base(files[0]);

			// Check some properties of other files for compatibility.
			for (unsigned int i = 1; i < files.size(); ++i) {
				Raster check(files[i]);
				if (check.props().resolutionX() != base.props().resolutionX() ||
						check.props().resolutionY() != base.props().resolutionY())
				g_argerr("Resolution of " << files[i] << " doesn't match base.")
			}

			// Create the destination file to modify it.
			Raster output(outfile, base.props());
			{
				GridProps pr(base.props());
				pr.setSize(base.props().cols(), 1000);
				MemRaster buf(pr);
				g_debug(" -- writing base file to output.");
				// Copy by block to avoid memory problems on big rasters.
				// TODO: Configurable block size.
				for (int r = 0; r < base.props().rows(); r += pr.rows()) {
					base.writeToBlock(buf, 0, 0, 0, r);
					buf.writeToBlock(output, 0, 0, 0, r);
				}
			}

			for (unsigned int i = 1; i < files.size(); ++i) {

				if (m_callbacks)
					m_callbacks->overallCallback((i - 0.5) / (files.size() - 1));

				Raster input(files[i]);

				Bounds bounds = base.props().bounds().intersection(input.props().bounds());

				// Compute the column/row bounds of the intersection.
				const GridProps& ipr = input.props();
				const GridProps& bpr = base.props();
				int iStartCol = ipr.toCol(ipr.resolutionX() > 0 ? bounds.minx() : bounds.maxx());
				int iStartRow = ipr.toRow(ipr.resolutionY() > 0 ? bounds.miny() : bounds.maxy());
				int iEndCol = ipr.toCol(ipr.resolutionX() > 0 ? bounds.maxx() : bounds.minx());
				int iEndRow = ipr.toRow(ipr.resolutionY() > 0 ? bounds.maxy() : bounds.miny());
				int oStartCol = bpr.toCol(bpr.resolutionX() > 0 ? bounds.minx() : bounds.maxx());
				int oStartRow = bpr.toRow(bpr.resolutionY() > 0 ? bounds.miny() : bounds.maxy());

				// Get the minimum absolute resolution.
				double res = g_min(g_abs(base.props().resolutionX()),
						g_abs(base.props().resolutionY()));
				g_debug(" -- resolution is: " << res);

				// Snap the distance to resolution
				distance = std::floor(distance / res) * res;
				g_debug(" -- snapped distance is " << distance);

				// The number of rows required to accomodate the fade distance.
				int buffer = (int) distance / res + 1;
				while (tileSize <= buffer) {
					g_warn(
							"The tileSize (" << tileSize
									<< ") must be larger than the buffer distance ("
									<< buffer << "). Doubling.");
					tileSize *= 2;
				}
				g_debug(" -- tile size is: " << tileSize);

				std::vector<std::unique_ptr<Tile> > tiles;

				for (int ir = iStartRow, rr = oStartRow; ir <= iEndRow;
						ir += tileSize, rr += tileSize) {
					for (int ic = iStartCol, oc = oStartCol; ic <= iEndCol; ic +=
							tileSize, oc += tileSize) {
						std::unique_ptr<Tile> t(
								new Tile(tileSize, buffer, ic, ir, oc, rr));
						tiles.push_back(std::move(t));
					}
				}

				std::atomic<int> tileStatus(0);

				#pragma omp parallel shared(tileStatus)
				{
					// Initialize the grids.
					GridProps pr;
					pr.setSize(tileSize + buffer * 2, tileSize + buffer * 2);
					MemRaster inGrid(pr);
					MemRaster outGrid(pr);
					MemRaster alphaGrid(pr);
					float outNodata = output.props().nodata();
					float inNodata = input.props().nodata();

					#pragma omp for nowait
					for (size_t t = 0; t < tiles.size(); ++t) {

						tileStatus++;

						if (m_callbacks)
							m_callbacks->stepCallback((tileStatus - 0.5) / tiles.size());

						std::unique_ptr<Tile> tile = std::move(tiles[t]);

						g_debug(" -- tile " << t << ", " << tile->id << "; thread: " << omp_get_thread_num());

						// Set to zero or nodata depending on the band.
						outGrid.fillFloat(outNodata);
						inGrid.fillFloat(inNodata);
						alphaGrid.fillFloat(1.0);

						// Load foreground
						if (!tile->readInput(inGrid, input))
							continue;

						// Feather returns true if it did any work.
						if (!feather(inGrid, alphaGrid, distance, inNodata, res))
							continue;

						// Load background
						if (!tile->readOutput(outGrid, output))
							continue;

						// Blend the fore/background images with alpha
						blend(inGrid, outGrid, alphaGrid, inNodata, outNodata, buffer);

						// Write to output.
						tile->writeOutput(outGrid, output);

						if (m_callbacks)
							m_callbacks->stepCallback(
									(float) tileStatus / tiles.size());

					}

				} // parallel

				if (m_callbacks) {
					m_callbacks->stepCallback(1.0);
					m_callbacks->overallCallback((float) i / (files.size() - 1));
				}
			}

			g_debug(" -- done.");

		}

	} // raster

} // geotools

