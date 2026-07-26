#ifndef __GRID_HPP__
#define __GRID_HPP__

#include <iostream>
#include <string>
#include <concepts>
#include <inttypes.h>
#include <cstdlib>
#include <unordered_set>
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/mman.h>
#endif

#include <geos_c.h>

#include <gdal/gdal_priv.h>
#include <gdal/ogr_spatialref.h>
#include <gdal/ogr_geometry.h>
#include <gdal/ogr_feature.h>
#include <gdal/ogrsf_frmts.h>

#include "util.hpp"

#define GRID_MMAP_THRESHOLD 1000000 // Change to mmap when the array size is larger than this.

using namespace tt::util::vec;

namespace tt {
namespace grid {

    /**
     * Stores a square grid of data.
     */
    template<class T>
    class Tile {
    private: 
        T* m_data;          // Allocated data.
        int m_size;    // The length of one size of the tile, in cells.

        /**
         * Check for initialization and coordinate validity.
         */
        void check(int col, int row) {
            if(!m_data)
                throw "This tile is not initialized.";
            if(col < 0 || row < 0 || col >= m_size || row >= m_size)
                throw "Coordinate is out of bounds.";
        }

    public:

        /**
         * Construct. resize must be called before use or an exception is thrown.
         */
        Tile() :
            m_data(nullptr) {
        }

        /**
         * Construct and initialize with the given size.
         */
        Tile(int size) :
            Tile() {
            resize(size);
        }

        /**
         * (Re)size the grid by allocating memory of the right size.
         */
        void resize(int size) {
            if(size == m_size)
                return;
            if(m_data)
                free(m_data);
            m_data = malloc(sizeof(T) * m_size * m_size);
            if(!m_data)
                throw "Failed to allocate memory for tile.";
        }

        /**
         * Get the value at a coordinate.
         */
        T get(int col, int row) {
            check(col, row);
            return m_data[row * m_size + col];
        }

        /**
         * Set the value at a coordinate.
         */
        void set(T value, int col, int row) {
            check(col, row);
            m_data[row * m_size + col] = value;
        }

        /**
         * Free the tile data.
         */
        ~Tile() {
            if(m_data)
                free(m_data);
        }
    };

    /**
     * Represents a raster grid. Loaded from and saved to a rater file.
     * Decides whether to use mapped or online memory depending on the size of the raster relative to
     * GRID_MMAP_THRESHOLD.
     * 
     * Provides methods for smoothing, etc.
     */
    template <class T>
    class Grid {
    private:
        T* m_grid;                              // The main storage location. May be in-memory or mapped. Mapping is used when the grid is larger than GRID_MMAP_THRESHOLD.
        bool m_mapped;                          // True if mapped.
        int m_cols;                             // The number of columns.
        int m_rows;                             // The number of rows.
        int m_band;                             // Raster band. Starts with 1.
        GDALDataType m_type;                    // The GDAL data type of the raster.
        std::vector<double> m_transform;        // The GDAL transform.
        std::string m_crs;                      // The CRS as a WKT string.
        std::vector<std::string> m_bandMeta;    // Band metadata.
        T m_nodata;                             // Value to use for nodata.
        std::string m_driver;                   // The raster driver (e.g., "GTiff").

        /**
         * Free or unmap the allocated memory.
         */
        void freeGrid() {
            if(m_grid) {
                if(m_mapped) {
#if defined(_WIN32)
                    std::free(m_grid);
#else
                    munmap(m_grid, sizeof(T) * m_rows * m_cols);
#endif
                } else {
                    std::free(m_grid);
                }
            }
            m_cols = 0;
            m_rows = 0;
        }

        /**
         * Initialize the memory for the grid. Free previously allocated memory 
         * if needed. If the size of the data segment is larger than GRID_MMAP_THRESHOLD,
         * memory is mapped, otherwise allocated.
         */
        void initGrid(int cols, int rows) {
            freeGrid();
            // Check if the size threshold is exceeded. If so, use mmap.
            int size = sizeof(T) * cols * rows;
#if defined(_WIN32)
            m_mapped = false;
#else
            m_mapped = size > GRID_MMAP_THRESHOLD;
#endif
            if(m_mapped) {
                // Map the data segment.
                m_grid = (T*) mmap(
                    nullptr,
                    sizeof(T) * cols * rows, 
                    PROT_READ|PROT_WRITE, 
                    MAP_PRIVATE|MAP_ANONYMOUS, 
                    -1, 
                    0
                );
                if(m_grid == MAP_FAILED) {
                    m_grid = nullptr;
                    throw "Failed to allocate data for grid.";
                }
            } else {
                // Allocate the data segment.
                m_grid = (T*) std::malloc(size);
                if(m_grid == nullptr)
                    throw "Failed to allocate data for grid.";
            }
            m_cols = cols;
            m_rows = rows;
            fill(m_nodata);
        }

    public:

        /**
         * Construct an empty grid with size 0.
         */
        Grid() :
                m_grid((T*) nullptr),
                m_mapped(false),
                m_cols(0),
                m_rows(0),
                m_band(1),
                m_nodata(-9999),
                m_type(GDALDataType::GDT_Float32) {
            
            m_transform.resize(6);
            
            if(std::is_same<T, float>::value) {
                m_type = GDALDataType::GDT_Float32;
            } else if(std::is_same<T, int>::value) {
                m_type = GDALDataType::GDT_Int32;
            } else {
                throw std::runtime_error("Only float and int are accepted types.");
            }
        }

        /**
         * Copy the properties of the other grid to this one, but no data.
         */
        Grid(const Grid<T>& other) : Grid() {
            copyOther(other);
        }

        T nodata() const {
            return m_nodata;
        }
        
        const std::vector<double>& transform() {
            return m_transform;
        }

        const std::string& crs() const {
            return m_crs;
        }

        void fill(T v) {
            for(int i = 0; i < m_cols * m_rows; ++i)
                m_grid[i] = v;
        }

        double xRes() {
            return m_transform[1];
        }

        double yRes() {
            return m_transform[5];
        }

        double toX(int col) {
            return m_transform[0] + col * m_transform[1];
        }

        double toY(int row) {
            return m_transform[3] + row * m_transform[5];
        }

        void copyBounds(std::vector<float>& bounds) {
            bounds.resize(4);
            int x1 = 0, x2 = 2, y1 = 3, y2 = 1;
            if(xRes() < 0) {
                x1 = 2;
                x2 = 0;
            }
            if(yRes() < 0) {
                y1 = 1;
                y2 = 3;
            }
            bounds[x1] = m_transform[0];
            bounds[x2] = bounds[0] + cols() * xRes();
            bounds[y1] = m_transform[3];
            bounds[y2] = bounds[1] + rows() * yRes();

        }

        /**
         * Copy this grid's affine transform to another array. The other array
         * must be initialized with 6 elements.
         */
        void copyTransform(std::vector<double>& trans) const {
            trans.resize(6);
            for(int i = 0; i < 6; ++i)
                trans[i] = m_transform[i];
        }

        /**
         * Copy the properties of the other grid to this one, but no data.
         */
        template <class U>
        void copyOther(const Grid<U>& other) {
            m_crs = other.crs();
            other.copyTransform(m_transform);
            initGrid(other.cols(), other.rows());
        }

        /**
         * Initialize a grid of the given size.
         */
        Grid(int cols, int rows) :  Grid() {
            initGrid(cols, rows);
        }

        /**
         * Initialize the grid and load the file.
         */
        Grid(const std::string& filename, int band=1) : Grid() {
            load(filename, band);
        }

        ~Grid() {
            freeGrid();
        }

        int cols() const {
            return m_cols;
        }

        int rows() const {
            return m_rows;
        }

        /**
         * Set the value at the given cell in the grid.
         */
        void set(int col, int row, T v) {
            if(col >= 0 && col < m_cols && row >= 0 && row < m_rows)
                m_grid[row * m_cols + col] = v;
        }

        /**
         * Get the value at the given cell in the grid. If it's out of bounds return nodata.
         */
        T get(int col, int row) {
            if(col >= 0 && col < m_cols && row >= 0 && row < m_rows)
                return m_grid[row * m_cols + col];
            return m_nodata;
        }

        /**
         * Calculate an array of weights for Gaussian smoothing.
         */
        void gaussianWeights(std::vector<float>& weights, int window, double sigma, double mean = 0) const {
            if (sigma <= 0)
                throw std::runtime_error("Sigma must be > 0.");
            if (window < 3)
                throw std::runtime_error("Kernel size must be 3 or larger.");
            if (window % 2 == 0) {
                ++window;
                std::cerr << "Gaussian kernel size must be an odd number >=3. Bumping up to " << window;
            }
            for (int r = 0; r < window; ++r) {
                for (int c = 0; c < window; ++c) {
                    int x = c - window / 2;
                    int y = r - window / 2;
                    weights[r * window + c] = (1 / (2.0 * M_PI * sigma * sigma)) * std::pow(M_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
                }
            }
        }
        
        /**
         * Read a rectangular region of cells into the given vector.
         */
        void read(std::vector<T>& tile, int col, int row, int w, int h) {
            tile.resize(w, h);
            // If the column + width is larger than the available width, trim it.
            if(col + w > m_cols || w < 1 || col < 0)
                _runerr("Invalid column or width: " << col << "; " << w);
            if(row + h > m_rows || h < 1 || row < 0)
                _runerr("Invalid row or height: " << row << "; " << h);
            for(int r = row; r < row + h; ++r) {
                for(int c = col; c < col + w; ++c) {
                    float v = m_nodata;
                    if(r >= 0 && r < m_rows && c >= 0 && c < m_cols)
                        v = *(m_grid + (r * m_cols + c));
                    tile[(r - row) * w + (c - col)] = v;
                }
            }
        }

        /**
         * Reset all pixels to nodata.
         */
        void clear() {
            for(int i = 0; i < m_rows * m_cols; ++i)
                m_grid[i] = m_nodata;
        }

        /**
         * Smooth the grid into the given instance using the given parameters to the Gaussian kernel.
         */
        void smooth(Grid<T>& smoothed, int window, float sigma) const {

            // Compute the weights for Gaussian smoothing.
            std::vector<float> weights(window * window);
            gaussianWeights(weights, window, sigma);
            
            smoothed.clear();

            T v, n, s = 0;
            T norm = 0;
            for(int r = 0; r < m_rows; ++r) {
                for(int c = 0; c < m_cols; ++c) {
                    s = 0;
                    norm = 0;
                    for(int rr = -window / 2; rr <= window / 2; ++rr) {
                        for(int cc = -window / 2; cc <= window / 2; ++cc) {
                            v = m_nodata;
                            if((r + rr) >= 0 && (r + rr) < m_rows && (c + cc) >= 0 && (c + cc) < m_cols) {
                                v = m_grid[(r + rr) * m_cols + (c + cc)];
                                if(v != m_nodata) {
                                    s += v * (n = weights[(rr + window / 2) * window + (cc + window / 2)]);
                                    norm += n;
                                }
                            }
                        }
                    }
                    if(norm > 0) {
                        smoothed.set(c, r, s / norm);
                    } else {
                        smoothed.set(c, r, smoothed.nodata());
                    }
                }
            }
        }

        /**
         * Check that the given GDAL data type corresponds to the type of this grid.
         */
        void checkType(GDALDataType type) {
            if(!(std::is_same<T, float>::value && type == GDALDataType::GDT_Float32)
                && (std::is_same<T, int>::value && type == GDALDataType::GDT_Int32)) {
                throw std::runtime_error("Raster type and template type must match. Currently float32 and int32 are accepted.");
            }
        }

        /**
         * Load a grid from a raster file. The raster type must correspond to the template type.
         */
        void save(const std::string& filename) {
            if (filename.empty())
                throw std::runtime_error("Filename must be given.");

            // Attempt to open the dataset.
            
            GDALDriverManager* dm = GetGDALDriverManager();
            GDALDriver* drv = dm->GetDriverByName("GTiff");
            char** opts = nullptr;
            const char* fn = filename.c_str();
            GDALDataset* ds = (GDALDataset *) drv->Create(fn, m_cols, m_rows, 1, m_type, opts);
            if (ds == NULL)
                throw std::runtime_error("Failed to open raster.");

            ds->SetGeoTransform(m_transform.data());
            ds->SetProjection(m_crs.c_str());

            GDALRasterBand* bnd = ds->GetRasterBand(1);
            bnd->SetNoDataValue(m_nodata);
            
            if(CE_None != bnd->RasterIO(GF_Write, 0, 0, m_cols, m_rows,
                    m_grid, m_cols, m_rows, m_type, 0, 0, nullptr)) { 
                std::cerr << "Failed to write raster.";
            }
            GDALClose(ds);

        }

        /**
         * Load a grid from a raster file. The raster type must correspond to the template type.
         */
        void load(const std::string& filename, int band) {
            if (filename.empty())
                throw std::runtime_error("Filename must be given.");

            // Attempt to open the dataset.

            GDALDataset* ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
            if (ds == NULL)
                throw std::runtime_error("Failed to open raster.");

            ds->GetGeoTransform(m_transform.data());
            if(band > ds->GetRasterCount()) {
                //std::stringstream ss;
                //ss << std::string("Invalid band ") << band << "; only " << ds->GetRasterCount() << " in raster."
                GDALClose(ds);
                throw std::runtime_error("Invalid band.");
            }

            ds->GetGeoTransform(m_transform.data());
            char** interleave = ds->GetMetadata("INTERLEAVE");
            GDALRasterBand* bnd = ds->GetRasterBand(band);
            GDALDataType type = bnd->GetRasterDataType();
            try {
                checkType(type);
            } catch(const std::runtime_error& e) {
                GDALClose(ds);
                throw e;
            }

            int cols = ds->GetRasterXSize();
            int rows = ds->GetRasterYSize();
            
            initGrid(cols, rows);

            ds->GetGeoTransform(m_transform.data());
            m_crs = ds->GetProjectionRef();
            m_nodata = ds->GetRasterBand(band)->GetNoDataValue();
            m_type = type;
            m_band = band;
            //grid.m_driver = driver;

            GDALRasterIOExtraArg arg;
            INIT_RASTERIO_EXTRA_ARG(arg);

            bnd->RasterIO(GF_Read, 0, 0, m_cols, m_rows,
                    m_grid, m_cols, m_rows, m_type, 0, 0, nullptr);
            
            GDALClose(ds);
        }

        /**
         * \brief Vectorizes the raster by consuming pixels and calling back with completed GEOS polygon objects.
         *
         * The callback is a functor which accepts an int for the ID, a GEOSGeometry* and a GEOSContextHandle_t.
         * The caller is responsible for disposing of the geometry object.
         *
         * \param callback The callback functor.
         * \param pc An option PolygonContext containing configuration information.
         */
        void polygonize(PolyMergeCallback* callback, PolyCtx* pc) {

            if(!pc)
                throw std::runtime_error("A PolygonContext is requried.");

            // Extract some grid properties.
            pc->cols = cols();
            pc->rows = rows();
            pc->resX = this->xRes();
            pc->resY = this->yRes();
            copyBounds(pc->bounds);

            // The starting corner coordinates. The bounds already respect the sign of the resolution.
            pc->startX = pc->bounds[0];
            pc->startY = pc->bounds[1];

            // "Epsilon" for snapping geometries.
            double eps = 0.0001;

            // Thread control features.
            pc->merging = true;

            // Row buffer.
            std::vector<T> buf(pc->cols);
            // Lists of geoms under construction.
            std::unordered_map<int, std::vector<GEOSGeometry*>> geomParts;
            // The list of geometries currently being built.
            std::unordered_set<int> activeIds;

            // Process raster.
            for(int r = 0; r < pc->rows; ++r) {

                read(buf, 0, r, pc->cols, 1);

                // Initialize the corner coordinates.
                double x0 = pc->startX;
                double y0 = pc->startY + r * pc->resY;
                double x1 = x0;
                double y1 = y0 + pc->resY;

                // For tracking cell values. TODO: An unsigned int is possible here: overflow.
                int v0 = buf[0];
                int v1 = -1;

                // Reset the list of IDs extant in the current row.
                activeIds.clear();

                // Note: Counts past the end to trigger writing the last cell.
                for(int c = 1; c < pc->cols; ++c) {

                    // If the current cell value differs from the previous one...
                    if(c == pc->cols - 1 || (v1 = buf[c]) != v0) {
                        // Update the right x coordinate.
                        x1 = pc->startX + c * pc->resX;
                        // If the value is a valid ID, create and the geometry and save it for writing.
                        if(v0 > 0) {
                            GEOSGeometry* geom = polyMakeGeom(pc->gctx, x0, y0, x1, y1, eps, 3);
                            geomParts[v0].push_back(geom);
                            activeIds.insert(v0);
                        }
                        // Update values for next loop.
                        v0 = v1;
                        x0 = x1;
                    }
                }

                // IDs that are in the geoms array and not in the current row are ready to be finalized.
                std::vector<int> rem;
                for(const auto& it : geomParts) {
                    if(activeIds.find(it.first) == activeIds.end()) {
                        pc->geomBuf.push_back(std::make_pair(it.first, std::move(geomParts[it.first])));
                        rem.push_back(it.first);
                    }
                }
                for(int i : rem)
                    geomParts.erase(i);

            }

            // Finalize all remaining geometries.
            for(const auto& it : geomParts) {
                pc->geomBuf.push_back(std::make_pair(it.first, std::move(geomParts[it.first])));
            }
            geomParts.clear();

            // Start merge. TODO: Threading starts here.
            polyMerge(callback, pc);

            // Let the threads shut down when they run out of geometries.
            pc->merging = false;

        }        
    };


}
}

#endif // __GRID_HPP__
