#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <iostream>
#include <string>
#include <concepts>
#include <inttypes.h>
#include <sys/mman.h>

#include <geos_c.h>

#include <gdal/gdal_priv.h>
#include <gdal/ogr_spatialref.h>
#include <gdal/ogr_geometry.h>
#include <gdal/ogr_feature.h>
#include <gdal/ogrsf_frmts.h>


#define GRID_MMAP_THRESHOLD 1000000 // Change to mmap when the array size is larger than this.


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
        double m_transform[6];                  // The GDAL transform.
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
                    munmap(m_grid, sizeof(T) * m_rows * m_cols);
                } else {
                    free(m_grid);
                }
            }
        }

        /**
         * Initialize the memory for the grid. Free previously allocated memory 
         * if needed. If the size of the data segment is larger than GRID_MMAP_THRESHOLD,
         * memory is mapped, otherwise allocated.
         */
        void initGrid(int cols, int rows) {
            freeGrid();
            // Check if the size threshold is exceeded. If so, use mmap.
            size_t size = sizeof(T) * cols * rows;
            m_mapped = size > GRID_MMAP_THRESHOLD;
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
                m_grid = (T*) malloc(size);
                if(m_grid == nullptr)
                    throw "Failed to allocate data for grid.";
            }
            m_cols = cols;
            m_rows = rows;
        }

        /**
         * Copy this grid's affine transform to another array. The other array
         * must be initialized with 6 elements.
         */
        void copyTransform(double* trans) const {
            for(int i = 0; i < 6; ++i)
                trans[i] = m_transform[i];
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
                m_nodata(-9999),
                m_type(GDALDataType::GDT_Float32) {
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
            m_crs = other.m_crs;
            m_nodata = other.m_nodata;
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
            for(int r = row - h / 2; r <= row + h / 2 + 1; ++r) {
                for(int c = col - w / 2; c <= col + w / 2 + 1; ++c) {
                    float v = m_nodata;
                    if(r >= 0 && r < m_rows && c >= 0 && c < m_cols)
                        v = *(m_grid + (r * m_cols + c));
                    tile[(r - row + h / 2) * w + (c - col + w / 2)] = v;
                }
            }
        }

        /**
         * Reset all pixels to nodata.
         */
        void clear() {
            for(size_t i = 0; i < m_rows * m_cols; ++i)
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
                        smoothed.set(c, r, m_nodata);
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

            ds->SetGeoTransform(m_transform);
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

            if(band > ds->GetRasterCount()) {
                //std::stringstream ss;
                //ss << std::string("Invalid band ") << band << "; only " << ds->GetRasterCount() << " in raster."
                GDALClose(ds);
                throw std::runtime_error("Invalid band.");
            }

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

            ds->GetGeoTransform(m_transform);
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
    };


}
}

#endif // _GRID_HPP_
