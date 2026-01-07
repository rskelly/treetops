#ifndef _GRID_HPP_
#define _GRID_HPP_

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

    template <class T>
    class Grid {
    private:
        Tile<T> m_tile;                 // A Tile instance to return windows of data.
        T* m_grid;                      // The main storage location. May be in-memory or mapped. Mapping is used when the grid is larger than GRID_MMAP_THRESHOLD.
        bool m_mapped;                  // True if mapped.
        int m_cols;                // The number of columns.
        int m_rows;                // The number of rows.
        int m_band;                     // Raster band. Starts with 1.
        GDALDataType m_type;
        double m_transform[6];           // The GDAL transform.
        std::string m_crs;              // The CRS as a WKT string.
        std::vector<std::string> m_bandMeta;
        float m_nodata;
        std::string m_driver;

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
            m_mapped = sizeof(T) * cols * rows > GRID_MMAP_THRESHOLD;
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
                m_grid = (T*) malloc(sizeof(T) * cols * rows);
                if(!m_grid)
                    throw "Failed to allocate data for grid.";
            }
        }

        /**
         * Convert a map unit to a grid unit.
         */
        float map2Cell(float v) {
            return 0;
        }

        /**
         * Convert a cell unit to a map unit.
         */
        float cell2Map(float v) {
            return 0;
        }

    public:

        /**
         * Initialize a grid of the given size.
         */
        Grid(int cols, int rows) :  
                m_grid((T*) nullptr),
                m_mapped(false),
                m_cols(cols),
                m_rows(rows) {
            initGrid(cols, rows);
        }

        ~Grid() {
            freeGrid();
        }

        /**
         * Retrieve a square tile of the given size centered on the given coordinate. Tile
         * size must be odd, and will be incremented if it is even. The tile is owned by 
         * the grid and cannot be updated or destroyed. Cells outside the bounds of the grid 
         * are set to zero.
         */
        const Tile<T>& read(int col, int row, int side) {
            if(side % 2 == 0)
                ++side;
            m_tile.resize(side, side);
            for(int r = row, j = 0; r < row + side + 1; ++r, ++j) {
                for(int c = col, i = 0; c < col + side + 1; ++c, ++i) {
                    T v = 0;
                    if(r >= 0 && c >= 0 && r < m_rows && c < m_cols)
                        v = m_grid[r * m_cols + c];
                    m_tile.set(i, j, 0);
                }
            }
            return m_tile;
        }

        /**
         * Clone the current grid into the given instance.
         */
        void clone(Grid<T>& grid) {

        }

        /**
         * Smooth the grid in-place using the given parameters to the Gaussian kernel.
         */
        void smooth(int window, float sigma) {

        }

        static void checkType(GDALDataType type) {
            if(!(std::is_same<T, float>::value && type == GDALDataType::GDT_Float32)
                && (std::is_same<T, int>::value && type == GDALDataType::GDT_Int32)) {
                throw std::runtime_error("Raster type and template type must match. Currently float32 and int32 are accepted.");
            }
        }

        /**
         * Load a grid from a raster file. The raster type must correspond to the template type.
         */
        static Grid<T> load(const std::string& filename, int band) {
            if (filename.empty())
                throw std::runtime_error("Filename must be given.");

            // Attempt to open the dataset.

            GDALDataset* ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_Update);
            if (ds == NULL)
                throw std::runtime_error("Failed to open raster.");

            if(band > ds->GetRasterCount()) {
                //std::stringstream ss;
                //ss << std::string("Invalid band ") << band << "; only " << ds->GetRasterCount() << " in raster."
                GDALClose(ds);
                throw std::runtime_error("Invalid band.");
            }

            GDALDriver *drv = ds->GetDriver();
            if(drv == NULL) {
                GDALClose(ds);
                throw std::runtime_error("Failed to retrieve driver.");
            }

            char** interleave = ds->GetMetadata("INTERLEAVE");
            std::string driver = drv->GetDescription();
            GDALDataType type = ds->GetRasterBand(band)->GetRasterDataType();
            try {
                Grid<T>::checkType(type);
            } catch(const std::runtime_error& e) {
                GDALClose(ds);
                throw e;
            }

            int cols = ds->GetRasterXSize();
            int rows = ds->GetRasterYSize();
            
            Grid<T> grid(cols, rows);
            ds->GetGeoTransform(grid.m_transform);
            grid.m_crs = ds->GetProjectionRef();
            grid.m_nodata = ds->GetRasterBand(band)->GetNoDataValue();
            grid.m_type = type;
            grid.m_driver = driver;
            grid.m_band = band;

            GDALRasterBand* bnd = ds->GetRasterBand(band);
            GDALRasterIOExtraArg arg;
            INIT_RASTERIO_EXTRA_ARG(arg);
            //struct gdalprg prg;
            //prg.p = 0;
            //arg.pfnProgress = gdalProgress;
            //arg.pProgressData = &prg;

            if(CE_None != bnd->RasterIO(GF_Read, 0, 0, grid.m_cols, grid.m_rows,
                    grid.m_grid, grid.m_cols, grid.m_rows, grid.m_type, 0, 0, nullptr)) { //&arg)) {
                // If the load was deliberately canceled, don't raise an error.
                /*
                if(Monitor::get().canceled()) {
                    throw std::runtime_error("Failed to copy raster row.");
                } else {
                    Monitor::get().status(0.0f, "Load canceled.");
                }
                */
            }
            GDALClose(ds);
            GDALDestroy();

            return grid;

        }
    };


}
}

#endif // _GRID_HPP_
