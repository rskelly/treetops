#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <string>
#include <inttypes.h>
#include <sys/mman.h>

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
        uint32_t m_size;    // The length of one size of the tile, in cells.

        /**
         * Check for initialization and coordinate validity.
         */
        void check(uint32_t col, uint32_t row) {
            if(!m_data)
                throw "This tile is not initialized.";
            if(col < 0 || row < 0 || col >= m_size || row >= m_size)
                throw "Coordinate is out of bounds.";
        }

    public:

        /**
         * Construct. resize must be called before use or an exception is thrown.
         */
        Tile() {}

        /**
         * Construct and initialize with the given size.
         */
        Tile(uint32_t size) {
            resize(size);
        }

        /**
         * (Re)size the grid by allocating memory of the right size.
         */
        void resize(uint32_t size) {
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
        T get(uint32_t col, uint32_t row) {
            check(col, row);
            return m_data[row * m_size + col];
        }

        /**
         * Set the value at a coordinate.
         */
        void set(T value, uint32_t col, uint32_t row) {
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
        uint32_t m_cols;                // The number of columns.
        uint32_t m_rows;                // The number of rows.
        uint8_t m_gdal_type;            // The GDAL data type.
        float m_gdal_transform[6];      // The GDAL transform.
        std::string m_crs;              // The CRS as a WKT string.
        
        /**
         * Free or unmap the allocated memory.
         */
        void freeGrid() {
            if(m_grid) {
                if(m_mapped) {
                    munmap(m_grid);
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
        void initGrid(uint32_t cols, uint32_t rows) {
            freeGrid();
            // Check if the size threshold is exceeded. If so, use mmap.
            m_mapped = sizeof(T) * cols * rows > GRID_MMAP_THRESHOLD;
            if(m_mapped) {
                // Map the data segment.
                m_grid = mmap(
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
                m_grid = malloc(sizeof(T) * cols * rows);
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
        Grid<T>(uint32_t cols, uint32_t rows) :  
                m_grid((T*) nullptr),
                m_mapped(false),
                m_cols(cols),
                m_rows(rows) {
            initGrid(cols, rows);
        }

        /**
         * Retrieve a square tile of the given size centered on the given coordinate. Tile
         * size must be odd, and will be incremented if it is even. The tile is owned by 
         * the grid and cannot be updated or destroyed. Cells outside the bounds of the grid 
         * are set to zero.
         */
        const Tile<T>& read(uint32_t col, uint32_t row, uint32_t side) {
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
        void smooth(Grid<T>& grid, float std, float radius) {

        }

        /**
         * Load a float grid from a raster file. Corresponds to
         * GDAL type Float32.
         */
        static Grid<float> loadAsFloat(const std::string& path) {
            return Grid<float>(0, 0);
        }

        /**
         * Load an integer grid from a raster file. Corresponds to 
         * GDAL type Int32.
         */
        static Grid<int> loadAsInt(const std::string& path) {
            return Grid<int>(0, 0);
        }
    };


}
}

#endif // _GRID_HPP_
