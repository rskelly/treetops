#include <queue>
#include <string>
#include <fstream>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

#include <gdal_alg.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>
#include <cpl_string.h>

#include "geotools.hpp"
#include "util.hpp"
#include "raster.hpp"

using namespace geotools::util;
using namespace geotools::raster;

bool _cancel = false;

// Implementations for Cell

Cell::Cell(uint16_t col, uint16_t row) :
	col(col), row(row) {
}

// Implementations for TargetOperator (for flood fill)

template<class T>
FillOperator<T>::~FillOperator() {}

template<class T>
TargetOperator<T>::TargetOperator(T match) :
		m_match(match) {
}

template<class T>
bool TargetOperator<T>::fill(T value) const {
	return value == m_match;
}

template<class T>
TargetOperator<T>::~TargetOperator() {}

template<class T>
Grid<T>::Grid() :
		m_min(0), m_max(0), m_mean(0), m_stddev(0), m_variance(0), m_sum(0), m_count(
				0), m_stats(false) {
}

// Implementations forthe Grid class

template<class T>
void Grid<T>::gaussianWeights(double *weights, uint16_t size, double sigma) {
	// If size is an even number, bump it up.
	if (size % 2 == 0) {
		++size;
		g_warn(
				"Gaussian kernel size must be an odd number >=3. Bumping up to "
						<< size);
	}
	for (int32_t r = 0; r < size; ++r) {
		for (int32_t c = 0; c < size; ++c) {
			int32_t x = size / 2 - c;
			int32_t y = size / 2 - r;
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma))
					* pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
		}
	}
}

template<class T>
void Grid<T>::computeStats() {
	uint32_t i;
	for (i = 0; i < size(); ++i) {
		if (!isNoData(i)) {
			m_min = m_max = get(i);
			break;
		}
	}
	m_sum = 0;
	m_count = 0;
	T m = 0;
	T s = 0;
	int32_t k = 1;
	// Welford's method for variance.
	// i has the index of the first non-nodata element.
	for (; i < size(); ++i) {
		T v = get(i);
		if (v != nodata()) {
			T oldm = m;
			m = m + (v - m) / k;
			s = s + (v - m) * (v - oldm);
			m_sum += v;
			m_min = g_min(m_min, v);
			m_max = g_max(m_max, v);
			++m_count;
			++k;
		}
	}
	m_mean = m_sum / m_count;
	m_variance = s / m_count;
	m_stddev = (T) std::sqrt(m_variance);
	m_stats = true;
}

template<class T>
void Grid<T>::normalize() {
	double sum = 0.0;
	for (uint32_t i = 0; i < size(); ++i) {
		double v = (double) get(i);
		if (v != nodata() && !std::isnan(v) && v < G_DBL_MAX_POS)
			sum += v;
	}
	double mean = sum / size();
	sum = 0.0;
	for (uint32_t i = 0; i < size(); ++i) {
		double v = (double) get(i);
		if (v != nodata() && !std::isnan(v) && v < G_DBL_MAX_POS)
			sum += std::pow(v - mean, 2.0);
	}
	double stdDev = std::sqrt(sum);
	for (uint32_t i = 0; i < size(); ++i) {
		double v = (double) get(i);
		if (v != nodata() && !std::isnan(v) && v < G_DBL_MAX_POS) {
			set(i, (const T) ((v - mean) / stdDev));
		} else {
			set(i, nodata());
		}
	}
}

template<class T>
void Grid<T>::logNormalize() {
	T n = min();
	T x = max();
	double e = std::exp(1.0) - 1.0;
	for(uint32_t i = 0; i < size(); ++i)
		set(i, (T) std::log(1.0 + e * (get(i) - n) / (x - n)));
}

template<class T>
T Grid<T>::max() {
	if (!m_stats)
		computeStats();
	return m_max;
}

template<class T>
T Grid<T>::min() {
	if (!m_stats)
		computeStats();
	return m_min;
}

template<class T>
T Grid<T>::mean() {
	if (!m_stats)
		computeStats();
	return m_mean;
}

template<class T>
T Grid<T>::stddev() {
	if (!m_stats)
		computeStats();
	return m_stddev;
}

template<class T>
T Grid<T>::variance() {
	if (!m_stats)
		computeStats();
	return m_variance;
}

template<class T>
void Grid<T>::floodFill(int32_t col, int32_t row, T target, T fill, bool d8,
		uint16_t *outminc, uint16_t *outminr,
		uint16_t *outmaxc, uint16_t *outmaxr,
		uint32_t *outarea) {
	TargetOperator<T> op(target);
	return floodFill(col, row, op, *this, fill, d8, outminc, outminr,
			outmaxc, outmaxr, outarea);
}

template<class T>
void Grid<T>::floodFill(int32_t col, int32_t row,
		FillOperator<T> &op, T fill, bool d8,
		uint16_t *outminc, uint16_t *outminr,
		uint16_t *outmaxc, uint16_t *outmaxr,
		uint32_t *outarea) {
	return floodFill(col, row, op, *this, fill, d8, outminc, outminr,
			outmaxc, outmaxr, outarea);
}

template<class T>
void Grid<T>::voidFillIDW(double radius, uint32_t count, double exp) {

	if (radius <= 0.0)
		throw std::invalid_argument("Radius must be larger than 0.");

	if (count <= 0)
		throw std::invalid_argument("Count must be larger than 0.");

	if (exp <= 0.0)
		throw std::invalid_argument("Exponent must be larger than 0.");

	MemRaster<T> tmp(cols(), rows());
	tmp.setNodata(nodata());
	tmp.fill(nodata());

	for (int32_t r = 0; r < rows(); ++r) {
		for (int32_t c = 0; c < cols(); ++c) {

			if (get(c, r) != nodata())
				continue;

			double rad = radius;
			bool found = false;

			do {

				double d = g_sq(rad);
				double a = 0.0;
				double b = 0.0;
				uint32_t cnt = 0;

				for (int32_t r0 = (int32_t) g_max(0, r - rad);
						r0 < (int32_t) g_min(rows(), r + rad + 1); ++r0) {
					for (int32_t c0 = (int32_t) g_max(0, c - rad);
							c0 < (int32_t) g_min(cols(), c + rad + 1); ++c0) {
						double d0 = g_sq((double) c0 - c)
								+ g_sq((double) r0 - r);
						if (d0 <= d && get(c0, r0) != nodata()) {
							double dp = 1.0 / std::pow(d0, exp);
							a += dp * get(c0, r0);
							b += dp;
							++cnt;
						}
					}
				}

				if (cnt >= count) {
					tmp.set(c, r, (T) (a / b));
					found = true;
					break;
				}

				rad += 1.0;

			} while (rad < g_min(cols(), rows()));

			if (!found)
				g_warn("Pixel not filled at " << c << "," << r
						<< ". Consider larger radius or smaller count.");
		}
	}

	writeBlock(tmp);
}

template<class T>
void Grid<T>::smooth(Grid<T> &smoothed, double sigma, uint16_t size,
		Callbacks *status, bool *cancel) {

	if (!cancel)
		cancel = &_cancel;

	if (status) {
		status->stepCallback(0.01f);
		status->statusCallback("Preparing...");
	}

	if (sigma <= 0)
		g_argerr("Sigma must be > 0.");
	if (size < 3)
		g_argerr("Kernel size must be 3 or larger.");
	if (size % 2 == 0) {
		g_warn("Kernel size must be odd. Rounding up.");
		size++;
	}

	double nd = (double) nodata();

	MemRaster<double> weights(size, size);
	gaussianWeights(weights.grid(), size, sigma);

	uint16_t bufSize = g_max(size, 1024);

	MemRaster<T> buf(cols(), bufSize + size, false);
	buf.setNodata((T) nd);
	buf.fill((T) nd);

	MemRaster<T> smooth(cols(), bufSize + size, false);
	smooth.setNodata((T) nd);
	smooth.fill((T) nd);

	if (status)
		status->stepCallback(0.02f);

	uint16_t curRow = 0;
	for (int32_t b = 0; b < rows() - size; b += bufSize) {
		if (*cancel) continue;

		if(status)
			status->statusCallback("Reading...");

		buf.fill((T) nd);
		smooth.fill((T) nd);

		uint16_t readOffset = b > 0 ? b - size / 2 : 0;
		uint16_t writeOffset = b > 0 ? 0 : size / 2;
		readBlock(0, readOffset, buf, 0, writeOffset);

		if(status)
			status->statusCallback("Processing...");

		for (int32_t r = 0; r < buf.rows() - size; ++r) {
			if (*cancel) continue;
			for (int32_t c = 0; c < buf.cols() - size; ++c) {
				double v, t = 0.0;
				bool foundNodata = false;
				for (int32_t gr = 0; !foundNodata && gr < size; ++gr) {
					for (int32_t gc = 0; !foundNodata && gc < size; ++gc) {
						v = (double) buf.get(c + gc, r + gr);
						if (v == nd) {
							foundNodata = true;
						} else {
							t += weights.get(gr * size + gc) * v;
						}
					}
				}
				if (!foundNodata)
					smooth.set(c + size / 2, r + size / 2, (T) t);
			}
		}

		if (status) {
			curRow += bufSize;
			status->stepCallback(0.2f + (float) b / rows() * 0.96f);
			status->statusCallback("Writing...");
		}

		writeOffset = b > 0 ? b - size / 2 : b;
		readOffset = size / 2;
		smoothed.writeBlock(0, b, smooth, 0, size / 2);
	}

	if (status) {
		status->stepCallback(1.0);
		status->statusCallback("Closing...");
	}
}

template<class T>
Grid<T>::~Grid() {
}

// Implementations for MemRaster

template<class T>
void MemRaster<T>::checkInit() const {
	if (m_grid == nullptr)
		g_runerr("This instance has not been initialized.");
}

template<class T>
MemRaster<T>::MemRaster() :
		m_grid(nullptr), m_cols(-1), m_rows(-1), m_item_dealloc(nullptr), m_nodata(
				0), m_mmapped(false), m_size(0) {
}

template<class T>
MemRaster<T>::MemRaster(uint16_t cols, uint16_t rows, bool mapped) {
	init(cols, rows, mapped);
}

template<class T>
MemRaster<T>::~MemRaster() {
	if (m_item_dealloc) {
		for (uint32_t i = 0; i < (uint32_t) m_cols * m_rows; ++i)
			m_item_dealloc(m_grid[i]);
	}
	freeMem();
}

template<class T>
void MemRaster<T>::setDeallocator(void (*item_dealloc)(T)) {
	m_item_dealloc = item_dealloc;
}

template<class T>
T *MemRaster<T>::grid() {
	return m_grid;
}

template<class T>
bool MemRaster<T>::hasGrid() const {
	return true;
}

template<class T>
uint16_t MemRaster<T>::rows() const {
	return m_rows;
}

template<class T>
uint16_t MemRaster<T>::cols() const {
	return m_cols;
}

template<class T>
uint32_t MemRaster<T>::size() const {
	return (uint32_t) m_rows * m_cols;
}

template<class T>
void MemRaster<T>::freeMem() {
	if (m_grid) {
		if (m_mmapped) {
			m_mappedFile.release();
		} else {
			free(m_grid);
		}
	}
}

template<class T>
void MemRaster<T>::init(uint16_t cols, uint16_t rows, bool mapped) {

	m_grid = nullptr;
	m_cols = 0;
	m_rows = 0;
	m_item_dealloc = nullptr;
	m_nodata = 0;
	m_mmapped = false;
	m_size = 0;

	if (cols != m_cols || rows != m_rows) {
		freeMem();
		m_cols = cols;
		m_rows = rows;
		m_mmapped = mapped;
		m_grid = nullptr;
		m_size = sizeof(T) * cols * rows;
		m_mappedFile.release();
		if (mapped) {
			const std::string filename = Util::tmpFile("/tmp");
			m_mappedFile = Util::mapFile(filename, m_size);
			m_grid = (T *) m_mappedFile->data();
		} else {
			m_grid = (T *) malloc(m_size);
		}
		if (!m_grid)
			g_runerr("Failed to allocate memory for MemRaster.");
	}
}

template<class T>
void MemRaster<T>::fill(const T value) {
	checkInit();
	for(uint64_t i = 0; i < size(); ++i)
		m_grid[i] = value;
}

template<class T>
T MemRaster<T>::get(uint32_t idx) {
	checkInit();
	if (idx >= size())
		g_argerr("Index out of bounds: " << idx << "; size: " << size());
	return m_grid[idx];
}

template<class T>
T MemRaster<T>::get(int32_t col, int32_t row) {
	uint32_t idx = (uint32_t) row * m_cols + col;
	return get(idx);
}

template<class T>
bool MemRaster<T>::isNoData(int32_t col, int32_t row) {
	return get(col, row) == m_nodata;
}

template<class T>
bool MemRaster<T>::isNoData(uint32_t idx) {
	return get(idx) == m_nodata;
}

template<class T>
void MemRaster<T>::set(int32_t col, int32_t row, const T value) {
	uint32_t idx = (uint32_t) row * m_cols + col;
	set(idx, value);
}

template<class T>
void MemRaster<T>::set(uint32_t idx, const T value) {
	checkInit();
	if (idx >= size())
		g_argerr(
				"Index out of bounds: " << idx << "; size: " << size()
						<< "; value: " << value << "; col: " << (idx % m_cols)
						<< "; row: " << (idx / m_cols));
	m_grid[idx] = value;
}

template<class T>
bool MemRaster<T>::has(int32_t col, int32_t row) const {
	return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
}

template<class T>
bool MemRaster<T>::has(uint32_t idx) const {
	return idx < (uint32_t) m_cols * m_rows;
}

template<class T>
bool MemRaster<T>::isSquare() const {
	return cols() == rows();
}

template<class T>
void MemRaster<T>::toMatrix(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
	for (int32_t r = 1; r < rows(); ++r) {
		for (int32_t c = 0; c < cols(); ++c)
			mtx(r, c) = get(c, r);
	}
}

template<class T>
void MemRaster<T>::fromMatrix(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx) {
	for (int32_t r = 1; r < rows(); ++r) {
		for (int32_t c = 0; c < cols(); ++c)
			set(c, r, mtx(r, c));
	}
}

template<class T>
T MemRaster<T>::nodata() const {
	return m_nodata;
}

template<class T>
void MemRaster<T>::setNodata(T nodata) {
	m_nodata = nodata;
}

template<class T>
void MemRaster<T>::readBlock(int32_t col, int32_t row, Grid<T> &block,
		int32_t dstCol, int32_t dstRow,
		int32_t xcols, int32_t xrows) {
	if (&block == this)
		g_argerr("Recursive call to readBlock.");
	if (dstCol < 0 || dstRow < 0 || dstCol >= block.cols() || dstRow >= block.rows())
		g_argerr("Invalid destination column or row: row: " << dstRow
			<< "; col: " << dstCol << "; block: " << block.rows() << "," << block.rows());
	if(col == 0 && row == 0 && block.cols() == cols() && block.rows() == rows() && block.hasGrid() && hasGrid()) {
		// If these are equivalent blocks and they have a grid, copy the memory.
		std::memcpy(block.grid(), grid(), size());
	} else {
		uint16_t cols = g_min(m_cols - col, block.cols() - dstCol);
		uint16_t rows = g_min(m_rows - row, block.rows() - dstRow);
		if (block.hasGrid()) {
			for (int32_t r = 0; r < rows; ++r) {
				std::memcpy(
					(block.grid() + (dstRow + r) * block.cols() + dstCol),
					(m_grid + (row + r) * m_cols + col),
					cols * sizeof(T)
				);
			}
		} else {
			for (int32_t r = 0; r < rows; ++r) {
				for (int32_t c = 0; c < cols; ++c)
					block.set(c + dstCol, r + dstRow, get(c + col, r + row));
			}
		}
	}
}

template<class T>
void MemRaster<T>::writeBlock(int32_t col, int32_t row, Grid<T> &block,
		int32_t srcCol, int32_t srcRow,
		int32_t xcols, int32_t xrows) {
	if (&block == this)
		g_argerr("Recursive call to writeBlock.");
	if (srcCol < 0 || srcRow < 0 || srcCol >= block.cols() || srcRow >= block.rows())
		g_argerr("Invalid source column or row: row: " << srcRow << "; col: "
				<< srcCol << "; block: " << block.rows() << ","	<< block.rows());
	if(col == 0 && row == 0 && block.cols() == cols() && block.rows() == rows() && block.hasGrid() && hasGrid()) {
		// If these are equivalent blocks and they have a grid, copy the memory.
		std::memcpy(grid(), block.grid(), size());
	} else {
		uint16_t cols = g_min(m_cols - col, block.cols() - srcCol);
		uint16_t rows = g_min(m_rows - row, block.rows() - srcRow);
		if (xcols > 0)
			cols = g_min(xcols, cols);
		if (xrows > 0)
			rows = g_min(xrows, rows);
		if (block.hasGrid()) {
			for (int32_t r = 0; r < rows; ++r) {
				std::memcpy(
						(m_grid + (r + row) * m_cols + col),
						(block.grid() + (r + srcRow) * block.cols() + srcCol),
						cols * sizeof(T)
				);
			}
		} else {
			for (int32_t r = 0; r < rows; ++r) {
				for (int32_t c = 0; c < cols; ++c)
					set(c + col, r + row, block.get(c + srcCol, r + srcRow));
			}
		}
	}
}

template<class T>
void MemRaster<T>::writeBlock(Grid<T> &block) {
	writeBlock(0, 0, block, 0, 0, block.cols(), block.rows());
}

template<class T>
void MemRaster<T>::readBlock(Grid<T> &block) {
	readBlock(0, 0, block, 0, 0, block.cols(), block.rows());
}

// Implementations for BlockCache

template <class T>
Block<T>::Block(uint64_t idx, uint16_t band,
		uint16_t col, uint16_t row,
		uint32_t size, uint64_t time) :
	idx(idx),
	band(band), col(col), row(row),
	size(size),
	time(time),
	dirty(false),
	data(nullptr) {
	data = (T*) malloc(sizeof(T) * size);
}

template <class T>
Block<T>::~Block() {
	free(data);
}

template<class T>
void BlockCache<T>::flushBlock(uint64_t idx) {
	if(m_blocks.find(idx) != m_blocks.end())
		flushBlock(m_blocks[idx]);
}

template<class T>
void BlockCache<T>::flushBlock(Block<T> *blk) {
	if (blk && m_ds->GetAccess() == GA_Update) {
		#pragma omp critical(__gdal_io)
		{
			g_debug("Flushing");
			if (m_ds->GetRasterBand(blk->band)->WriteBlock(blk->col, blk->row,
					blk->data) != CE_None)
				g_runerr("Failed to flush block.");
		}
		blk->dirty = false;
	}
}

template<class T>
uint64_t BlockCache<T>::toIdx(uint16_t band, uint16_t col, uint16_t row) const {
	return ((uint64_t) band << 32) | ((uint64_t) (col / m_bw) << 16) | (row / m_bh);
}

template<class T>
uint16_t BlockCache<T>::toCol(uint64_t idx) const {
	return (uint16_t) ((idx >> 16) & 0xffff) * m_bw;
}

template<class T>
uint16_t BlockCache<T>::toRow(uint64_t idx) const {
	return (uint16_t) (idx & 0xffff) * m_bh;
}

template<class T>
uint16_t BlockCache<T>::toBand(uint64_t idx) const {
	return ((idx >> 32) & 0xffff);
}

template<class T>
Block<T>* BlockCache<T>::freeOldest() {
	auto it = m_time_idx.rbegin();
	uint64_t time = it->first;
	uint64_t idx = it->second;
	Block<T> *blk = m_blocks[idx];
	if(!blk)
		g_runerr("Found a null block.");
	if(blk->dirty)
		flushBlock(blk);
	m_blocks.erase(idx);
	m_time_idx.erase(time);
	return blk;
}

template<class T>
Block<T>* BlockCache<T>::freeOne() {
	Block<T> *blk = nullptr;
	while (m_blocks.size() >= m_size) {
		if (blk)
			delete blk;
		blk = freeOldest();
	}
	return blk;
}

template<class T>
BlockCache<T>::BlockCache() :
	m_ds(nullptr),
	m_size(0), m_time(0),
	m_bw(0), m_bh(0) {
}

template<class T>
uint16_t BlockCache<T>::blockWidth() const {
	return m_bw;
}

template<class T>
uint16_t BlockCache<T>::blockHeight() const {
	return m_bh;
}

template<class T>
uint32_t BlockCache<T>::toBlockIdx(uint16_t col, uint16_t row) const {
	return (row % m_bh) * m_bw + (col % m_bw);
}

template<class T>
void BlockCache<T>::setDataset(GDALDataset *ds) {
	m_ds = ds;
	m_ds->GetRasterBand(1)->GetBlockSize(&m_bw, &m_bh);
}

template<class T>
bool BlockCache<T>::hasBlock(uint16_t band, uint16_t col, uint16_t row) const {
	return hasBlock(toIdx(band, col, row));
}

template<class T>
bool BlockCache<T>::hasBlock(uint64_t idx) const {
	bool has = false;
	#pragma omp critical(__blockcache)
	{
		has = m_blocks.find(idx) != m_blocks.end();
	}
	return has;
}

template<class T>
void BlockCache<T>::setSize(uint32_t size) {
	#pragma omp critical(__blockcache)
	{
		while (m_blocks.size() > size)
			freeOne();
	}
	m_size = size;
}

template<class T>
uint32_t BlockCache<T>::getSize() {
	return m_size;
}

template<class T>
Block<T>* BlockCache<T>::getBlock(uint16_t band, uint16_t col, uint16_t row, bool forWrite) {
	uint64_t idx = toIdx(band, col, row);
	Block<T> *blk = nullptr;
	if (m_blocks.find(idx) == m_blocks.end()) {
		blk = freeOne();
		if (!blk) {
			blk = new Block<T>(idx, band, col, row, sizeof(T) * m_bw * m_bh, ++m_time);
		} else {
			blk->idx = idx;
			blk->band = band;
			blk->col = col;
			blk->row = row;
			blk->time = ++m_time;
		}
		#pragma omp critical(__gdal_io)
		{
			if (m_ds->GetRasterBand(band)->ReadBlock(col / m_bw, row / m_bh, (void *) blk->data) != CE_None)
				g_runerr("Failed to read block.");
		}
		m_blocks[idx] = blk;
		m_time_idx[blk->time] = idx;
	} else {
		blk = m_blocks[idx];
		if(m_time_idx.find(blk->time) != m_time_idx.end())
			m_time_idx.erase(blk->time);
		blk->time = ++m_time;
		m_time_idx[blk->time] = idx;
	}
	return blk;
}

template<class T>
T BlockCache<T>::get(uint16_t band, int32_t col, int32_t row) {
	T v;
	#pragma omp critical(__blockcache)
	{
		Block<T> *blk = getBlock(band, col, row, false);
		v = blk->data[toBlockIdx(col, row)];
	}
	return v;
}

template<class T>
void BlockCache<T>::set(uint16_t band, int32_t col, int32_t row, T value) {
	#pragma omp critical(__blockcache)
	{
		Block<T> *blk = getBlock(band, col, row, true);
		blk->data[toBlockIdx(col, row)] = value;
	}
}

template<class T>
void BlockCache<T>::flush() {
	#pragma omp critical(__blockcache)
	{
		for (const auto &it : m_blocks)
			flushBlock(it.second);
	}
}

template<class T>
void BlockCache<T>::close() {
	flush();
	#pragma omp critical(__blockcache)
	{
		for(auto &it : m_blocks)
			if(it.second) delete it.second;
		m_blocks.clear();
	}
}

template<class T>
BlockCache<T>::~BlockCache() {
	close();
}

// Implementations for Raster

template<class T>
GDALDataType Raster<T>::getType(uint64_t v) {
	(void) v;
	g_runerr("Raster with 64 bit integral type requested. Not implemented.");
}

template<class T>
GDALDataType Raster<T>::getType(int64_t v) {
	(void) v;
	g_runerr("Raster with 64 bit integral type requested. Not implemented.");
}

template<class T>
GDALDataType Raster<T>::getType(double v) {
	(void) v;
	return GDT_Float64;
}

template<class T>
GDALDataType Raster<T>::getType(float v) {
	(void) v;
	return GDT_Float32;
}

template<class T>
GDALDataType Raster<T>::getType(uint32_t v) {
	(void) v;
	return GDT_UInt32;
}

template<class T>
GDALDataType Raster<T>::getType(int32_t v) {
	(void) v;
	return GDT_Int32;
}

template<class T>
GDALDataType Raster<T>::getType(uint16_t v) {
	(void) v;
	return GDT_UInt16;
}

template<class T>
GDALDataType Raster<T>::getType(int16_t v) {
	(void) v;
	return GDT_Int16;
}

template<class T>
GDALDataType Raster<T>::getType(uint8_t v) {
	(void) v;
	return GDT_Byte;
}

template<class T>
GDALDataType Raster<T>::getType(int8_t v) {
	(void) v;
	return GDT_Byte;
}

template<class T>
GDALDataType Raster<T>::getType() {
	return getType((T) 0);
}

template<class T>
T Raster<T>::getDefaultNodata() {
	return 0;
}

template<class T>
Raster<T>::Raster() :
	m_cols(-1), m_rows(-1),
	m_bands(1),
	m_writable(false),
	m_ds(nullptr),
	m_type(getType()),
	m_inited(false) {
}

template<class T>
Raster<T>::Raster(const std::string &filename, uint16_t bands, const Raster<T> &tpl) :
	m_cols(-1), m_rows(-1),
	m_bands(1),
	m_writable(false),
	m_ds(nullptr),
	m_type(getType()),
	m_inited(false) {
	init(filename, bands, tpl);
}

template<class T>
Raster<T>::Raster(const std::string &filename, uint16_t bands, double minx,
		double miny, double maxx, double maxy, double resolutionX,
		double resolutionY, const std::string &proj) :
		m_cols(-1), m_rows(-1),
		m_bands(1),
		m_writable(false),
		m_ds(nullptr),
		m_type(getType()),
		m_inited(false) {
	init(filename, bands, minx, miny, maxx, maxy, resolutionX, resolutionY, proj);
}

template<class T>
Raster<T>::Raster(const std::string &filename, uint16_t bands, const Bounds &bounds,
		double resolutionX, double resolutionY, uint16_t crs) :
		m_cols(-1), m_rows(-1),
		m_bands(1),
		m_writable(false),
		m_ds(nullptr),
		m_type(getType()),
		m_inited(false) {
	std::string proj = epsg2ProjText(crs);
	init(filename, bands, bounds.minx(), bounds.miny(), bounds.maxx(),
			bounds.maxy(), resolutionX, resolutionY, proj);
}

template<class T>
Raster<T>::Raster(const std::string &filename, uint16_t bands, double minx,
		double miny, double maxx, double maxy, double resolutionX,
		double resolutionY, uint16_t crs) :
		m_cols(-1), m_rows(-1),
		m_bands(1),
		m_writable(false),
		m_ds(nullptr),
		m_type(getType()),
		m_inited(false) {
	std::string proj = epsg2ProjText(crs);
	init(filename, bands, minx, miny, maxx, maxy, resolutionX, resolutionY, proj);
}

template<class T>
Raster<T>::Raster(const std::string &filename, bool writable) :
		m_cols(-1), m_rows(-1),
		m_bands(1),
		m_writable(false),
		m_ds(nullptr),
		m_type(getType()),
		m_inited(false) {
	init(filename, writable);
}

template<class T>
void Raster<T>::init(const std::string &filename, uint16_t bands,
		const Bounds &bounds, double resolutionX, double resolutionY,
		const std::string &proj) {
	init(filename, bands, bounds.minx(), bounds.miny(), bounds.maxx(),
			bounds.maxy(), resolutionX, resolutionY, proj);
}

template<class T>
void Raster<T>::init(const std::string &filename, uint16_t bands, double minx,
		double miny, double maxx, double maxy, double resolutionX,
		double resolutionY, const std::string &proj) {
	if (resolutionX == 0 || resolutionY == 0)
		g_argerr("Resolution must be larger or smaller than 0.");
	if (maxx < minx)
		g_argerr("Minimum x must be smaller than or equal to maximum x.");
	if (maxy < miny)
		g_argerr("Minimum y must be smaller than or equal to maximum y.");
	if (filename.empty())
		g_argerr("Filename must be given.");

	m_filename.assign(filename);
	m_type = getType();

	// Compute columns/rows
	int32_t width = (int) std::ceil((maxx - minx) / g_abs(resolutionX));
	int32_t height = (int) std::ceil((maxy - miny) / g_abs(resolutionY));

	// Create GDAL dataset.
	char **opts = NULL;
	opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
	opts = CSLSetNameValue(opts, "BIGTIFF", "YES");
	GDALAllRegister();
	m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
			filename.c_str(), width, height, bands, m_type, opts);
	if (m_ds == nullptr)
		g_runerr("Failed to create file.");

	// Initialize geotransform.
	double trans[6] = {
		resolutionX < 0 ? maxx : minx, resolutionX, 0.0,
		resolutionY < 0 ? maxy : miny, 0.0, resolutionY
	};
	for (int i = 0; i < 6; ++i)
		m_trans[i] = trans[i];
	m_ds->SetGeoTransform(m_trans);

	// Set projection.
	if (!proj.empty())
		m_ds->SetProjection(proj.c_str());

	// Save some dataset properties.
	m_rows = m_ds->GetRasterYSize();
	m_cols = m_ds->GetRasterXSize();
	m_bands = bands;
	m_cache.setSize(100);
	m_cache.setDataset(m_ds);
	m_writable = true;
	m_inited = true;
}

template<class T>
void Raster<T>::init(const std::string &filename, bool writable) {

	if (filename.empty())
		g_argerr("Filename must be given.");

	m_filename.assign(filename);
	m_type = getType();

	// Attempt to open the dataset.
	GDALAllRegister();
	m_ds = (GDALDataset *) GDALOpen(filename.c_str(),
			writable ? GA_Update : GA_ReadOnly);
	if (m_ds == NULL)
		g_runerr("Failed to open raster.");

	// Save some raster
	m_ds->GetGeoTransform(m_trans);
	m_rows = m_ds->GetRasterYSize();
	m_cols = m_ds->GetRasterXSize();
	m_cache.setSize(100);
	m_cache.setDataset(m_ds);
	m_writable = writable;
	m_inited = true;
}

template<class T>
int32_t Raster<T>::getType(const std::string &filename) {
	GDALDataset *ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
	int32_t type = ds->GetRasterBand(1)->GetRasterDataType();
	GDALClose(ds);
	switch (type) {
	case GDT_Byte:
		return Raster<T>::BYTE;
	case GDT_UInt16:
		return Raster<T>::UINT16;
	case GDT_UInt32:
		return Raster<T>::UINT32;
	case GDT_Int16:
		return Raster<T>::INT16;
	case GDT_Int32:
		return Raster<T>::INT32;
	case GDT_Float32:
		return Raster<T>::FLOAT32;
	case GDT_Float64:
		return Raster<T>::FLOAT64;
	default:
		g_argerr("Unknown data type: " << type);
	}
}

template<class T>
void Raster<T>::setCacheSize(uint32_t size) {
	m_cache.setSize(size);
}

template<class T>
std::string Raster<T>::filename() const {
	return m_filename;
}

template<class T>
uint16_t Raster<T>::bandCount() const {
	return m_ds->GetRasterCount();
}

template<class T>
std::string Raster<T>::epsg2ProjText(uint16_t crs) const {
	OGRSpatialReference ref;
	char *wkt;
	ref.importFromEPSG(crs);
	ref.exportToWkt(&wkt);
	return std::string(wkt);
}

template<class T>
bool Raster<T>::inited() const {
	return m_inited;
}

template<class T>
void Raster<T>::fill(T value, uint16_t band) {
	m_cache.flush();
	MemRaster<T> grd(m_cache.blockWidth(), m_cache.blockHeight());
	grd.fill(value);
	#pragma omp critical(__gdal_io)
	{
		GDALRasterBand *bnd = m_ds->GetRasterBand(band);
		for (int32_t r = 0; r < rows() / grd.rows(); ++r) {
			for (int32_t c = 0; c < cols() / grd.cols(); ++c) {
				if (bnd->WriteBlock(c, r, grd.grid()) != CE_None)
					g_runerr("Fill error.");
			}
		}
	}
}

template<class T>
void Raster<T>::fill(T value) {
	fill(value, 1);
}

template<class T>
void Raster<T>::readBlock(uint16_t band, int32_t col, int32_t row, Grid<T> &grd,
		int32_t dstCol, int32_t dstRow,
		int32_t xcols, int32_t xrows) {
	if (&grd == this)
		g_runerr("Recursive call to readBlock.");
	m_cache.flush();
	uint16_t cols = g_min(m_cols - col, grd.cols() - dstCol);
	uint16_t rows = g_min(m_rows - row, grd.rows() - dstRow);
	if (xcols > 0)
		cols = g_min(xcols, cols);
	if (xrows > 0)
		rows = g_min(xrows, rows);
	if (cols < 1 || rows < 1)
		g_argerr("Zero read size.");
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	if (grd.hasGrid()) {
		if (cols != grd.cols() || rows != grd.rows()) {

			MemRaster<T> g(cols, rows);
			#pragma omp critical(__gdal_io)
			{
				if (bnd->RasterIO(GF_Read, col, row, cols, rows, g.grid(),
						cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to read from band. (1)");
				bnd->FlushCache();
			}
			grd.writeBlock(dstCol, dstRow, g);
		} else {
			#pragma omp critical(__gdal_io)
			{
				if (bnd->RasterIO(GF_Read, col, row, cols, rows, grd.grid(),
						cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to read from band. (2)");
				bnd->FlushCache();
			}
		}
	} else {
		MemRaster<T> mr(cols, rows);
		#pragma omp critical(__gdal_io)
		{
			if (bnd->RasterIO(GF_Read, col, row, cols, rows, mr.grid(), cols,
					rows, getType(), 0, 0) != CE_None)
				g_runerr("Failed to read from band. (3)");
			bnd->FlushCache();
		}
		grd.writeBlock(dstCol, dstRow, mr);
	}
}

template<class T>
void Raster<T>::readBlock(int32_t col, int32_t row, Grid<T> &grd,
		int32_t dstCol, int32_t dstRow,
		int32_t xcols, int32_t xrows) {
	readBlock(1, col, row, grd, dstCol, dstRow, xcols, xrows);
}

template<class T>
void Raster<T>::writeBlock(uint16_t band, int32_t col, int32_t row, Grid<T> &grd,
		int32_t srcCol, int32_t srcRow, int32_t xcols, int32_t xrows) {
	if (&grd == this)
		g_runerr("Recursive call to writeBlock.");
	m_cache.flush();
	uint16_t cols = g_min(m_cols - col, grd.cols() - srcCol);
	uint16_t rows = g_min(m_rows - row, grd.rows() - srcRow);
	if (xcols > 0)
		cols = g_min(xcols, cols);
	if (xrows > 0)
		rows = g_min(xrows, rows);
	if (cols < 1 || rows < 1)
		g_argerr("Zero write size.");
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	if (grd.hasGrid()) {
		if (cols != grd.cols() || rows != grd.rows()) {
			MemRaster<T> g(cols, rows);
			grd.readBlock(srcCol, srcRow, g);
			#pragma omp critical(__gdal_io)
			{
				if (bnd->RasterIO(GF_Write, col, row, cols, rows, g.grid(),
						cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to write to band. (1)");
				bnd->FlushCache();
			}
		} else {
			#pragma omp critical(__gdal_io)
			{
				if (bnd->RasterIO(GF_Write, col, row, cols, rows, grd.grid(),
						cols, rows, getType(), 0, 0) != CE_None)
					g_runerr("Failed to write to band. (2)");
				bnd->FlushCache();
			}
		}
	} else {
		MemRaster<T> mr(cols, rows);
		grd.readBlock(srcCol, srcRow, mr);
		#pragma omp critical(__gdal_io)
		{
			if (bnd->RasterIO(GF_Write, col, row, cols, rows, mr.grid(),
					cols, rows, getType(), 0, 0) != CE_None)
				g_runerr("Failed to write to band. (3)");
			bnd->FlushCache();
		}
	}
}

template<class T>
void Raster<T>::writeBlock(int32_t col, int32_t row, Grid<T> &grd,
		int32_t srcCol, int32_t srcRow, int32_t xcols, int32_t xrows) {
	writeBlock(1, col, row, grd, srcCol, srcRow, xcols, xrows);
}

template<class T>
void Raster<T>::writeBlock(uint16_t band, Grid<T> &block) {
	writeBlock(band, 0, 0, block);
}

template<class T>
void Raster<T>::writeBlock(Grid<T> &block) {
	writeBlock(1, block);
}

template<class T>
void Raster<T>::readBlock(uint16_t band, Grid<T> &block) {
	readBlock(band, 0, 0, block);
}

template<class T>
void Raster<T>::readBlock(Grid<T> &block) {
	readBlock(1, block);
}

template<class T>
double Raster<T>::resolutionX() const {
	return m_trans[1];
}

template<class T>
double Raster<T>::resolutionY() const {
	return m_trans[5];
}

template<class T>
bool Raster<T>::positiveX() const {
	return resolutionX() > 0;
}

template<class T>
bool Raster<T>::positiveY() const {
	return resolutionY() > 0;
}

template<class T>
void Raster<T>::projection(std::string &proj) const {
	proj.assign(m_ds->GetProjectionRef());
}

template<class T>
GDALDataType Raster<T>::type() const {
	return m_ds->GetRasterBand(1)->GetRasterDataType();
}

template<class T>
Bounds Raster<T>::bounds() const {
	return Bounds(minx(), miny(), maxx(), maxy());
}

template<class T>
double Raster<T>::minx() const {
	return toX(0);
}

template<class T>
double Raster<T>::maxx() const {
	return toX(cols());
}

template<class T>
double Raster<T>::miny() const {
	return toY(rows());
}

template<class T>
double Raster<T>::maxy() const {
	return toY(0);
}

template<class T>
double Raster<T>::leftx() const {
	return resolutionX() > 0 ? minx() : maxx();
}

template<class T>
double Raster<T>::rightx() const {
	return resolutionX() < 0 ? minx() : maxx();
}

template<class T>
double Raster<T>::topy() const {
	return resolutionY() > 0 ? miny() : maxy();
}

template<class T>
double Raster<T>::bottomy() const {
	return resolutionY() < 0 ? miny() : maxy();
}

template<class T>
double Raster<T>::width() const {
	return maxx() - minx();
}

template<class T>
double Raster<T>::height() const {
	return maxy() - miny();
}

template<class T>
T Raster<T>::nodata(uint16_t band) const {
	return (T) m_ds->GetRasterBand(band)->GetNoDataValue();
}

template<class T>
T Raster<T>::nodata() const {
	return nodata(1);
}

template<class T>
void Raster<T>::setNodata(const T nodata, uint16_t band) {
	m_ds->GetRasterBand(band)->SetNoDataValue((double) nodata);
}

template<class T>
void Raster<T>::setNodata(const T nodata) {
	setNodata(nodata, 1);
}

template<class T>
uint16_t Raster<T>::cols() const {
	return m_cols;
}

template<class T>
uint16_t Raster<T>::rows() const {
	return m_rows;
}

template<class T>
int32_t Raster<T>::toRow(double y) const {
	return (int32_t) ((y - m_trans[3]) / m_trans[5]);
}

template<class T>
int32_t Raster<T>::toCol(double x) const {
	return (int32_t) ((x - m_trans[0]) / m_trans[1]);
}

template<class T>
double Raster<T>::toX(int32_t col) const {
	return m_trans[0] + col * m_trans[1];
}

template<class T>
double Raster<T>::toY(int32_t row) const {
	return m_trans[3] + row * m_trans[5];
}

template<class T>
double Raster<T>::toCentroidX(int32_t col) const {
	return m_trans[0] + col * m_trans[1] + m_trans[1] / 2.0;
}

template<class T>
double Raster<T>::toCentroidY(int32_t row) const {
	return m_trans[3] + row * m_trans[5] + m_trans[5] / 2.0;
}

template<class T>
uint32_t Raster<T>::size() const {
	return (uint32_t) (m_cols * m_rows);
}

template<class T>
bool Raster<T>::isNoData(int32_t col, int32_t row, uint16_t band) {
	return get(col, row, band) == nodata(band);
}

template<class T>
bool Raster<T>::isNoData(int32_t col, int32_t row) {
	return isNoData(col, row, 1);
}

template<class T>
bool Raster<T>::isNoData(uint32_t idx, uint16_t band) {
	return get(idx, band) == nodata(band);
}

template<class T>
bool Raster<T>::isNoData(uint32_t idx) {
	return isNoData(idx, (uint16_t) 1);
}

template<class T>
bool Raster<T>::isNoData(double x, double y, uint16_t band) {
	return isNoData(toCol(x), toRow(y), band);
}

template<class T>
bool Raster<T>::isNoData(double x, double y) {
	return isNoData(x, y, 1);
}

template<class T>
bool Raster<T>::isValid(int32_t c, int32_t r, uint16_t band) {
	return getOrNodata(c, r, band) != nodata(band);
}

template<class T>
bool Raster<T>::isValid(double x, double y, uint16_t band) {
	return getOrNodata(x, y, band) != nodata(band);
}

template<class T>
T Raster<T>::getOrNodata(double x, double y, uint16_t band) {
	if (!has(x, y)) {
		return nodata(band);
	} else {
		return get(x, y, band);
	}
}

template<class T>
T Raster<T>::getOrNodata(int32_t col, int32_t row, uint16_t band) {
	if (!has(col, row)) {
		return nodata(band);
	} else {
		return get(col, row, band);
	}
}

template<class T>
T *Raster<T>::grid() {
	g_implerr("grid() Not implemented in Raster.");
}

template<class T>
bool Raster<T>::hasGrid() const {
	return false;
}

template<class T>
T Raster<T>::get(double x, double y, uint16_t band) {
	return get(toCol(x), toRow(y), band);
}

template<class T>
T Raster<T>::get(double x, double y) {
	return get(x, y, 1);
}

template<class T>
T Raster<T>::get(int32_t col, int32_t row, uint16_t band) {
	return m_cache.get(band, col, row);
}

template<class T>
T Raster<T>::get(int32_t col, int32_t row) {
	return get(col, row, 1);
}

template<class T>
T Raster<T>::get(uint32_t idx, uint16_t band) {
	if (idx >= size())
		g_argerr("Index out of bounds.");
	return get(idx % m_cols, (int) idx / m_cols, band);
}

template<class T>
T Raster<T>::get(uint32_t idx) {
	return get(idx, (uint16_t) 1);
}

template<class T>
void Raster<T>::set(int32_t col, int32_t row, T v, uint16_t band) {
	if (!m_writable)
		g_runerr("This raster is not writable.");
	m_cache.set(band, col, row, v);
}

template<class T>
void Raster<T>::set(int32_t col, int32_t row, T v) {
	set(col, row, v, 1);
}

template<class T>
void Raster<T>::set(uint32_t idx, T v, uint16_t band) {
	if (idx >= size())
		g_argerr("Index out of bounds.");
	set(idx % m_cols, (int) idx / m_cols, v, band);
}

template<class T>
void Raster<T>::set(uint32_t idx, T v) {
	set(idx, v, (uint16_t) 1);
}

template<class T>
void Raster<T>::set(double x, double y, T v, uint16_t band) {
	set(toCol(x), toRow(y), v, band);
}

template<class T>
void Raster<T>::set(double x, double y, T v) {
	set(x, y, v, 1);
}

template<class T>
bool Raster<T>::isSquare() const {
	return cols() == rows();
}

template<class T>
bool Raster<T>::has(int32_t col, int32_t row) const {
	return col >= 0 && col < m_cols && row >= 0 && row < m_rows;
}

template<class T>
bool Raster<T>::has(double x, double y) const {
	return has(toCol(x), toRow(y));
}

template<class T>
bool Raster<T>::has(uint32_t idx) const {
	return idx < (uint32_t) (m_cols * m_rows);
}

// Keeps track of polygonization progress and provides a cancellation mechanism.
class __PolyProgressData {
public:
	Callbacks *cb;
	bool *cancel;
	bool _cancel;
	__PolyProgressData(Callbacks *cb, bool *cancel) :
		cb(cb), cancel(cancel), _cancel(false) {
		if (cancel == nullptr)
			cancel = &_cancel;
	}
};

// Callback for polygonization.
int __polyProgress(double dfComplete, const char *pszMessage, void *pProgressArg) {
	__PolyProgressData *data = (__PolyProgressData *)pProgressArg;
	if (data) {
		data->cb->stepCallback((float)dfComplete);
		// TODO: The message doesn't really do anything.
		// if(pszMessage != NULL)
		//	data->cb->statusCallback(std::string(pszMessage));
		return *(data->cancel) ? 0 : 1;
	}
	return 1;
}

template<class T>
void Raster<T>::polygonize(const std::string &filename, uint16_t band, Callbacks *callbacks, bool *cancel) {
	Util::rm(filename);
	GDALAllRegister();
	GDALDriver *drv = GetGDALDriverManager()->GetDriverByName("SQLite");
	GDALDataset *ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	ds->SetProjection(m_ds->GetProjectionRef());
	double trans[6];
	m_ds->GetGeoTransform(trans);
	ds->SetGeoTransform(trans);
	char **opts = nullptr;
	opts = CSLSetNameValue(opts, "FORMAT", "SPATIALITE");
	opts = CSLSetNameValue(opts, "GEOMETRY_NAME", "geom");
	opts = CSLSetNameValue(opts, "SPATIAL_INDEX", "YES");
	OGRLayer *layer = ds->CreateLayer("boundary", NULL, wkbMultiPolygon, opts);
	OGRFieldDefn field( "id", OFTInteger);
	layer->CreateField(&field);
	if (callbacks) {
		__PolyProgressData pd(callbacks, cancel);
		GDALPolygonize(m_ds->GetRasterBand(1), NULL, layer, 0, NULL, &__polyProgress, &pd);
	} else {
		GDALPolygonize(m_ds->GetRasterBand(1), NULL, layer, 0, NULL, NULL, NULL);
	}
	GDALClose(ds);
}

template<class T>
void Raster<T>::flush() {
	if (m_writable)
		m_cache.flush();
}

template<class T>
Raster<T>::~Raster() {
	m_cache.close();
	if(m_ds)
		GDALClose(m_ds);
}

template class geotools::raster::BlockCache<float>;
template class geotools::raster::BlockCache<double>;
template class geotools::raster::BlockCache<uint64_t>;
template class geotools::raster::BlockCache<uint32_t>;
template class geotools::raster::BlockCache<uint16_t>;
template class geotools::raster::BlockCache<uint8_t>;
template class geotools::raster::BlockCache<int64_t>;
template class geotools::raster::BlockCache<int32_t>;
template class geotools::raster::BlockCache<int16_t>;
template class geotools::raster::BlockCache<int8_t>;
template class geotools::raster::BlockCache<char>;

template class geotools::raster::Grid<float>;
template class geotools::raster::Grid<double>;
template class geotools::raster::Grid<uint64_t>;
template class geotools::raster::Grid<uint32_t>;
template class geotools::raster::Grid<uint16_t>;
template class geotools::raster::Grid<uint8_t>;
template class geotools::raster::Grid<int64_t>;
template class geotools::raster::Grid<int32_t>;
template class geotools::raster::Grid<int16_t>;
template class geotools::raster::Grid<int8_t>;
template class geotools::raster::Grid<char>;

template class geotools::raster::Raster<float>;
template class geotools::raster::Raster<double>;
template class geotools::raster::Raster<uint64_t>;
template class geotools::raster::Raster<uint32_t>;
template class geotools::raster::Raster<uint16_t>;
template class geotools::raster::Raster<uint8_t>;
template class geotools::raster::Raster<int64_t>;
template class geotools::raster::Raster<int32_t>;
template class geotools::raster::Raster<int16_t>;
template class geotools::raster::Raster<int8_t>;
template class geotools::raster::Raster<char>;

template class geotools::raster::MemRaster<float>;
template class geotools::raster::MemRaster<double>;
template class geotools::raster::MemRaster<uint64_t>;
template class geotools::raster::MemRaster<uint32_t>;
template class geotools::raster::MemRaster<uint16_t>;
template class geotools::raster::MemRaster<uint8_t>;
template class geotools::raster::MemRaster<int64_t>;
template class geotools::raster::MemRaster<int32_t>;
template class geotools::raster::MemRaster<int16_t>;
template class geotools::raster::MemRaster<int8_t>;
template class geotools::raster::MemRaster<char>;

template class geotools::raster::TargetOperator<float>;

template class geotools::raster::FillOperator<float>;
