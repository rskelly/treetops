#include <queue>
#include <string>
#include <fstream>
#include <atomic>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

#include <gdal_alg.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_port.h>

#include "geotools.hpp"
#include "util.hpp"
#include "raster.hpp"

using namespace geotools::util;
using namespace geotools::raster;

bool _cancel = false;

int _getTypeSize(DataType type) {
	switch(type) {
	case DataType::Byte: return sizeof(uint8_t);
	case DataType::Float32: return sizeof(float);
	case DataType::Float64: return sizeof(double);
	case DataType::Int16: return sizeof(int16_t);
	case DataType::Int32: return sizeof(int32_t);
	case DataType::UInt16: return sizeof(uint16_t);
	case DataType::UInt32: return sizeof(uint32_t);
	default:
		g_runerr("No size for type: " << type);
	}
}

GDALDataType _dataType2GDT(DataType type) {
	switch(type) {
	case DataType::Byte:  	return GDT_Byte;
	case DataType::UInt16: 	return GDT_UInt16;
	case DataType::UInt32:	return GDT_UInt32;
	case DataType::Int16:	return GDT_Int16;
	case DataType::Int32:	return GDT_Int32;
	case DataType::Float64:	return GDT_Float64;
	case DataType::Float32:	return GDT_Float32;
	case DataType::None:
	default:
		break;
	}
	return GDT_Unknown;
}

DataType _gdt2DataType(GDALDataType type) {
	switch(type) {
	case GDT_Byte:	  	return DataType::Byte;
	case GDT_UInt16: 	return DataType::UInt16;
	case GDT_UInt32:	return DataType::UInt32;
	case GDT_Int16:		return DataType::Int16;
	case GDT_Int32:		return DataType::Int32;
	case GDT_Float64:	return DataType::Float64;
	case GDT_Float32:	return DataType::Float32;
	case GDT_Unknown:
	case GDT_CInt16:
	case GDT_CInt32:
	case GDT_CFloat32:
	case GDT_CFloat64:
	case GDT_TypeCount:
	default:
		break;
	}
	return DataType::None;
}

GridProps::GridProps() :
		m_cols(0), m_rows(0),
		m_vsrid(0), m_hsrid(0),		// Vertical and horizontal srid
		m_bands(0),           		// The number of bands
		m_writable(false),			// True if the grid is writable
		m_nodata(0),
		m_type(DataType::None) {	// The data type.
}

void GridProps::setNoData(double nodata) {
	m_nodata = nodata;
}

bool GridProps::hasCell(int col, int row) const {
	return !(col < 0 || row < 0 || row >= m_rows || col >= m_cols);
}

bool GridProps::hasCell(double x, double y) const {
	return hasCell(toCol(x), toRow(y));
}

int GridProps::toRow(double y) const {
	return (int) ((y - m_trans[3]) / m_trans[5]);
}

int GridProps::toCol(double x) const {
	return (int) ((x - m_trans[0]) / m_trans[1]);
}

double GridProps::toX(int col) const {
	return m_trans[0] + col * m_trans[1];
}


double GridProps::toY(int row) const {
	return m_trans[3] + row * m_trans[5];
}

void GridProps::setResolution(double resolutionX, double resolutionY) {
	m_trans[1] = resolutionX;
	m_trans[5] = resolutionY;
}

double GridProps::resolutionX() const {
	return m_trans[1];
}

double GridProps::resolutionY() const {
	return m_trans[5];
}

void GridProps::setDataType(DataType type) {
	m_type = type;
}

DataType GridProps::dataType() const {
	return m_type;
}

void GridProps::setSize(int cols, int rows) {
	m_cols = cols;
	m_rows = rows;
}

int GridProps::cols() const {
	return m_cols;
}

int GridProps::rows() const {
	return m_rows;
}

void GridProps::setSrid(int hsrid, int vsrid) {
	m_vsrid = vsrid;
	m_hsrid = hsrid;
}

int GridProps::vsrid() const {
	return m_vsrid;
}

int GridProps::hsrid() const {
	return m_vsrid;
}

void GridProps::setProjection(const std::string &proj) {
	m_projection = proj;
}

std::string GridProps::projection() const {
	if(m_projection.empty() && m_hsrid > 0) {
		OGRSpatialReference ref;
		ref.importFromEPSG(m_hsrid);
		char *proj;
		ref.exportToWkt(&proj);
		return std::string(proj);
	} else {
		return m_projection;
	}
}

void GridProps::setTrans(double trans[6]) {
	for(int i = 0; i < 6; ++i)
		m_trans[i] = trans[i];
}

void GridProps::trans(double trans[6]) const {
	for(int i = 0; i < 6; ++i)
		trans[i] = m_trans[i];
}

void GridProps::setBands(int bands) {
	m_bands = bands;
}

int GridProps::bands() const {
	return m_bands;
}

void GridProps::setWritable(bool writable) {
	m_writable = writable;
}

bool GridProps::writable() const {
	return m_writable;
}


// Implementations for Cell

Cell::Cell(int col, int row) :
	col(col), row(row) {
}

// Implementations for TargetOperator (for flood fill)

FillOperator::~FillOperator() {}

TargetOperator::TargetOperator(int match) :
		m_match(match) {
}

bool TargetOperator::fill(int value) const {
	return value == m_match;
}

TargetOperator::~TargetOperator() {}

Grid::Grid() {
}

// Implementations forthe Grid class

void Grid::gaussianWeights(double *weights, int size, double sigma) {
	// If size is an even number, bump it up.
	if (size % 2 == 0) {
		++size;
		g_warn("Gaussian kernel size must be an odd number >=3. Bumping up to " << size);
	}
	for (int r = 0; r < size; ++r) {
		for (int c = 0; c < size; ++c) {
			int x = size / 2 - c;
			int y = size / 2 - r;
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma))
					* pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
		}
	}
}

GridStats Grid::stats() {
	GridStats st;
	long i;
	double v, nodata = props().nodata();
	for (i = 0; i < props().size(); ++i) {
		if ((v = getFloat(i)) != nodata) {
			st.min = st.max = v;
			break;
		}
	}
	st.sum = 0;
	st.count = 0;
	double m = 0;
	double s = 0;
	int k = 1;
	// Welford's method for variance.
	// i has the index of the first non-nodata element.
	for (; i < props().size(); ++i) {
		if ((v = getFloat(i)) != nodata) {
			double oldm = m;
			m = m + (v - m) / k;
			s = s + (v - m) * (v - oldm);
			st.sum += v;
			st.min = g_min(st.min, v);
			st.max = g_max(st.max, v);
			++st.count;
			++k;
		}
	}
	st.mean = st.sum / st.count;
	st.variance = s / st.count;
	st.stdDev = std::sqrt(st.variance);
	return st;
}

void Grid::normalize(int band) {
	GridStats st = stats();
	double v, nodata = props().nodata();
	double mean = st.mean;
	double stdDev = st.stdDev;
	for (long i = 0; i < props().size(); ++i) {
		if ((v = getFloat(i, band)) != nodata && !std::isnan(v) && v < G_DBL_MAX_POS) {
			setFloat(i, ((v - mean) / stdDev), band);
		} else {
			setFloat(i, nodata, band);
		}
	}
}

void Grid::logNormalize(int band) {
	GridStats st = stats();
	double n = st.min;
	double x = st.max;
	double e = std::exp(1.0) - 1.0;
	for(long i = 0; i < props().size(); ++i)
		setFloat(i, std::log(1.0 + e * (getFloat(i) - n) / (x - n)));
}


void Grid::floodFill(int col, int row,
    FillOperator &op, Grid &other, int fill, bool d8,
	int *outminc, int *outminr,	int *outmaxc, int *outmaxr,
	int *outarea) {

	if(!props().isInt() || !other.props().isInt())
		g_argerr("Flood fill is only implemented for integer rasters.");

	int cols = props().cols();
	int rows = props().rows();
	int size = props().size();
	int minc = cols + 1;
	int minr = props().rows() + 1;
	int maxc = -1;
	int maxr = -1;
	int area = 0;
	std::queue<std::unique_ptr<Cell> > q;
	q.push(std::unique_ptr<Cell>(new Cell(col, row)));

	std::vector<bool> visited(size, false); // Tracks visited pixels.

	while (q.size()) {

		std::unique_ptr<Cell> cel = std::move(q.front());
		row = cel->row;
		col = cel->col;
		q.pop();

		uint64_t idx = (uint64_t) row * cols + col;

		if (!visited[idx] && op.fill(getInt(col, row))) {

			minc = g_min(col, minc);
			maxc = g_max(col, maxc);
			minr = g_min(row, minr);
			maxr = g_max(row, maxr);
			++area;
			other.setInt(col, row, fill);
			visited[idx] = true;

			if (row > 0)
				q.push(std::unique_ptr<Cell>(new Cell(col, row - 1)));
			if (row < rows - 1)
				q.push(std::unique_ptr<Cell>(new Cell(col, row + 1)));

			int c;
			for (c = col - 1; c >= 0; --c) {
				idx = (uint64_t) row * cols + c;
				if (!visited[idx] && op.fill(getInt(c, row))) {
					minc = g_min(c, minc);
					++area;
					other.setInt(c, row, fill);
					visited[idx] = true;
					if (row > 0)
						q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
					if (row < rows - 1)
						q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
				} else {
					break;
				}
			}
			if(d8) {
				if (row > 0)
					q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
				if (row < rows - 1)
					q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
			}
			for (c = col + 1; c < cols; ++c) {
				idx = (uint64_t) row * cols + c;
				if (!visited[idx] && op.fill(getInt(c, row))) {
					maxc = g_max(c, maxc);
					++area;
					other.setInt(c, row, fill);
					visited[idx] = true;
					if (row > 0)
						q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
					if (row < rows - 1)
						q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
				} else {
					break;
				}
			}
			if(d8) {
				if (row > 0)
					q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
				if (row < rows - 1)
					q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
			}
		}
	}
	if(outminc != nullptr)
		*outminc = minc;
	if(outminr != nullptr)
		*outminr = minr;
	if(outmaxc != nullptr)
		*outmaxc = maxc;
	if(outmaxr != nullptr)
		*outmaxr = maxr;
	if(outarea != nullptr)
		*outarea = area;
}

void Grid::floodFill(int col, int row, int target, int fill, bool d8,
		int *outminc, int *outminr,
		int *outmaxc, int *outmaxr,
		int *outarea) {
	TargetOperator op(target);
	return floodFill(col, row, op, *this, fill, d8, outminc, outminr,
			outmaxc, outmaxr, outarea);
}

void Grid::floodFill(int col, int row,
		FillOperator &op, int fill, bool d8,
		int *outminc, int *outminr,
		int *outmaxc, int *outmaxr,
		int *outarea) {
	return floodFill(col, row, op, *this, fill, d8, outminc, outminr,
			outmaxc, outmaxr, outarea);
}

void Grid::convert(Grid &g, int srcBand, int dstBand) {
	if(g.props().isInt()) {
		for (long i = 0; i < props().size(); ++i)
			g.setInt(i, getInt(i, srcBand), dstBand);
	} else {
		for (long i = 0; i < props().size(); ++i)
			g.setFloat(i, getFloat(i, srcBand), dstBand);
	}
}

void Grid::voidFillIDW(double radius, int count, double exp, int band) {

	if(!props().isFloat())
		g_runerr("IDW fill only implemented for float rasters.");

	if (radius <= 0.0)
		throw std::invalid_argument("Radius must be larger than 0.");

	if (count <= 0)
		throw std::invalid_argument("Count must be larger than 0.");

	if (exp <= 0.0)
		throw std::invalid_argument("Exponent must be larger than 0.");

	MemRaster tmp(props());
	double nodata = props().nodata();
	int rows = props().rows();
	int cols = props().cols();
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < cols; ++c) {
			if (getFloat(c, r) != nodata)
				continue;
			double rad = radius;
			bool found = false;
			do {
				double d = g_sq(rad);
				double a = 0.0;
				double b = 0.0;
				int cnt = 0;
				for (int r0 = (int) g_max(0, r - rad); r0 < (int) g_min(rows, r + rad + 1); ++r0) {
					for (int c0 = (int) g_max(0, c - rad); c0 < (int) g_min(cols, c + rad + 1); ++c0) {
						double d0 = g_sq((double) c0 - c) + g_sq((double) r0 - r);
						if (d0 <= d && getFloat(c0, r0, band) != nodata) {
							double dp = 1.0 / std::pow(d0, exp);
							a += dp * getFloat(c0, r0, band);
							b += dp;
							++cnt;
						}
					}
				}

				if (cnt >= count) {
					tmp.setFloat(c, r, (a / b), band);
					found = true;
					break;
				}

				rad += 1.0;

			} while (rad < g_min(cols, rows));

			if (!found)
				g_warn("Pixel not filled at " << c << "," << r << ". Consider larger radius or smaller count.");
		}
	}
	writeBlock(tmp);
}

void Grid::smooth(Grid &smoothed, double sigma, int size, int band,
		Callbacks *status, bool *cancel) {

	if(!props().isFloat())
		g_runerr("Smoothing only implemented for float rasters.");

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

	double nd = props().nodata();

	Buffer weights(size * size * _getTypeSize(DataType::Float64));
	Grid::gaussianWeights((double *) weights.buf, size, sigma);

	std::atomic<int> curRow(0);

	if (status)
		status->stepCallback(0.02f);

	int cols = props().cols();
	int rows = props().rows();

	#pragma omp parallel
	{

		int bufSize = g_max(size, 1024);

		GridProps pr;
		pr.setSize(cols, bufSize + size);
		pr.setNoData(nd);
		pr.setDataType(DataType::Float64);
		MemRaster buf(pr, false);
		buf.fillFloat(nd);

		MemRaster smooth(pr, false);
		smooth.fillFloat(nd);

		#pragma omp for
		for (int i = 0; i < (rows - size) / bufSize + 1; ++i) {
			if (*cancel) continue;

			int b = i * bufSize;

			if(status)
				status->statusCallback("Reading...");

			buf.fillFloat(nd);
			smooth.fillFloat(nd);

			int readOffset = b > 0 ? b - size / 2 : 0;  // If this is the first row, read from zero, otherwise -(size / 2)
			int writeOffset = b > 0 ? 0 : size / 2;     // If this is the first row, write to (size / 2), otherwise 0.
			#pragma omp critical(__smooth_read)
			readBlock(buf, pr.cols(), pr.rows(), 1, 0, readOffset, 0, writeOffset);

			if(status)
				status->statusCallback("Processing...");

			// Process the entire block, even the buffer parts.
			for (int r = 0; r < pr.rows() - size; ++r) {
				for (int c = 0; c < pr.cols() - size; ++c) {
					double v, t = 0.0;
					bool foundNodata = false;
					for (int gr = 0; gr < size; ++gr) {
						for (int gc = 0; gc < size; ++gc) {
							v = buf.getFloat(c + gc, r + gr);
							if (v == nd) {
								foundNodata = true;
								break;
							} else {
								t += *(((double *) weights.buf) + gr * size + gc) * v;
							}
						}
						if(foundNodata) break;
					}
					if (!foundNodata)
						smooth.setFloat(c + size / 2, r + size / 2, t);
				}
			}

			if (status) {
				curRow += bufSize;
				status->stepCallback(0.2f + (float) curRow / rows * 0.97f);
				status->statusCallback("Writing...");
			}

			#pragma omp critical(__smooth_write)
			smoothed.writeBlock(smooth, cols, g_min(bufSize, rows - b - size), 0, b, 0, size / 2); // Always write to b and read from (size / 2)
		}
	}
	if (*cancel)
		return;

	if (status) {
		status->stepCallback(1.0);
		status->statusCallback("Closing...");
	}
}

Grid::~Grid() {
}

// Implementations for MemRaster

void MemRaster::checkInit() const {
	if (m_grid == nullptr)
		g_runerr("This instance has not been initialized.");
}

MemRaster::MemRaster() :
	m_grid(nullptr),
	m_mmapped(false) {
}

MemRaster::MemRaster(const GridProps &props, bool mapped) :
	m_grid(nullptr),
	m_mmapped(mapped) {
	init(props);
}

MemRaster::~MemRaster() {
	freeMem();
}

void* MemRaster::grid() {
	return m_grid;
}

void MemRaster::freeMem() {
	if (m_grid) {
		if (m_mmapped) {
			m_mappedFile.release();
		} else {
			free(m_grid);
		}
	}
}

void MemRaster::init(const GridProps &pr, bool mapped) {
	m_grid = nullptr;
	m_mmapped = false;
	if (pr.cols() != props().cols() || pr.rows() != props().rows()) {
		freeMem();
		m_props = pr;
		m_mmapped = mapped;
		m_grid = nullptr;
		int size = _getTypeSize(pr.dataType()) * pr.cols() * pr.rows();
		m_mappedFile.release();
		if (mapped) {
			const std::string filename = Util::tmpFile("/tmp");
			m_mappedFile = Util::mapFile(filename, size);
			m_grid = m_mappedFile->data();
		} else {
			m_grid = malloc(size);
		}
		if (!m_grid)
			g_runerr("Failed to allocate memory for MemRaster.");
	}
}

void MemRaster::fillFloat(double value, int band) {
	checkInit();
	for(long i = 0; i < props().size(); ++i)
		*(((double *) m_grid) + i) = value;
}

void MemRaster::fillInt(int value, int band) {
	checkInit();
	for(long i = 0; i < props().size(); ++i)
		*(((int *) m_grid) + i) = value;
}

double MemRaster::getFloat(long idx, int band) {
	checkInit();
	if (idx < 0 || idx >= props().size())
		g_argerr("Index out of bounds: " << idx << "; size: " << props().size());
	return *(((double *) m_grid) + idx);
}

double MemRaster::getFloat(int col, int row, int band) {
	long idx = (long) row * props().cols() + col;
	return getFloat(idx, band);
}

void MemRaster::setFloat(int col, int row, double value, int band) {
	long idx = (long) row * props().cols() + col;
	setFloat(idx, value, band);
}

void MemRaster::setFloat(long idx, double value, int band) {
	checkInit();
	if (idx >= props().size())
		g_argerr("Index out of bounds: " << idx << "; size: " << props().size()
						<< "; value: " << value << "; col: " << (idx % props().cols())
						<< "; row: " << (idx / props().cols()));
	*(((double *) m_grid) + idx) = value;
}

void MemRaster::toMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band) {
	int cols = props().cols();
	int rows = props().rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			mtx(r, c) = getFloat(c, r, band);
	}
}

void MemRaster::fromMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band) {
	int cols = props().cols();
	int rows = props().rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			setFloat(c, r, (double) mtx(r, c), band);
	}
}

void MemRaster::readBlock(Grid &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {

	if (&grd == this)
		g_argerr("Recursive call to readBlock.");

	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	g_runerr("Not implemented.");
}


void MemRaster::writeBlock(Grid &grd,
            		int cols, int rows,
            		int srcCol, int srcRow,
					int dstCol, int dstRow,
					int srcBand, int dstBand) {

	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	g_runerr("Not implemented.");
}



// Implementations for Raster

Raster::Raster(const std::string &filename, const GridProps &props) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_block(nullptr),
		m_type(GDT_Unknown) {

	if (props.resolutionX() == 0 || props.resolutionY() == 0)
		g_argerr("Resolution must not be zero.");
	if (props.cols() <= 0 || props.rows() <= 0)
		g_argerr("Columns and rows must be larger than zero.");
	if (filename.empty())
		g_argerr("Filename must be given.");

	m_props = props;
	m_filename = filename;

	// Create GDAL dataset.
	char **opts = NULL;
	// TODO: Compress option in props
	// opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
	opts = CSLSetNameValue(opts, "BIGTIFF", "YES");
	GDALAllRegister();
	m_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
			filename.c_str(), m_props.cols(), m_props.rows(), m_props.bands(),
			_dataType2GDT(m_props.dataType()), opts);
	if (!m_ds)
		g_runerr("Failed to create file.");

	// Initialize geotransform.
	double trans[6];
	m_props.trans(trans);
	m_ds->SetGeoTransform(trans);

	// Set projection.
	std::string proj = m_props.projection();
	if (!proj.empty())
		m_ds->SetProjection(proj.c_str());

	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
 	m_block = malloc(m_bcols * m_brows * _getTypeSize(m_props.dataType()));
}

Raster::Raster(const std::string &filename, bool writable) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_block(nullptr),
		m_type(GDT_Unknown) {

	if (filename.empty())
		g_argerr("Filename must be given.");

	m_filename = filename;

	// Attempt to open the dataset.
	GDALAllRegister();
	m_ds = (GDALDataset *) GDALOpen(filename.c_str(),
			writable ? GA_Update : GA_ReadOnly);
	if (m_ds == NULL)
		g_runerr("Failed to open raster.");

	// Save some raster properties
	double trans[6];
	m_ds->GetGeoTransform(trans);
	m_props.setTrans(trans);
	m_props.setSize(m_ds->GetRasterXSize(), m_ds->GetRasterYSize());
	m_props.setDataType(_gdt2DataType(m_ds->GetRasterBand(1)->GetRasterDataType()));
	m_props.setBands(m_ds->GetRasterCount());
	m_props.setWritable(writable);
	m_props.setProjection(std::string(m_ds->GetProjectionRef()));

	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
 	m_block = malloc(m_bcols * m_brows * _getTypeSize(m_props.dataType()));
}

const GridProps& Raster::props() const {
	return m_props;
}

DataType Raster::getFileDataType(const std::string &filename) {
	GDALDataset *ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
	DataType type = _gdt2DataType(ds->GetRasterBand(1)->GetRasterDataType());
	GDALClose(ds);
	return type;
}

std::string Raster::filename() const {
	return m_filename;
}

void Raster::fill(double value, int band) {
	Buffer buf(m_bcols * m_brows * _getTypeSize(props().dataType()));
	for(int i = 0; i < m_bcols * m_brows; ++i)
		*(((double *) buf.buf) + i) = value;
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	for (int r = 0; r < m_brows; ++r) {
		for (int c = 0; c < m_bcols; ++c) {
			if (bnd->WriteBlock(c, r, buf.buf) != CE_None)
				g_runerr("Fill error.");
		}
	}
}

void Raster::readBlock(Raster &grd, int cols, int rows,
		int srcBand, int dstBand,
		int srcCol, int srcRow,
		int dstCol, int dstRow) {
	if (&grd == this)
		g_runerr("Recursive call to readBlock.");
	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);
	if(dstBand < 1 || dstBand > grd.props().bands())
		g_argerr("Invalid destination band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	// A buffer for the entire copy.
	int typeSize = _getTypeSize(props().dataType());
	Buffer buf(cols * rows * typeSize);

	if(CPLE_None != m_ds->GetRasterBand(srcBand)->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
			buf.buf, cols, rows, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to read from: " << filename());
	if(CPLE_None != grd.m_ds->GetRasterBand(dstBand)->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
			buf.buf, cols, rows, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to write to: " << grd.filename());
}


void Raster::readBlock(MemRaster &grd, int cols, int rows,
		int srcBand, int srcCol, int srcRow, int dstCol, int dstRow) {

	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	if(CPLE_None != m_ds->GetRasterBand(srcBand)->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
			grd.grid(), dstCol, dstRow, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to read from: " << filename());
}

void Raster::writeBlock(Raster &grd, int cols, int rows,
		int srcBand, int dstBand,
		int srcCol, int srcRow,
		int dstCol, int dstRow) {
	if (&grd == this)
		g_runerr("Recursive call to readBlock.");
	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);
	if(dstBand < 1 || dstBand > grd.props().bands())
		g_argerr("Invalid destination band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	// A buffer for the entire copy.
	int typeSize = _getTypeSize(props().dataType());
	Buffer buf(cols * rows * typeSize);

	if(CPLE_None != grd.m_ds->GetRasterBand(dstBand)->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
			buf.buf, cols, rows, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to read from: " << grd.filename());
	if(CPLE_None != m_ds->GetRasterBand(srcBand)->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
			buf.buf, cols, rows, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to write to: " << filename());
}

void Raster::writeBlock(MemRaster &grd, int cols, int rows,
		int srcBand, int srcCol, int srcRow, int dstCol, int dstRow) {

	if(props().dataType() != grd.props().dataType())
		g_runerr("Rasters must have the same data type. Use convert if different.");
	if(srcBand < 1 || srcBand > props().bands())
		g_argerr("Invalid source band: " << srcBand);

	cols = g_min(g_min(cols, props().cols() - srcCol), grd.props().cols() - dstCol);
	rows = g_min(g_min(rows, props().rows() - srcRow), grd.props().rows() - dstRow);

	if(CPLE_None != m_ds->GetRasterBand(srcBand)->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
			grd.grid(), srcCol, srcRow, _dataType2GDT(props().dataType()), 0, 0, 0))
		g_runerr("Failed to write to: " << filename());
}

void _readFromBlock(void *block, GDALDataType type, double *value, int idx) {
	switch(type) {
	case GDT_Float32:
		*value = (double) *(((float *) block) + idx);
		break;
	case GDT_Float64:
		*value = (double) *(((double *) block) + idx);
		break;
	case GDT_UInt32:
		*value = (double) *(((uint32_t *) block) + idx);
		break;
	case GDT_UInt16:
		*value= (double) *(((uint16_t *) block) + idx);
		break;
	case GDT_Int32:
		*value = (double) *(((int32_t *) block) + idx);
		break;
	case GDT_Int16:
		*value = (double) *(((int16_t *) block) + idx);
		break;
	case GDT_Byte:
		*value = (double) *(((uint8_t *) block) + idx);
		break;
	default:
		g_runerr("Data type not implemented: " << type);
		break;
	}
}

void _readFromBlock(void *block, GDALDataType type, int *value, int idx) {
	switch(type) {
	case GDT_Float32:
		*value = (int) *(((float *) block) + idx);
		break;
	case GDT_Float64:
		*value = (int) *(((double *) block) + idx);
		break;
	case GDT_UInt32:
		*value = (int) *(((uint32_t *) block) + idx);
		break;
	case GDT_UInt16:
		*value= (int) *(((uint16_t *) block) + idx);
		break;
	case GDT_Int32:
		*value = (int) *(((int32_t *) block) + idx);
		break;
	case GDT_Int16:
		*value = (int) *(((int16_t *) block) + idx);
		break;
	case GDT_Byte:
		*value = (int) *(((uint8_t *) block) + idx);
		break;
	default:
		g_runerr("Data type not implemented: " << type);
		break;
	}
}

double Raster::getFloat(int col, int row, int band) {
	if (!props().writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	GDALRasterBand *rb = m_ds->GetRasterBand(band);
	if(!rb)
		g_argerr("No such band: " << band);
	if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
		g_runerr("Failed to read from: " << filename());
	int idx = (row - brow * m_brows) * m_bcols + (col - bcol * m_bcols);
	double v;
	_readFromBlock(m_block, getGDType(), &v, idx);
	return v;
}

double Raster::getFloat(long idx, int band) {
	return getFloat((int) idx % props().cols(), (int) idx / props().rows(), band);
}


double Raster::getFloat(double x, double y, int band) {
	return getFloat(props().toCol(x), props().toRow(y), band);
}

int Raster::getInt(int col, int row, int band) {
	if (!props().writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	GDALRasterBand *rb = m_ds->GetRasterBand(band);
	if(!rb)
		g_argerr("No such band: " << band);
	if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
		g_runerr("Failed to read from: " << filename());
	int idx = (row - brow * m_brows) * m_bcols + (col - bcol * m_bcols);
	int v;
	_readFromBlock(m_block, getGDType(), &v, idx);
	return v;
}

int Raster::getInt(long idx, int band) {
	return getInt((int) idx % props().cols(), (int) idx / props().rows(), band);
}


int Raster::getInt(double x, double y, int band) {
	return getInt(props().toCol(x), props().toRow(y), band);
}

void _writeToBlock(void *block, GDALDataType type, double value, int idx) {
	switch(type) {
	case GDT_Float32:
		*(((float *) block) + idx) = (float) value;
		break;
	case GDT_Float64:
		*(((double *) block) + idx) = (double) value;
		break;
	case GDT_UInt32:
		*(((uint32_t *) block) + idx) = (uint32_t) value;
		break;
	case GDT_UInt16:
		*(((uint16_t *) block) + idx) = (uint16_t) value;
		break;
	case GDT_Int32:
		*(((int32_t *) block) + idx) = (int32_t) value;
		break;
	case GDT_Int16:
		*(((int16_t *) block) + idx) = (int16_t) value;
		break;
	case GDT_Byte:
		*(((uint8_t *) block) + idx) = (uint8_t) value;
		break;
	default:
		g_runerr("Data type not implemented: " << type);
		break;
	}
}

void _writeToBlock(void *block, GDALDataType type, int value, int idx) {
	switch(type) {
	case GDT_Float32:
		*(((float *) block) + idx) = (float) value;
		break;
	case GDT_Float64:
		*(((double *) block) + idx) = (double) value;
		break;
	case GDT_UInt32:
		*(((uint32_t *) block) + idx) = (uint32_t) value;
		break;
	case GDT_UInt16:
		*(((uint16_t *) block) + idx) = (uint16_t) value;
		break;
	case GDT_Int32:
		*(((int32_t *) block) + idx) = (int32_t) value;
		break;
	case GDT_Int16:
		*(((int16_t *) block) + idx) = (int16_t) value;
		break;
	case GDT_Byte:
		*(((uint8_t *) block) + idx) = (uint8_t) value;
		break;
	default:
		g_runerr("Data type not implemented: " << type);
		break;
	}
}

void Raster::setInt(int col, int row, int v, int band) {
	if (!props().writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	GDALRasterBand *rb = m_ds->GetRasterBand(band);
	if(!rb)
		g_argerr("No such band: " << band);
	if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
		g_runerr("Failed to read from: " << filename());
	int idx = (row - brow * m_brows) * m_bcols + (col - bcol * m_bcols);
	_writeToBlock(m_block, getGDType(), v, idx);
	if(CPLE_None != rb->WriteBlock(bcol, brow, m_block))
		g_runerr("Failed to write to: " << filename());
}

void Raster::setInt(long idx, int v, int band) {
	setInt((int) idx % props().cols(), (int) idx / props().rows(), v, band);
}


void Raster::setInt(double x, double y, int v, int band) {
	setInt(props().toCol(x), props().toRow(y), v, band);
}

void Raster::setFloat(int col, int row, double v, int band) {
	if (!props().writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	GDALRasterBand *rb = m_ds->GetRasterBand(band);
	if(!rb)
		g_argerr("No such band: " << band);
	if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
		g_runerr("Failed to read from: " << filename());
	int idx = (row - brow * m_brows) * m_bcols + (col - bcol * m_bcols);
	_writeToBlock(m_block, getGDType(), v, idx);
	if(CPLE_None != rb->WriteBlock(bcol, brow, m_block))
		g_runerr("Failed to write to: " << filename());
}

void Raster::setFloat(long idx, double v, int band) {
	setFloat((int) idx % props().cols(), (int) idx / props().rows(), v, band);
}


void Raster::setFloat(double x, double y, double v, int band) {
	setFloat(props().toCol(x), props().toRow(y), v, band);
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
		return *(data->cancel) ? 0 : 1;
	}
	return 1;
}


void Raster::polygonize(const std::string &filename, int srid, int band,
		Callbacks *callbacks, bool *cancel) {
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
	OGRSpatialReference sr;
	sr.importFromEPSG(srid);
	OGRLayer *layer = ds->CreateLayer("boundary", &sr, wkbMultiPolygon, opts);
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


void Raster::flush() {
	m_ds->FlushCache();
}


Raster::~Raster() {
	GDALClose(m_ds);
}
