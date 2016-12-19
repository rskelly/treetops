#include <unordered_set>
#include <unordered_map>
#include <list>
#include <memory>
#include <iostream>
#include <iomanip>

#include "sqlite.hpp"
#include "geotools.hpp"
#include "raster.hpp"
#include "spectral.hpp"

using namespace geotools::spectral;
using namespace geotools::spectral::config;
using namespace geotools::db;
using namespace geotools::raster;

/**
 * Represents a single pixel across all bands.
 * Contains the x/y position of the pixel center, and a list
 * of digital numbers corresponding to the bands, in order.
 */
class Px {
public:
	double x;
	double y;
	unsigned int id;
	std::map<int, unsigned short> dn;

	/**
	 * Construct the Px with the position.
	 */
	Px(unsigned int id, double x, double y) :
			x(x), y(y), id(id) {
	}

	/**
	 * Print to the output stream. Prints the x, y and
	 * all band in order, comma-delimited.
	 */
	void print(std::ostream &out) {
		out << id << "," << x << "," << y;
		for (const auto &it : dn)
			out << "," << it.second;
		out << std::endl;
	}

	void save(SQLite &db) {
		std::map<std::string, std::string> fields;
		fields["id"] = std::to_string(id);
		for (const auto &it : dn)
			fields[std::to_string(it.first)] = std::to_string(it.second);
		db.addPoint(x, y, 0, fields);
	}

	/**
	 * Returns true if the pixel should be abandoned because it contains nodata.
	 */
	bool isNoData(unsigned int nodata) {
		// TODO: Consult on the best way to determine whether a pixel is
		// abandoned because of nodata.
		return dn[0] == nodata;
	}
};

/**
 * Computes a band list. If the given list is empty, populates the output list with
 * all of the avaible bands. Otherwise returns a list of the desired bands.
 */
std::vector<int> computeBandList(const std::set<int> &bands,
		const std::string &specFilename) {
	std::vector<int> b;
	b.assign(bands.begin(), bands.end());
	if (b.empty()) {
		Raster<unsigned int> tmp(specFilename);
		for (int i = 1; i <= tmp.bandCount(); ++i)
			b.push_back(i);
		g_warn("No bands selected; using all " << b.size() << " bands.");
	}
	return b;
}

/**
 * Computes the bounding box which describes the overlapping area between the index
 * raster and the spectral raster.
 */
Bounds computeOverlapBounds(Raster<unsigned int> &idxRaster,
		const std::string &specFilename) {
	Raster<float> tmp(specFilename);
	return idxRaster.bounds().intersection(tmp.bounds());
}

/**
 * Process a single spectral file.
 * config       -- A SpectralConfig object with operational parameters.
 * idxRaster    -- A Raster containing the polygon IDs that are used to group the extracted values.
 * specFilename -- The filename of the current spectral file.
 */
void processSpectralFile(const SpectralConfig &config,
		const std::vector<int> &bands, const std::string &specFilename,
		SQLite &db) {

	g_debug("processSpectralFile [config] [raster] " << specFilename);

	uint16_t startRow, endRow, startCol, endCol;
	uint32_t idxNodata = config.idxNodata;
	uint16_t specNodata = config.specNodata;
	Bounds bounds;

	{
		Raster<uint32_t> idxRaster(config.indexFilename);
		// Get the start and end cols/rows from the bounds intersection.
		bounds.collapse();
		bounds.extend(computeOverlapBounds(idxRaster, specFilename));
		startRow = idxRaster.toRow(
				idxRaster.resolutionY() > 0 ? bounds.miny() : bounds.maxy());
		endRow = idxRaster.toRow(
				idxRaster.resolutionY() > 0 ? bounds.maxy() : bounds.miny());
		startCol = idxRaster.toCol(bounds.minx());
		endCol = idxRaster.toCol(bounds.maxx());
		// If nodata is not set, get it from the raster.
		if (!config.hasIdxNodata)
			idxNodata = idxRaster.nodata();
	}

	// If nodata is not set, get it from the raster.
	if (!config.hasSpecNodata) {
		Raster<uint16_t> specRaster(specFilename);
		specNodata = specRaster.nodata();
	}

	//std::unordered_map<unsigned int, std::unique_ptr<Poly> > polys;
	std::unordered_map<uint64_t, std::unique_ptr<Px> > px;
	std::unordered_set<uint64_t> skip;
	std::unordered_set<uint64_t> finalize;

	// Iterate over bands to populate Polys
	#pragma omp parallel for
	for (int block = startRow; block < endRow; block += 100) {

		Raster<unsigned int> idxRaster(config.indexFilename);
		Raster<unsigned short> specRaster(specFilename);

		for (int32_t row = block; row < g_min(endRow, block + 100); ++row) {
			for (uint32_t b = 0; b < bands.size(); ++b) {
				uint16_t band = bands[b];
				for (uint32_t col = startCol; col < endCol; ++col) {
					uint32_t id = idxRaster.get(col, row, 1);
					if (id == idxNodata)
						continue;
					double x = idxRaster.toX(col)
							+ idxRaster.resolutionX() / 2.0;
					double y = idxRaster.toY(row)
							+ idxRaster.resolutionY() / 2.0;
					int scol = specRaster.toCol(x);
					int srow = specRaster.toRow(y);
					if (specRaster.has(scol, srow)) {
						unsigned int v = specRaster.get(scol, srow);
						if (v == specNodata)
							continue;
						uint64_t idx = ((uint64_t) col << 32) | row;
						#pragma omp critical (PX)
						{
							if (skip.find(idx) == skip.end()) {
								if (px.find(idx) == px.end()) {
									std::unique_ptr<Px> pp(new Px(id, x, y));
									px[idx] = std::move(pp);
								}
								px[idx]->dn[band] = v;
								//g_debug("px id: " << id << "; idx: " << idx << "; dn size: " << p.dn.size());
								if (px[idx]->dn.size() == bands.size())
									finalize.insert(idx);
							} else {
								skip.emplace(idx);
							}
						}
					}
				}
			}
		}
		// Remove polys that have IDs that are not in the current row.	
		std::cout << std::fixed << std::setprecision(3);
		#pragma omp critical (PX)
		{
			g_debug("Save " << finalize.size());
			db.begin();
			for (const uint32_t &idx : finalize) {
				if (px[idx]->dn.size() == bands.size())
					px[idx]->save(db);
			}
			db.commit();
		}
		#pragma omp critical (PX)
		{
			g_debug("Finalize.");
			for (const uint32_t &idx : finalize)
				px.erase(idx);
			finalize.clear();
		}

	}
}

Spectral::Spectral() :
		m_callbacks(nullptr),
		m_cancel(nullptr) {
}

bool __spec_cancel = false;

void Spectral::extractSpectra(const SpectralConfig &config,
		Callbacks *callbacks, bool *cancel) {

	// Check the config for problems.
	config.check();

	m_cancel = cancel == nullptr ? &__spec_cancel : cancel;
	m_callbacks = callbacks;


	// Check the bands list; fill it if necessary.
	std::vector<int> bands = computeBandList(config.bands,
			config.spectralFilenames[0]);

	// Develop the fields list for the DB.
	std::map<std::string, int> fields;
	fields["id"] = SQLite::INTEGER;
	for (const int &band : bands)
		fields["b" + std::to_string(band)] = SQLite::INTEGER;

	// Initialize the SQLite table.
	SQLite db(config.outputFilename, SQLite::POINT, config.srid, fields);

	// Do the work.
	for (const std::string &specFilename : config.spectralFilenames) {
		g_debug(" - file: " << specFilename);
		processSpectralFile(config, bands, specFilename, db);
	}
}

