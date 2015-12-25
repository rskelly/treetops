/**
 * "Feathers" the edges of a data region (i.e. not nodata) to the specified
 * distance in map units, using the specified curve.
 */

//#include <float.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#define PI 3.14159265358979323846

float _min(float a, float b) {
	return a > b ? b : a;
}

float _max(float a, float b) {
	return a < b ? b : a;
}

/**
 * Returns a value between 0 and 1 following the tan curve.
 * The range of step is expected to be 0 -> steps and is clamped.
 */
float _tanCurve(float step, float steps) {
	if(step <= 0.0) return 0.0;
	if(step >= steps) return 1.0;
	return tanh(((step - steps / 2.0) / (steps / 2.0)) * PI) * 0.5 + 0.5;
}

/**
 * Returns true if the given pixel is next to a nodata pixel, or the edge of the grid.
 */
bool isEdgePixel(char *fillGrid, int col, int row, int cols, int rows) {
	if(fillGrid[row * cols + col] == 0)
		return false;
	for(int r = row - 1; r < row + 2; ++r) {
		for(int c = col - 1; c < col + 2; ++c) {
			if(c <= 0 || r <= 0 || c >= cols-1 || r >= rows-1 || fillGrid[r * cols + c] == 0)
				return true;
		}
	}
	return false;
}

/**
 * Feathers the edges of data regions in a raster grid by setting the alpha
 * value for a pixel in proportion to its distance from the nearest null.
 */
void feather(float *srcGrid, float *dstGrid, int cols, int rows, float distance, float nodata, float resolution) {
	// Fill grid is used to keep track of where the edges are as they're "snowed in"
	// Starts out as a mask of non-nodata pixels from the source.
	char *fillGrid = (char *) malloc(cols * rows * sizeof(char));
	for(int i = 0; i < rows * cols; ++i)
		fillGrid[i] = srcGrid[i] == nodata ? 0 : 1;
	// The number of steps is just the number of pixels needed to 
	// cover the distance of the fade.
	float step = 0.0;
	float steps = _min(1.0, distance / resolution);
	bool found = false;
	// "Snow in" the alpha mask. The loop will exit when no more edge pixels can be found.
	do {
		found = false;
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				if(isEdgePixel(fillGrid, col, row, cols, rows)) {
					fillGrid[row * cols + col] = 2; // Set edge to dirty.
					dstGrid[row * cols + col] = _tanCurve(step, steps); // TODO: Configurable curves.
					found = true;
				}
			}
		}
		// Reset dirty edges to 0
		for(int i = 0; i < rows * cols; ++i) 
			if(fillGrid[i] == 2) fillGrid[i] = 0;
		++step;
	} while(found);
	free(fillGrid);
}

/**
 * Blends two rasters together using the alpha grid for blending.
 */
void blend(float *imgGrid, float *bgGrid, float *alpha, int cols, int rows, float imNodata, float bgNodata) {
	for(int i = 0; i < cols * rows; ++i) {
		if(!(bgGrid[i] == bgNodata || imgGrid[i] == imNodata))
			bgGrid[i] = bgGrid[i] * (1.0 - alpha[i]) + imgGrid[i] * alpha[i];
	}
}

/**
 * Mosaic the given files together using the first as the base. The base file will serve as a clipping
 * mask for the others. The spatial reference system and resolution of all layers must match.
 * The distance argument determines the distance over which the edges are feathered. The
 * overviews argument, if true, forces the construction of new overviews. If this is not
 * used, the existing overviews may obscure the fact that the image has changed (when zoomed out.)
 */
void mosaic(std::vector<std::string> &files, std::string &outfile, float distance, bool overviews) {
	
	if(distance <= 0.0)
		throw "Invalid distance.";

	if(outfile.size() == 0)
		throw "No output file given.";

	if(files.size() < 2)
		throw "Less than 2 files. Nothing to do.";

	GDALAllRegister();

	// Copy the background file to the destination.
	std::cout << "Copying background file." << std::endl;
	std::ifstream src(files[0].c_str(), std::ios::binary);
	std::ofstream dst(outfile.c_str(), std::ios::binary);
	dst << src.rdbuf();

	GDALDataset *outDS = NULL, *imDS = NULL;
	float *imGrid = NULL, *outGrid = NULL, *alphaGrid = NULL;
	double outTrans[6], imTrans[6];
	float imNodata, outNodata;
	int cols, rows, col, row;
	int rowHeight, rowOffset, bufRows;
	int addrOffset, bufRows0, bufRow0, rowHeight0, rowOffset0;

	try {

		// Open the destination file to modify it.
		outDS = (GDALDataset *) GDALOpen(outfile.c_str(), GA_Update);
		if(outDS == NULL)
			throw "Failed to open destination image.";

		// Load the transform parameters.
		outDS->GetGeoTransform(outTrans);

		// Iterate over the files, adding each one to the background.
		for(unsigned int i = 1; i < files.size(); ++i) {

			std::cout << "Processing file: " << files[i] << std::endl;

			imDS = (GDALDataset *) GDALOpen(files[i].c_str(), GA_ReadOnly);
			if(!imDS)
				throw "Failed to open mosaic image.";

			imDS->GetGeoTransform(imTrans);
			if(outTrans[1] != imTrans[1] || outTrans[5] != imTrans[5]) {
				GDALClose(imDS);
				throw "Raster's resolution does not match the background.";
			}

			cols = imDS->GetRasterXSize();
			rows = imDS->GetRasterYSize();
			col = (int) ((imTrans[0] - outTrans[0]) / outTrans[1]);
			row = (int) ((imTrans[3] - outTrans[3]) / outTrans[5]);
			
			if(col + cols > outDS->GetRasterXSize())
				cols = outDS->GetRasterXSize() - col;
			if(row + rows > outDS->GetRasterYSize())
				rows = outDS->GetRasterYSize() - row;

			imNodata = imDS->GetRasterBand(1)->GetNoDataValue();
			outNodata = outDS->GetRasterBand(1)->GetNoDataValue();
			
			// TODO: Configurable row height.
			rowHeight = 500 > rows ? rows : 500;
			rowOffset = (int) (distance / imTrans[1]) + 1;
			bufRows = rowHeight + rowOffset * 2;

			// Initialize the grids.
			imGrid = (float *) malloc(cols * bufRows * sizeof(float));
			alphaGrid = (float *) calloc(cols * bufRows, sizeof(float));
			outGrid = (float *) malloc(cols * bufRows * sizeof(float));

			// Move the window down by rowHeight and one rowOffset *not* two.
			for(int bufRow = -rowOffset; bufRow < rows; bufRow += rowHeight) { 

				// Set to zero or nodata depending on the band.
				std::fill_n(outGrid, cols * bufRows, outNodata);
				std::fill_n(imGrid, cols * bufRows, imNodata);
				std::fill_n(alphaGrid, cols * bufRows, 1.0);

				addrOffset = 0;
				bufRows0 = bufRows;
				bufRow0 = bufRow;
				rowHeight0 = rowHeight;
				rowOffset0 = rowOffset;
				if(bufRow0 < 0) {
					bufRow0 = 0;
					addrOffset = rowOffset * cols;
					bufRows0 = bufRows - rowOffset;
					rowOffset0 = 0;
				}

				// The last row might have to be smaller.
				if(bufRow0 + bufRows0 > rows) {
					bufRows0 = rows - bufRow0;
					rowHeight0 = bufRows0 - rowOffset0;
				}

				std::cout << bufRow0 << " - " << bufRows0 << " - " << rowHeight0 << " - " << rows << std::endl;

				// Load the overlay.
				CPLErr err = imDS->RasterIO(GF_Read, 0, bufRow0, cols, bufRows0, (imGrid + addrOffset), cols, bufRows0,
					GDT_Float32, 1, NULL, 0, 0, 0, NULL);
				if(err != CPLE_None)
					throw "Failed to read raster.";

				std::cout << "Feathering" << std::endl;
				// Feather the overlay
				feather(imGrid, alphaGrid, cols, bufRows, distance, imNodata, imTrans[1]);

				// Read background data.
				err = outDS->RasterIO(GF_Read, col, row + bufRow0, cols, bufRows0, (outGrid + addrOffset), cols, bufRows0,
					GDT_Float32, 1, NULL, 0, 0, 0, NULL);
				if(err != CPLE_None)
					throw "Failed to read raster.";

				std::cout << "Blending" << std::endl;
				// Blend the overlay into the output.
				blend(imGrid, outGrid, alphaGrid, cols, bufRows0, imNodata, outNodata);

				// Write back to the output.
				// We are extracting a slice out of the buffer, not writing the whole thing.
				err = outDS->RasterIO(GF_Write, col, row + bufRow0 + rowOffset0, cols, rowHeight0, (outGrid + addrOffset + rowOffset0 * cols), cols, rowHeight0,
					GDT_Float32, 1, NULL, 0, 0, 0, NULL);
				if(err != CPLE_None)
					throw "Failed to write raster.";
			}

			free(imGrid);
			free(outGrid);
			free(alphaGrid);

			GDALClose(imDS);

			if(overviews) {
				std::cout << "Building overviews..." << std::endl;
				// Must close and open the DS to avoid errors.
				GDALClose(outDS);
				outDS = (GDALDataset *) GDALOpen(outfile.c_str(), GA_Update);
				if(outDS == NULL)
					throw "Failed to open destination image for overview creation.";
				// TODO: Configure overviews and method.
				int overViews[] = { 2, 4, 8, 16 };
				outDS->BuildOverviews("NEAREST", 4, overViews, 0, NULL, NULL, NULL);
			}
		}
	} catch(const char *e) {

		free(imGrid);
		free(outGrid);
		free(alphaGrid);

		GDALClose(outDS);
		GDALClose(imDS);
		throw e;
	}

	GDALClose(outDS);

}

void usage() {
	std::cout << "Usage: mosaic [options] -o <output file> <file [file [file [...]]]>" << std::endl;
	std::cout << "    -o <file>     -- The output file." << std::endl;
	std::cout << "    -d <distance> -- The feather distance in map units (default 100.)" << std::endl;
	std::cout << "    -b            -- Build overviews (default: old overviews are preserved.)" << std::endl;
	std::cout << "    <file [...]>  -- A list of files. The first is used as the background " << std::endl;
	std::cout << "                     and determines the size of the output. Subsequent " << std::endl;
	std::cout << "                     images are layered on top and feathered." << std::endl;
}

int main(int argc, char **argv) {

 	try {

	 	float distance = 100.0;
	 	std::vector<std::string> files;
	 	std::string outfile;
	 	bool overviews = false;

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-d") {
	 			distance = atof(argv[++i]);
	 		} else if(arg == "-o") {
	 			outfile = argv[++i];
	 		} else if(arg == "-b") {
	 			overviews = true;
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

 		mosaic(files, outfile, distance, overviews);

 	} catch(const char *e) {
 		std::cerr << e << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }