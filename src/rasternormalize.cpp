#include <iostream>
#include "raster.hpp"

using namespace geotools::raster;

int main(int argc, char **argv) {
	
	if(argc < 3) {
		std::cerr << "Usage: rasternormalize <input> <output>\n";
		return 1;
	}

	char *input = argv[1];
	char *output = argv[2];

	{
		Raster<float> inrast(std::string(input), false);
		{
			MemRaster<float> tmp(inrast.cols(), inrast.rows());
			tmp.setNodata(inrast.nodata());
			tmp.writeBlock(inrast);
			tmp.normalize();
			Raster<float> outrast(std::string(output), 1, inrast);
			outrast.setNodata(inrast.nodata());
			outrast.writeBlock(tmp);
		}
	}
	return 0;
}