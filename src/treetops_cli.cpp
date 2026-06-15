#include <iostream>
#include <string>
#include <vector>

#include <gdal/gdal_priv.h>

#include "config.hpp"
#include "process.hpp"
#include "util.hpp"

namespace {

void printUsage(const char* prog) {
	std::cerr
		<< "Usage: " << prog << " [options]\n"
		<< "\n"
		<< "Run treetop detection and crown delineation on a canopy height model.\n"
		<< "\n"
		<< "Options:\n"
		<< "  -c, --config FILE     JSON settings file\n"
		<< "  -i, --input FILE      Input CHM/DSM raster (overrides config)\n"
		<< "  -o, --output-dir DIR  Output directory (derives output paths from input name)\n"
		<< "      --no-smooth        Skip Gaussian smoothing\n"
		<< "      --no-crowns        Detect treetops only\n"
		<< "  -h, --help            Show this help\n"
		<< "\n"
		<< "Example:\n"
		<< "  " << prog << " -c _data/settings.json\n"
		<< "  " << prog << " -i _data/J5_10cm_CHM.tif -o _data\n";
}

struct Options {
	std::string configFile;
	std::string inputFile;
	std::string outputDir;
	bool noSmooth = false;
	bool noCrowns = false;
	bool showHelp = false;
};

Options parseArgs(int argc, char** argv) {
	Options opts;
	for(int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-h" || arg == "--help") {
			opts.showHelp = true;
		} else if((arg == "-c" || arg == "--config") && i + 1 < argc) {
			opts.configFile = argv[++i];
		} else if((arg == "-i" || arg == "--input") && i + 1 < argc) {
			opts.inputFile = argv[++i];
		} else if((arg == "-o" || arg == "--output-dir") && i + 1 < argc) {
			opts.outputDir = argv[++i];
		} else if(arg == "--no-smooth") {
			opts.noSmooth = true;
		} else if(arg == "--no-crowns") {
			opts.noCrowns = true;
		} else {
			throw std::runtime_error("Unknown argument: " + arg);
		}
	}
	return opts;
}

} // namespace

int main(int argc, char** argv) {
	GDALAllRegister();

	try {
		Options opts = parseArgs(argc, argv);
		if(opts.showHelp) {
			printUsage(argv[0]);
			return 0;
		}

		tt::config::Config config;

		if(!opts.configFile.empty())
			config.load(opts.configFile);

		if(!opts.inputFile.empty()) {
			std::string outputDir = opts.outputDir.empty()
				? tt::util::parent(opts.inputFile)
				: opts.outputDir;
			config.deriveOutputPaths(opts.inputFile);
			if(!opts.outputDir.empty()) {
				std::string stem = tt::util::basename(opts.inputFile);
				config.set("smoothedCHM", tt::util::join(outputDir, stem + "_smooth.tif"));
				config.set("treetopsDatabase", tt::util::join(outputDir, "tops.sqlite"));
				config.set("crownsRaster", tt::util::join(outputDir, stem + "_crowns.tif"));
				config.set("crownsDatabase", tt::util::join(outputDir, stem + "_crowns.sqlite"));
				config.set("topsWindowsRaster", tt::util::join(outputDir, "tops_windows.tif"));
				config.set("topsIdsRaster", tt::util::join(outputDir, "tops_ids.tif"));
			}
		} else if(!opts.outputDir.empty()) {
			std::string chm = config.get("originalCHM", "");
			if(chm.empty())
				throw std::runtime_error("No input raster. Use --input or set originalCHM in the config file.");
			config.deriveOutputPaths(chm);
		}

		if(opts.noSmooth)
			config.set("doSmoothing", false);
		if(opts.noCrowns)
			config.set("doCrowns", false);

		if(!config.canRun())
			throw std::runtime_error("No input raster configured. Use --input or set originalCHM in the config file.");

		std::cout << "Input:  " << config.get("originalCHM", "") << "\n";
		std::cout << "Smooth: " << (config.get("doSmoothing", false) ? config.get("smoothedCHM", "") : "(skipped)") << "\n";
		std::cout << "Tops:   " << config.get("treetopsDatabase", "") << "\n";
		if(config.get("doCrowns", true)) {
			std::cout << "Crowns: " << config.get("crownsRaster", "") << "\n";
			std::cout << "Poly:   " << config.get("crownsDatabase", "") << "\n";
		}

		tt::proc::Processor processor(config);
		processor.run();

		std::cout << "Done.\n";
		return 0;
	} catch(const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}
}
