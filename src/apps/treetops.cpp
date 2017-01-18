#include <string>
#include <vector>
#include <map>

#include "omp.h"

#include "geotools.hpp"
#include "treetops.hpp"

#ifdef WITH_GUI
#include "treetops_ui.hpp"
#include "cpl_conv.h"
#endif

#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")

void usage() {
	std::cerr << "Usage: treetops <options>\n"
			<< "This program finds tree tops by locating the maximum value in a \n"
			<< "window as it moves across a raster. The tops are just the maxima at \n"
			<< "the center of the window.\n"
			<< " -so <filename>     The input raster to be smoothed.\n"
			<< " -sd <int>          The standard deviation for Gaussian smoothing.\n"
			<< " -sw <int>          The window size for Gaussian smoothing. Will be \n"
			<< "                    bumped up to the next odd value if even given. Minimum 3.\n"
			<< " -ss <filename>     The smoothed output raster.\n\n"

			<< " -to <filename>     The input raster for treetop detection; a (non-smoothed)\n"
			<< "                    LiDAR-derived canopy height model.\n"
			<< " -ts <filename>     The input raster for treetop detection; a (smoothed)\n"
			<< "                    LiDAR-derived canopy height model.\n"
			<< " -tt <float int>    Threshold. A space-delimited pair consisting of a float and a positive\n"
			<< "                    byte representing a maximum height and window size for locating tops.\n"
			<< "                    Window size will be bumped up to the next odd value; minumum 3.\n"
			<< " -td <filename>     The treetop vector file. An sqlite file.\n\n"

			<< " -cs <filename>     The input raster for crown delineation; a (smoothed)\n"
			<< "                    LiDAR-derived canopy height model.\n"
			<< " -ct <filename>     The treetop vector file. An sqlite file.\n\n"
			<< " -cm <float>        The minimum height of pixels to consider for inclusion \n"
			<< "                    in a crown. Default 4.\n"
			<< " -cf <float>        The minimum height of pixels, as a proportion of the \n"
			<< "                    treetop height to consider for inclusion in a crown.\n"
			<< "                    0 < n < 1; default 0.65.\n"
			<< " -cd <float>        The radius beyond which pixels will not be considered for \n"
			<< "                    inclusion in a crown.\n"
			<< " -cr <filename>     The crowns raster. Each crown is represented with pixels\n"
			<< "                    having a unique ID.\n"
			<< " -cp <filename>     The crowns vector. A vectorized version of the crown \n"
			<< "                    raster. Optional.\n\n"

			<< " -threads           The number of threads to used for processing. Default 1.\n"
			<< " -v                 Print debug messages.\n";
}

int runWithGui(int argc, char **argv) {
#ifdef WITH_GUI

#ifdef _MSC_VER
	CPLSetConfigOption("GDAL_DATA", "./gdal-data");
#endif

	class TTApplication : public QApplication {
	public:
		TTApplication(int &argc, char **argv) : QApplication(argc, argv) {}
		bool notify(QObject *receiver, QEvent *e) {
			try {
				return QApplication::notify(receiver, e);
			} catch(const std::exception &ex) {
				QMessageBox err;
				err.setText("Error");
				err.setInformativeText(QString(ex.what()));
				err.exec();
				return false;
			}
		}
	};

	TTApplication q(argc, argv);
	geotools::ui::TreetopsForm f;
	f.show();
	return q.exec();
#else
	std::cerr << "GUI not enabled." << std::endl;
	return 1;
#endif
}

int main(int argc, char **argv) {

	using namespace geotools::treetops;
	using namespace geotools::treetops::util;
	using namespace geotools::treetops::config;

	try {

		TreetopsConfig config;

		int i = 1;
		for (; i < argc; ++i) {
			std::string arg = argv[i];

			if (arg == "-so") {
				config.smoothOriginalCHM = argv[++i];
				config.doSmoothing = true;
			} else if (arg == "-sd") {
				config.smoothSigma = atof(argv[++i]);
				config.doSmoothing = true;
			} else if (arg == "-ss") {
				config.smoothSmoothedCHM = argv[++i];
				config.doSmoothing = true;
			} else if (arg == "-to") {
				config.topsOriginalCHM = argv[++i];
				config.doTops = true;
			} else if (arg == "-ts") {
				config.topsSmoothedCHM = argv[++i];
				config.doTops = true;
			} else if (arg == "-tt") {
				if(argc < i + 2)
					g_argerr("Too few parameters for -tt.");
				float height = atof(argv[++i]);
				uint8_t window = atoi(argv[++i]);
				config.topsThresholds[height] = window;
				config.doTops = true;
			} else if (arg == "-td") {
				config.topsTreetopsDatabase = argv[++i];
				config.doTops = true;
			} else if (arg == "-cs") {
				config.crownsSmoothedCHM = argv[++i];
				config.doCrowns = true;
			} else if (arg == "-ct") {
				config.crownsTreetopsDatabase = argv[++i];
				config.doCrowns = true;
			} else if (arg == "-cm") {
				config.crownsMinHeight = atof(argv[++i]);
				config.doCrowns = true;
			} else if (arg == "-cf") {
				config.crownsHeightFraction = atof(argv[++i]);
				config.doCrowns = true;
			} else if (arg == "-cd") {
				config.crownsRadius = atof(argv[++i]);
				config.doCrowns = true;
			} else if (arg == "-cr") {
				config.crownsCrownsRaster = argv[++i];
				config.doCrowns = true;
			} else if (arg == "-cp") {
				config.crownsCrownsDatabase = argv[++i];
				config.doCrowns = true;
			} else if (arg == "-threads") {
				config.threads = atoi(argv[++i]);
				config.doCrowns = true;
			} else if (arg == "-v") {
				g_loglevel(0);
			}
		}

		Treetops tu;

		// Create tree tops.
		if (config.doSmoothing)
			tu.smooth(config);

		// Create treetops.
		if (config.doTops)
			tu.treetops(config);

		// Create crowns.
		if (config.doCrowns)
			tu.treecrowns(config);

		if (!config.doTops && !config.doSmoothing && !config.doCrowns)
			return runWithGui(argc, argv);

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}

