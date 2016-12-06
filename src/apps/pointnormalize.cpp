#include <iostream>

#include "geotools.hpp"
#include "pointnormalize.hpp"
#ifdef WITH_GUI
#include "pointnormalize_ui.hpp"
#endif

using namespace geotools::point;

void usage() {
	std::cerr
			<< "Usage: pointnormalize [options] <output dir> <point file [point file [point file ...]]>\n"
			<< " -v                          Verbose output.\n"
			<< " -h                          Print this message.\n"
			<< " --threads                   The number of threads to use for computing output.\n";
}

int runWithUI(int argc, char **argv) {
#ifdef WITH_GUI
	QApplication q(argc, argv);
	QWidget *w = new QWidget();
	geotools::ui::PointNormalizeForm f;
	f.setupUi(w);
	w->show();
	return q.exec();
#else
	std::cerr << "GUI not enabled." << std::endl;
	return 1;
#endif
}

int main(int argc, char **argv) {

	try {

		PointNormalizeConfig config;

		g_loglevel(0);

		for (int i = 1; i < argc; ++i) {
			std::string s(argv[i]);
			if (s == "-h") {
				usage();
				return 0;
			} else if (s == "-v") {
				g_loglevel (G_LOG_DEBUG);
			} else if (s == "--threads") {
				config.threads = atoi(argv[++i]);
			} else {
				if (config.outputDir.empty()) {
					config.outputDir = s;
				} else {
					config.sourceFiles.push_back(s);
				}
			}
		}

		if (config.sourceFiles.size() == 0) {
			return runWithUI(argc, argv);
		} else {
			PointNormalize pn;
			pn.normalize(config);
		}

	} catch (const std::exception &ex) {
		g_error(ex.what());
		usage();
		return 1;
	}

	return 0;
}
