
#include <set>
#include <vector>

#include "pointstats.hpp"

#ifdef WITH_GUI
#include "pointstats_ui.hpp"
#endif

using namespace geotools::util;
using namespace geotools::point;

void usage() {
    std::cerr << "Usage: pointstats <options> <file [file [file]]>\n"
            << " -o <output file>\n"
            << " -t <type>                   Output median, mean, max, min, variance (sample), pvariance (population),\n"
            << "                             count, density, stddev (sample), pstddev (population). Default mean.\n"
            << "                             Use a Comma-delimited list to do more than one. Requires a Comma-delimited\n"
            << "                             list of filenames.\n"
            << " -g <type>                   Gap Fraction type. These are IR, FR, RR, BLa and BLb, CCF and GAP.\n"
            << "                             the first five are adapted from Hopkins and Chasmer, 2009: Testing LiDAR Models "
            << "                             of Fractional Cover. The last two require a threshold value (-gt).\n"
            << " -gt <threshold>             A double representing the threshold height for gap fraction.\n"
            << " -r <resolution>             Resolution (default 2).\n"
            << " -s <srid>                   The EPSG ID of the CRS.\n"
            << " -c <classes>                Comma-delimited (e.g. '2,0' (ground and unclassified)).\n"
            << " -a <attribute>              Use height, intensity (default height).\n"
            << " -b <minx miny maxx maxy>    Extract points from the given box and create a raster of this size.\n"
            //		<< " -p                          Snap to the resolution.\n"
            << " -n                          Normalize the output so that one std. dev is represented as +-1.\n"
            << " -v                          Verbose output.\n"
            << " -h                          Print this message.\n"
            << " --threads                   The number of threads to use for computing output.\n"
            << " --angle-limit               Points located outside of this angle (devation from nadir) are excluded.\n"
            << " -gui                        Run the graphical user interface.\n";
}

int runWithUI(int argc, char **argv) {
#ifdef WITH_GUI
    QApplication q(argc, argv);
    QWidget *w = new QWidget();
    geotools::ui::PointStatsForm f;
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

        int crs = 0;
        int threads = 1;
        bool fill = false;
        bool gui = false;
        bool snap = true;
        bool rebuild = true;
        bool normalize = false;
        double resolution = 2.0;
        unsigned char angleLimit = 100;
        std::vector<std::string> dstFiles;
        std::vector<std::string> files;
        std::vector<std::string> types;
        std::string att = "height";
        std::string gap;
        std::set<unsigned char> classes;
        Bounds bounds;
        double gapThreshold;
        
        g_loglevel(0);

        for (int i = 1; i < argc; ++i) {
            std::string s(argv[i]);
            if (s == "-h") {
                usage();
                return 0;
            } else if (s == "-gui") {
                gui = true;
            } else if (s == "-o") {
                Util::splitString(argv[++i], dstFiles);
            } else if (s == "-s") {
                crs = atoi(argv[++i]);
            } else if (s == "-f") {
                fill = true;
            } else if (s == "-t") {
                Util::splitString(argv[++i], types);
            } else if (s == "-n") {
                normalize = true;
            } else if (s == "-r") {
                resolution = atof(argv[++i]);
            } else if (s == "-c") {
                Util::intSplit(classes, argv[++i]);
            } else if (s == "-a") {
                att = argv[++i];
            } else if (s == "-p") {
                snap = true;
            } else if (s == "-v") {
                g_loglevel(G_LOG_DEBUG);
            } else if (s == "-g") {
                gap = argv[++i];
            } else if(s == "-gt") {
                gapThreshold = atof(argv[++i]);
            } else if (s == "-C") {
                rebuild = false;
            } else if (s == "--angle-limit") {
                angleLimit = (unsigned char) atoi(argv[++i]);
            } else if (s == "--threads") {
                threads = atoi(argv[++i]);
            } else if (s == "-b") {
                bounds.extend(atof(argv[i + 1]), atof(argv[i + 2]));
                bounds.extend(atof(argv[i + 3]), atof(argv[i + 4]));
                i += 4;
            } else {
                files.push_back(argv[i]);
            }
        }

        if (gui) {
            return runWithUI(argc, argv);
        } else {
            PointStats lg;
            PointStatsConfig config;
            config.dstFiles = dstFiles;
            config.sourceFiles = files;
            config.classes = classes;
            config.hsrid = crs;
            config.attribute = config.parseAtt(att);
            config.types = config.parseTypes(types);
            config.resolution = resolution;
            config.bounds = bounds;
            config.angleLimit = angleLimit;
            config.fill = fill;
            config.snap = snap;
            config.threads = threads;
            config.gapFractionType = config.parseGap(gap);
            config.rebuild = rebuild;
            config.normalize = normalize;
            config.threshold = gapThreshold;
            lg.pointstats(config);
        }

    } catch (const std::exception &ex) {
        g_error(ex.what());
        usage();
        return 1;
    }

    return 0;
}
