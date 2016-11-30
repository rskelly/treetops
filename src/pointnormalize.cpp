#include <atomic>
#include <thread>
#include <vector>

#include <boost/filesystem.hpp>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <liblas/liblas.hpp>

#include "pointnormalize.hpp"
#include "raster.hpp"
#include "util.hpp"
#include "lasreader.hpp"

using namespace geotools::point;
using namespace geotools::raster;
using namespace geotools::las;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::All_faces_iterator All_faces_iterator;
typedef Delaunay::Face Face;
typedef Delaunay::Face_handle Face_handle;


PointNormalizeConfig::PointNormalizeConfig() :
    dropNegative(true),
    dropGround(true),
    threads(1),
    overwrite(false),
    buffer(10.0) {
}

class FileSorter {
private:
    double m_colSize;
    double m_rowSize;
public:

    FileSorter(double colSize, double rowSize) :
    m_colSize(colSize), m_rowSize(rowSize) {
    }

    bool operator()(const std::string &a, const std::string &b) {
        LASReader ar(a);
        LASReader br(b);
        Bounds ab = ar.bounds();
        Bounds bb = br.bounds();
        int idxa = ((int) (ab.miny() / m_rowSize)) * ((int) (ab.width() / m_colSize)) + ((int) (ab.minx() / m_colSize));
        int idxb = ((int) (bb.miny() / m_rowSize)) * ((int) (bb.width() / m_colSize)) + ((int) (bb.minx() / m_colSize));
        return idxa < idxb;
    }
};

class FileSearcher {
private:
    std::vector<std::string> m_files;
    std::vector<std::pair<std::string, Bounds> > m_bounds;
public:
    FileSearcher(const std::vector<std::string> &files) :
        m_files(files) {
            init();
    }
        
    void init() {
        for(const std::string &file : m_files) {
            LASReader r(file);
            Bounds b = r.bounds();
            m_bounds.push_back(std::make_pair(file, b));
        }
    }
    
    std::vector<std::string> getFiles(const Bounds &bounds) {
        std::vector<std::string> files;
        for(const auto &pr : m_bounds) {
            std::cerr << std::setprecision(9);
            g_debug(" -- check: " << pr.second.print() << "; " << bounds.print() << "; " << pr.second.intersects(bounds));
            if(pr.second.intersects(bounds))
                files.push_back(pr.first);
        }
        return files;
    }
};

bool _cancel = false;

void PointNormalize::normalize(const PointNormalizeConfig &config, 
        const Callbacks *callbacks, bool *cancel) {

    // If the cancel flag isn't given, take a reference to a dummy
    // that is always false.
    if(!cancel) cancel = &_cancel;
    
    // Sanity checks.
    if (config.outputDir.empty())
        g_argerr("A point output dir must be provided.");
    if (!Util::mkdir(config.outputDir))
        g_argerr("Couldn't create output dir " << config.outputDir);
    uint32_t threads = config.threads;
    uint32_t tmax = std::thread::hardware_concurrency();
    if(threads < 1) {
        g_warn("Too few threads specified. Using 1.");
        threads = 1;
    }
    if(threads > tmax && tmax > 0) {
        g_warn("Too many threads specified. Using 1.");
        threads = 1;
    }
    
    // Set the thread count.
    omp_set_dynamic(1);
    omp_set_num_threads(threads);
    
    using namespace boost::filesystem;

    path outDir(config.outputDir);
    liblas::ReaderFactory rf;

    std::vector<std::string> files(config.sourceFiles.begin(), config.sourceFiles.end());
    FileSearcher fileSearch(files);
    
    for (unsigned int i = 0; !*cancel && i < files.size(); ++i) {
        if(*cancel) return;
        
        if(callbacks)
            callbacks->overallCallback(((float) i + 0.5f) / files.size());
        
        const std::string &pointFile = files[i];
        path oldPath(pointFile);
        path newPath(outDir / oldPath.filename());

        g_debug(" -- normalize (point) processing " << pointFile << " to " << newPath);
        if (!config.overwrite && exists(newPath)) {
            g_warn("The output file " << newPath << " already exists.");
            continue;
        }

        std::vector<liblas::Point> objPts;
        std::list<Point_3> meshPts;
    
        // Build a mesh from the buffered ground points.
        LASReader lr(pointFile);
        Bounds bounds = lr.bounds();
        g_debug(" -- bounds: " << bounds.print());
        bounds.extend(bounds.minx() - 10.0, bounds.miny() - 10.0);
        bounds.extend(bounds.maxx() + 10.0, bounds.maxy() + 10.0);
        g_debug(" -- bounds: " << bounds.print());
        std::vector<std::string> bufFiles = fileSearch.getFiles(bounds);
        g_debug(" -- " << bufFiles.size() << " files to satisfy the buffer");
        for(const std::string &file : bufFiles) {            
            std::ifstream instr(file.c_str(), std::ios::binary);
            liblas::Reader lasReader = rf.CreateWithStream(instr);
            while(lasReader.ReadNextPoint()) {
                const liblas::Point &pt = lasReader.GetPoint();
                if(bounds.contains(pt.GetX(), pt.GetY()) && pt.GetClassification().GetClass() == 2) {
                    Point_3 p(pt.GetX(), pt.GetY(), pt.GetZ());
                    meshPts.push_back(std::move(p));                    
                }
                if(*cancel)
                    return;
            }
        }
        
        // Load the non-ground points.
        std::ifstream instr(pointFile.c_str(), std::ios::binary);
        liblas::Reader lasReader = rf.CreateWithStream(instr);
        liblas::Header lasHeader = lasReader.GetHeader();

        while (lasReader.ReadNextPoint()) {
            if(*cancel) return;
            const liblas::Point &opt = lasReader.GetPoint();
            if(!config.dropGround || opt.GetClassification().GetClass() != 2)
                objPts.push_back(opt);
            if(*cancel)
                return;
        }

        // Build the Delaunay triangulation.
        Delaunay dt(meshPts.begin(), meshPts.end());
        meshPts.clear();
        if(*cancel) return;
        
        Face_handle hint;
        uint64_t finalPtCount = 0;
        uint64_t finalPtCountByReturn[255];
        for(int i = 0; i < 5; ++i)
            finalPtCountByReturn[i] = 0;
        
        std::atomic<size_t> curPt(0);
        size_t ptCount = objPts.size();
        
        std::ofstream outstr(newPath.c_str(), std::ios::binary);
        liblas::Header outHeader(lasHeader);
        liblas::Writer lasWriter(outstr, outHeader);

        #pragma omp parallel for
        for(size_t i = 0; i < ptCount; ++i) {
            if(*cancel) continue;
            
            if(callbacks)
                callbacks->stepCallback((float) ++curPt / ptCount);
            
            liblas::Point &opt = objPts[i];
            if(opt.GetClassification().GetClass() == 2) {
                opt.SetZ(0.0);
            } else {
                Point_3 p(opt.GetX(), opt.GetY(), opt.GetZ());
                hint = dt.locate(p, hint);
                // TODO: The point is outside the mesh and can't be normalized.
                //       Ideally, this would build the mesh to the neighbouring
                //       files and normalize these points...
                if(dt.is_infinite(hint))
                    continue;
                double area = 0.0;
                double total = 0.0;
                for(int i = 0; i < 3; ++i) {
                    Point_3 p1 = hint->vertex(i)->point();
                    Point_3 p2 = hint->vertex((i + 1) % 3)->point();
                    Point_3 p3 = hint->vertex((i + 2) % 3)->point();
                    double h = Util::computeArea(
                        p1.x(), p1.y(), p1.z(), 
                        p2.x(), p2.y(), p2.z(), 
                        p.x(), p.y(), p.z()
                    );
                    area += h;
                    total += h * p3.z();
                }
                double z = p.z() - total / area;
                if(z < 0.0 && config.dropNegative)
                    continue;
                opt.SetZ(z);
            }
            #pragma omp critical
            {
                lasWriter.WritePoint(opt);
                finalPtCountByReturn[opt.GetReturnNumber() - 1]++;
                ++finalPtCount;
            }
        }

        g_debug(" -- setting original count " << outHeader.GetPointRecordsCount() << " to " << ptCount);
        outHeader.SetPointRecordsCount(finalPtCount);
        for(int i = 0; i < 5; ++i)
            outHeader.SetPointRecordsByReturnCount(i + 1, finalPtCountByReturn[i]);
        lasWriter.SetHeader(outHeader);
        
        instr.close();
        outstr.close();

        if(callbacks)
            callbacks->overallCallback(((float) i + 1.0f) / files.size());
    }


}