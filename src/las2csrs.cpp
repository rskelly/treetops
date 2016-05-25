#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <liblas/liblas.hpp>
#include <proj_api.h>
#include <gdal_priv.h>
#include <gdal.h>
#include <cpl_conv.h> // for CPLMalloc()

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

#define LAS2CSRS_DATA "LAS2CSRS_DATA"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
namespace las = liblas;

bool quiet = true;

// Binary interpolation.
static double _binterp(float *grid, double c, double r, int c0, int r0, int c1, int r1, int width) {
	double x1 = (c1 - c) / (c1 - c0) * grid[r0 * width + c0] + (c - c0) / (c1 - c0) * grid[r0 * width + c1];
	double x2 = (c1 - c) / (c1 - c0) * grid[r1 * width + c0] + (c - c0) / (c1 - c0) * grid[r1 * width + c1];
	return (r1 - r) / (r1 - r0) * x1 + (r - r0) / (r1 - r0) * x2;
}

/**
 * Convert projected distance in mm to latlon in radians.
 * dx, dy, dx  -- Distance in m.
 * lat         -- The latitude at which distances are computed
 * h           -- Ellipsoidal height.
 * a           -- Semi-major axis
 * e2          -- Eccentricity^2
 * count       -- Number of points.
 * dlat        -- The distance in rad (out).
 * dlon        -- The distance in rad (out).
 */
static void _shift2latlon(Grid<double> &gdx, Grid<double> &gdy, Grid<double> &glat, Grid<double> &gh,
	double a, double e2, int count, Grid<double> &gdlat, Grid<double> &gdlon) {

	double *dx = gdx.grid();
	double *dy = gdy.grid();
	double *lat = glat.grid();
	double *dlat = gdlat.grid();
	double *dlon = gdlon.grid();
	double *h = gh.grid();

	for(int i = 0; i < count; ++i) {
		double l = *(lat + i);
		double m = a * (1.0 - e2) / pow((1.0 - e2 * _sq(sin(l))), 3.0/2.0); 	// Meridional radius of curvature.
		double n = a / pow((1.0 - e2 * _sq(sin(l))), 1.0/2.0); 					// Parallel radius of curvature.
		double r = n * cos(l); 													// Radius of parallel.
		*(dlon + i) = *(dx + i) / (r + *(h + i));
		*(dlat + i) = *(dy + i) / (m + *(h + i));		
	}
}

// Miliarcseconds to radians.
inline double _mas2rad(double x) {
	return (x * 4.84813681 / 1000000000.0);
}


/**
 * Loads and interpolates the NAD83(CSRS) shift grid.
 */
class ShiftGrid {
private:
	MemRaster<float> xg;
	MemRaster<float> yg;
	MemRaster<float> zg;
	double tg[6];
	int width;
	int height;

public:

	/**
	 * Load the shift grid raster and store the x, y, z shifts internally.
	 */
	void load() {
		char path[PATH_MAX];
		sprintf(path, "%s%s", std::getenv("LAS2CSRS_DATA"), "/NAD83v6VG.tif");

		GDALAllRegister();
		GDALDataset *ds;
		GDALRasterBand *xband, *yband, *zband;

		ds = (GDALDataset *) GDALOpen(path, GA_ReadOnly);
		if(!ds)
			throw "Failed to load shift grid.";

		xband = ds->GetRasterBand(1);
		yband = ds->GetRasterBand(2);
		zband = ds->GetRasterBand(3);
		if(!xband || !yband || !zband)
			throw "Failed to retrieve shift bands.";

		width = xband->GetXSize();
		height = yband->GetYSize();
		if(width <= 0 || height <= 0)
			throw "The dimensions of the shift grid are invalid.";

		xg.init(width, height);
		yg.init(width, height);
		zg.init(width, height);

		CPLErr xe = xband->RasterIO(GF_Read, 0, 0, width, height, xg.grid(), width, height, GDT_Float32, 0, 0);
		CPLErr ye = yband->RasterIO(GF_Read, 0, 0, width, height, yg.grid(), width, height, GDT_Float32, 0, 0);
		CPLErr ze = zband->RasterIO(GF_Read, 0, 0, width, height, zg.grid(), width, height, GDT_Float32, 0, 0);

		if(xe == CE_Failure || ye == CE_Failure || ze == CE_Failure)
			throw "Failed to read shift grid.";

		ds->GetGeoTransform(tg);
	}

	/**
	 * Compute the the shifts in x, y and z at position x, y. x and y are radians,
	 * dx, dy and dz are m.
	 */
	void interpolate(Grid<double> &gx, Grid<double> &gy, Grid<double> &gdx, 
		Grid<double> &gdy, Grid<double> &gdz) {

		int count = gx.size();
		double *x = gx.grid();
		double *y = gy.grid();
		double *dx = gdx.grid();
		double *dy = gdy.grid();
		double *dz = gdz.grid();

		double c, r;
		int c0, c1, r0, r1;
		for(;count;--count, ++x, ++y, ++dx, ++dy, ++dz) {
			c = ((_deg(*x) - tg[0]) / tg[1]);
			r = ((_deg(*y) - tg[3]) / tg[5]);
			c0 = (int) c;
			r0 = (int) r;
			c1 = c0 + 1;
			r1 = r0 + 1;
			if(c0 < 0) c0 = 0;
			if(r0 < 0) r0 = 0;
			if(c1 >= width) c1 = width-1;
			if(r1 >= height) r1 = height-1;
			*dx = _binterp(xg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0; // Convert to m
			*dy = _binterp(yg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
			*dz = _binterp(zg.grid(), c, r, c0, r0, c1, r1, width) / 1000.0;
		}
	}

	~ShiftGrid() {
	}
};

class Params {
public:
	// Source reference frame.
	std::string ffrom;
	// From and to epochs.
	double efrom;
	double eto;
	// SRIDs of the from and to CRSes.
	int fsrid;
	int tsrid;
	std::string fgeoid;
	std::string tgeoid;
	// Transform parameters; loaded from the itrf file.
	double tx, ty, tz, dtx, dty, dtz;		// Shifts, rates.
	double rx, ry, rz, drx, dry, drz;		// Rotations, rates.
	double epoch;							// ITRF Transform epoch.
	//double dt; 							// Time delta
	double d, dd; 							// Scale, scale rate.

};

/**
 * Performs the work of transforming coordinates from a given reference frame to NAD83(CSRS).
 */
class Transformer {
private:
	// Projection objects.
	projPJ projFrom;
	projPJ projECEF;
	projPJ projTo;
	projPJ projGeog;
	// The shift grid.
	ShiftGrid shiftGrid;
	// Transformation params;
	Params params;     

	/**
		Transform the coordinate using the procedure listed in Craymer (2006).

		x, y, z -- 	The coordinate arrays.
		count 	-- 	The number of points.
	*/
	void epochTransform(Params &params, Grid<double> &gx, Grid<double> &gy, Grid<double> &gz, int count, double dt) {

		if(!quiet)
			std::cerr << "epochTransform" << std::endl;

		double *x = gx.grid();
		double *y = gy.grid();
		double *z = gz.grid();

		double Txt = params.tx + params.dtx * dt;            // Translation in X plus X velocity * time
		double Tyt = params.ty + params.dty * dt;
		double Tzt = params.tz + params.dtz * dt;
		double DSt = params.d + params.dd * dt;         	// Scale plus delta scale * time
		double Rxt = _mas2rad(params.rx + params.drx * dt); // Rotation in X plus velocity * time; seconds to radians.
		double Ryt = _mas2rad(params.ry + params.dry * dt);
		double Rzt = _mas2rad(params.rz + params.drz * dt);

		DSt += 1.0; // Scale is w/r/t 1.

		if(!quiet) {
			std::cerr << "  > Txt: " << Txt << "; Tyt: " << Tyt << "; Tzt: " << Tzt << std::endl;
			std::cerr << "  > Rxt: " << Rxt << "; Ryt: " << Ryt << "; Rzt: " << Rzt << std::endl;
			std::cerr << "  > DSt: " << DSt << std::endl;
			std::cerr << "  > drx: " << params.drx << "; dry: " << params.dry << "; drz: " << params.drz << std::endl;
		}

		for(;count;--count, ++x, ++y, ++z) {
			*x = Txt + (DSt * *x) + (-Rzt * *y) + (Ryt * *z);
			*y = Tyt + (Rzt * *x) + (DSt * *y) + (-Rxt * *z);
			*z = Tzt + (-Ryt * *x) + (Rxt * *y) + (DSt * *z);
		}
	}

	/**
	 * Initialize the proj4 projection objects.
	 */
	void initProjections() {

		if(!quiet)
			std::cerr << "initProjection" << std::endl;

		// Initialize projections.
		if(!params.fsrid || !params.tsrid)
			throw "SRIDs are not set.";
		char str[128];
		if(!params.fgeoid.empty()) {
			sprintf(str, "+init=epsg:%u +geoidgrids=%s", params.fsrid, params.fgeoid.c_str());
		} else {
			sprintf(str, "+init=epsg:%u", params.fsrid);
		}
		if(!quiet)
			std::cerr << "From projection: " << str << std::endl;
		projFrom = pj_init_plus(str);
		if(!projFrom)
			throw std::string(pj_strerrno(pj_errno));
		if(!params.tgeoid.empty()) {
			sprintf(str, "+init=epsg:%u +geoidgrids=%s", params.tsrid, params.tgeoid.c_str());
		} else {
			sprintf(str, "+init=epsg:%u", params.tsrid);
		}
		if(!quiet)
			std::cerr << "To projection: " << str << std::endl;
		projTo = pj_init_plus(str);
		if(!projTo)
			throw std::string(pj_strerrno(pj_errno));
		projECEF = pj_init_plus("+proj=geocent +ellps=GRS80 +units=m +no_defs");
		if(!projECEF)
			throw std::string(pj_strerrno(pj_errno));
		projGeog = pj_init_plus("+proj=latlon +ellps=GRS80 +no_defs");
		if(!projGeog)
			throw std::string(pj_strerrno(pj_errno));
	}

	/**
	 * Load the transformation database.
	 */
	void loadHelmert(Params &p) {

		if(!quiet)
			std::cerr << "loadHelmert " << p.ffrom << std::endl;

		char path[PATH_MAX];
		sprintf(path, "%s%s", std::getenv(LAS2CSRS_DATA), "/itrf.csv");

		char ffrom[64], fto[64];
		bool found = false;
		float epoch, tx, ty, tz, rx, ry, rz, d, dtx, dty, dtz, drx, dry, drz, dd;
		FILE *f = fopen(path, "r");
		if(f == NULL)
			throw "ITRF database file not found.";
		while(fscanf(f, " %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
						ffrom, fto, &epoch, &tx, &ty, &tz, &rx, &ry, &rz, &d,
						&dtx, &dty, &dtz, &drx, &dry, &drz, &dd) > 0) {
			std::string _ffrom(ffrom);
			if(!quiet)
				std::cerr << " -- Checking: " << _ffrom << std::endl;
			if(_ffrom == p.ffrom) {
				if(!quiet)
					std::cerr << " -- Found entry for " << _ffrom << std::endl;
				p.epoch = epoch;
				p.d = d / 1000000000.0; 		// Listed as ppb.
				p.dd = dd / 1000000000.0;
				p.tx = tx;
				p.ty = ty;
				p.tz = tz;
				p.rx = rx;
				p.ry = ry;
				p.rz = rz;
				p.dtx = dtx;
				p.dty = dty;
				p.dtz = dtz;
				p.drx = drx;
				p.dry = dry;
				p.drz = drz;
				found = true;
				break;
			}
		}
		fclose(f);

		if(!found)
			throw "Failed to find a transformation matching the parameters.";

		if(!quiet) {
			std::cerr << " -- Params: dtx: " << p.dtx << "; dty: " << p.dty << "; dtz: " << p.dtz << std::endl;
			std::cerr << "            drx: " << p.drx << "; dry: " << p.dry << "; drz: " << p.drz << std::endl;
		}

	}

public:

	/**
	 * Prepares the Helmert transformation parameters for the given transformation.
	 * These will be used in later method calls.
	 * The constructor will load an process the grid shift file and the transformation database.
	 * ffrom  -- The name of the reference frame, e.g. 'itrf90'
	 * efrom  -- The epoch of data collection (decimal years), e.g. 1994.2
	 * eto    -- The target epoch. One might select 1997.0 for BC or 2002.0 for Albera (etc.)
	 * fsrid  -- The SRID (EPSG code) of the source.
	 * tsrid  -- The SRID of the destination. This is the code for the UTM zone under
	 *           NAD83(CSRS), e.g. 2956 for UTM12N.
	 * fvsrid -- The input geoid model. This presumes input orthometric heights and the 
	 *           existence of a geoid model in the proj database. If the input is ellipsoidal, leave
	 *           as empty. (Example input: HT2_0.gtx for CGVD28)
	 * tvsrid -- The output geoid model. If ellipsoidal output is desired, leave as empty.
	 */
	Transformer(const std::string &ffrom, float efrom, float eto, int fsrid, int tsrid, const std::string &fgeoid, const std::string &tgeoid) {

		if(!quiet) {
			std::cerr << "Transformer:" << std::endl
				<< " -- src ref frame: " << ffrom << "; src epoch: " << efrom << "; dst epoch: " << eto << std::endl
				<< " -- src srid: " << fsrid << "; dst srid: " << tsrid << std::endl
				<< " -- src vsrid: " << fgeoid << "; dst vsrid: " << tgeoid << std::endl;
		}
		params.efrom = efrom;
		params.eto = eto;
		params.fsrid = fsrid;
		params.tsrid = tsrid;
		params.fgeoid.assign(fgeoid);
		params.tgeoid.assign(tgeoid);
		params.ffrom.assign(ffrom);

		initProjections();
		loadHelmert(params);
		shiftGrid.load();
	}

	/**
	 * Transforms coordinate(s) from one available reference frame to NAD83(CSRS).
	 * x, y, z -- Coordinate arrays.
	 * count   -- The number of coordinates.
	 * bounds  -- The bounds of the transformed coordinate list.
	 */
	void transformPoints(Grid<double> &x, Grid<double> &y, Grid<double> &z, int count, double bounds[6]) {

		// 1) Transform original points to their ECEF representation (Cartesian) at the original epoch.
		// 2) Deduce the transformation parameters to NAD83 @ 1997.0.
		// 3) Transform from NAD83 forward to target epoch.

		if(!quiet) {
			std::cerr << "transformPoints" << std::endl;
			std::cerr << " -- Original: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

			char str[128];
			sprintf(str, "+init=epsg:%u", params.fsrid);
			projPJ projTmp = pj_init_plus(str);

			MemRaster<double> tx(count, 1);
			MemRaster<double> ty(count, 1);
			MemRaster<double> tz(count, 1);
			memcpy(tx.grid(), x.grid(), sizeof(double) * count);
			memcpy(ty.grid(), y.grid(), sizeof(double) * count);
			memcpy(tz.grid(), z.grid(), sizeof(double) * count);
			pj_transform(projFrom, projTmp, count, 1, tx.grid(), ty.grid(), tz.grid());
			std::cerr << " -- Ellipsoidal: " << tx.grid()[0] << ", " << ty.grid()[0] << ", " << tz.grid()[0] << std::endl;
		}

		// 1) Project to ECEF in the original reference frame.
		pj_transform(projFrom, projECEF, count, 1, x.grid(), y.grid(), z.grid());

		if(!quiet)
			std::cerr << " -- ECEF (Original): " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

		// 2) Transform to NAD83 @ 1997.
		epochTransform(params, x, y, z, count, params.efrom - 1997.0);

		if(!quiet)
			std::cerr << " -- ECEF (CSRS): " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

		// Only use the grid shift if the epoch changes.
		if(params.efrom != params.eto) {

			// 3) Transform from NAD83 @ 1997 to NAD83 at target epoch.

			// Copy the coordinate arrays for transformation.
			MemRaster<double> x0(count, 1);
			MemRaster<double> y0(count, 1);
			MemRaster<double> z0(count, 1);
			memcpy(x0.grid(), x.grid(), sizeof(double) * count);
			memcpy(y0.grid(), y.grid(), sizeof(double) * count);
			memcpy(z0.grid(), z.grid(), sizeof(double) * count);

			// Transform from Cartesian (CSRS) to latlon.
			pj_transform(projECEF, projGeog, count, 1, x0.grid(), y0.grid(), z0.grid());

			if(!quiet)
				std::cerr << " -- Lat Lon (CSRS): " << _deg(x0.grid()[0]) << ", " << _deg(y0.grid()[0]) << ", " << z0.grid()[0] << std::endl;

			// Initalize shift arrays.
			MemRaster<double> dx(count, 1);
			MemRaster<double> dy(count, 1);
			MemRaster<double> dz(count, 1);
			
			// Interpolate shifts using latlon coords -- returns in mm. (d)
			shiftGrid.interpolate(x0, y0, dx, dy, dz);
			
			if(!quiet)
				std::cerr << " -- Grid Velocities: " << dx.grid()[0] << ", " << dy.grid()[0] << ", " << dz.grid()[0] << std::endl;

			// Transform mm shifts to latlon
			MemRaster<double> dlat(count, 1);
			MemRaster<double> dlon(count, 1);

			// Get projection's spheroid props.
			double a, e2;
			pj_get_spheroid_defn(projTo, &a, &e2);

			// Get angular shifts from grid shift (mm).
			_shift2latlon(dx, dy, y0, z0, a, e2, count, dlat, dlon);

			if(!quiet)
				std::cerr << " -- Lat lon shifts: " << dlat.grid()[0] << ", " << dlon.grid()[0] << std::endl;

			// Good to here...
			
			double dt = params.eto - params.efrom;
			// Apply shifts to latlon coords.
			for(int i = 0; i < count; ++i) {
				*(x0.grid() + i) += *(dlon.grid() + i) * dt;
				*(y0.grid() + i) += *(dlat.grid() + i) * dt;
				*(z0.grid() + i) += *(dz.grid() + i) * dt;
			}
			
			if(!quiet)
				std::cerr << " -- Lat lon shifted: " << _deg(x0.grid()[0]) << ", " << _deg(y0.grid()[0]) << ", " << z0.grid()[0] << std::endl;

			// Transform latlon to target proj
			pj_transform(projGeog, projTo, count, 1, x0.grid(), y0.grid(), z0.grid());
			
			// Assign the shifted coords to the output arrays.
			memcpy(x.grid(), x0.grid(), sizeof(double) * count);
			memcpy(y.grid(), y0.grid(), sizeof(double) * count);
			memcpy(z.grid(), z0.grid(), sizeof(double) * count);

		} else {

			// Reproject to the dest coordinates
			pj_transform(projECEF, projTo, count, 1, x.grid(), y.grid(), z.grid());

		}

		if(!quiet)
			std::cerr << std::setprecision(15) << " -- Final: " << x.grid()[0] << ", " << y.grid()[0] << ", " << z.grid()[0] << std::endl;

		// Expand the bounds for the new header.
		for(int i = 0; i < count; ++i) {
			if(x[i] < bounds[0]) bounds[0] = x[i];
			if(x[i] > bounds[1]) bounds[1] = x[i];
			if(y[i] < bounds[2]) bounds[2] = y[i];
			if(y[i] > bounds[3]) bounds[3] = y[i];
			if(z[i] < bounds[4]) bounds[4] = z[i];
			if(z[i] > bounds[5]) bounds[5] = z[i];
		}			
	}

	/**
	 * Transform the given LAS file(s) from the given reference frame to the one
	 * configured in this Transformer. 
	 * srcfile  -- The folder containing las files, or the path to a single las file.
	 * dstdir   -- The destination folder.
	 */
	int transformLas(std::string &srcfile, std::string &dstdir, bool overwrite) {
		
		if(!quiet)
			std::cerr << "transformLas" << std::endl;

		const fs::path src(srcfile);
		const fs::path dst(dstdir);

		// Check the files for sanity.
		if(src == dst)
			throw "Destination and source are the same: " + dst.string() + "==" + src.string();

		if(!fs::exists(dst)) {
			if(!fs::create_directory(dst))
				throw "Failed to create output directory: " + dst.string();
		}

		// Get the list of las files.
		if(!quiet)
			std::cerr << " -- Getting file list." << std::endl;

		std::list<fs::path> files;
		if(fs::is_regular_file(src)) {
			if(overwrite || !fs::exists(dst / src.leaf()))
				files.push_back(src);
		} else {
			std::string ext(".las");
			fs::directory_iterator end;
			fs::directory_iterator di(src);
			for(; di != end; ++di) {
				std::string p(di->path().string());
				alg::to_lower(p);
				if((overwrite || !fs::exists(dst / di->path().leaf())) && alg::ends_with(p, ext)) {
					files.push_back(di->path());
				} else {
					throw "The destination file exists and -o is not specified.";
				}
			}
		}

		if(files.size() == 0)
			throw "No matching files found.";

		if(!quiet)
			std::cerr << "Processing " << files.size() << " files." << std::endl;

		// Start
		las::WriterFactory wf;
		las::ReaderFactory rf;

		// The overall bounds: min x, max x, min y, max y, min z, max z
		double bounds[] = { DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG, DBL_MAX_POS, DBL_MAX_NEG };

		for(std::list<fs::path>::iterator it = files.begin(); it != files.end(); ++it) {

			const char * filename = it->c_str();

			std::cout << "Processing file " << filename << std::endl;

			// Open the source file.
			std::ifstream in(filename, std::ios::in | std::ios::binary);
			las::Reader r = rf.CreateWithStream(in);
			las::Header h = r.GetHeader();

			// Open the destination file.
			fs::path outfile(dst / it->leaf());
			las::Header dsth(h);

			int count = h.GetPointRecordsCount();

			MemRaster<double> x(count, 1);
			MemRaster<double> y(count, 1);
			MemRaster<double> z(count, 1);

			// Iterate over the points.
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				*(x.grid() + i) = pt.GetX();
				*(y.grid() + i) = pt.GetY();
				*(z.grid() + i) = pt.GetZ();
			}

			// Transform the points in-place.
			transformPoints(x, y, z, count, bounds);

			// Set bounds.
			dsth.SetMin(bounds[0], bounds[2], bounds[4]);
			dsth.SetMax(bounds[1], bounds[3], bounds[5]);
			std::ofstream out(outfile.c_str(), std::ios::out | std::ios::binary);
			las::Writer w(out, dsth);

			// Iterate over the points.
			r.Reset();
			for(int i = 0; r.ReadNextPoint(); ++i) {
				las::Point pt = r.GetPoint();
				pt.SetX(*(x.grid() + i));
				pt.SetY(*(y.grid() + i));
				pt.SetZ(*(z.grid() + i));
				// Write to the output
				w.WritePoint(pt);
			}

			// Close files.
			in.close();
			out.close();

		}

		return 0;

	}

	~Transformer() {
		pj_free(projFrom);
		pj_free(projECEF);
		pj_free(projTo);
		pj_free(projGeog);
	}

};

	
void usage() {
	std::cout << "Usage: las2csrs [options] <src file or dir> <dst dir> <src ref frame> <src epoch> <dst epoch> <src srid> <dst srid> [src geoid] [dst geoid]" << std::endl
		<< "This program converts coordinates from LAS files from any reference frame to NAD83(CSRS), " << std::endl
		<< "between any two epochs. If orthometric heights are used, provide the geoid file names (no path)," << std::endl
		<< "to ensure that the geoid models are used to convert to ellipsoidal heights for the calculations." << std::endl
		<< "Not doing this will result in small errors." << std::endl;

	std::cout << " -o     Overwrite existing files. Defaults to false." << std::endl;
	std::cout << " -v     Verbose output." << std::endl;
	std::cout << "Set LAS2CSRS_DATA to point to ITRF DB and grid shift file." << std::endl;
}

int test(int argc, char **argv) {

	if(argc < 10) {
		std::cerr << "Too few arguments." << std::endl;
		return 1;
	}
	try {

		int i = 1;
		std::string ffrom(argv[++i]);
		double efrom = atof(argv[++i]);
		double eto = atof(argv[++i]);
		int fsrid = atoi(argv[++i]);
		int tsrid = atoi(argv[++i]);
		double x = atof(argv[++i]);
		double y = atof(argv[++i]);
		double z = atof(argv[++i]);
		std::string fgeoid;
		std::string tgeoid;

		Transformer trans(ffrom, efrom, eto, fsrid, tsrid, fgeoid, tgeoid);

		MemRaster<double> xx(1,1);
		MemRaster<double> yy(1,1);
		MemRaster<double> zz(1,1);

		xx.set(0, x);
		yy.set(0, y);
		zz.set(0, z);

		double bounds[6];

		trans.transformPoints(xx, yy, zz, 1, bounds);
		std::cerr << xx[0] << " " << yy[0] << " " << zz[0] << std::endl;

	} catch(const char *e) {
		std::cerr << e << std::endl;
		return 1;
	}

	return 0;
}

int main(int argc, char **argv) {

	std::string comp("test");
	if(argc >= 2 && comp.compare(argv[1]) == 0) {
		return test(argc, argv);
	}

	if(argc < 8) {
		usage();
		return 1;
	}

	bool overwrite = false;

	try {

		int i = 1;
		for(; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg.c_str()[0] != '-') {
				break;
			} else if(arg == "-o") {
				overwrite = true;
			} else if(arg == "-v") {
				quiet = false;
			}
		}

		std::string srcfile(argv[i]);
		std::string dstdir(argv[++i]);
		std::string ffrom(argv[++i]);
		double efrom = atof(argv[++i]);
		double eto = atof(argv[++i]);
		int fsrid = atoi(argv[++i]);
		int tsrid = atoi(argv[++i]);
		std::string fgeoid;
		std::string tgeoid;
		std::cerr << argc << "," << i << std::endl;
		if(argc - 1 > i)
			fgeoid.assign(argv[++i]);
		std::cerr << argc << "," << i << std::endl;
		if(argc - 1 > i)
			tgeoid.assign(argv[++i]);

		std::cerr << std::setprecision(15);

		if(!std::getenv(LAS2CSRS_DATA))
			setenv(LAS2CSRS_DATA, "..", 1);

		Transformer trans(ffrom, efrom, eto, fsrid, tsrid, fgeoid, tgeoid);

		trans.transformLas(srcfile, dstdir, overwrite);

	} catch(const std::exception &err) {
		std::cerr << err.what() << std::endl;
		usage();
		return 1;
	} catch(const std::string &err) {
		std::cerr << err << std::endl;
		usage();
		return 1;
	} catch(const char *err) {
		std::cerr << err << std::endl;
		usage();
		return 1;
	} catch(...) {
		std::cerr << "An unknown exception occurred." << std::endl;
		usage();
		return 1;
	}

	return 0;

}

