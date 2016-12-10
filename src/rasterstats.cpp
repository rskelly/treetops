/*
 * This progam computes some basic statistics for differences between
 * pairs of rasters, aggregated by class.
 *
 * The classes come from a classification raster, which is the first argument.
 * Then, each of the given files is differenced with each other file
 * to produce pairwise mean, variance and standard deviation within
 * the class.
 * 
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <map>

#include "geotools.hpp"
#include "util.hpp"
#include "raster.hpp"

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

	namespace raster {

		namespace util {

		/**
		 * Represents a single class, for a single pair of files.
		 */
		class Stat {
		private:
			// The list of values.
			std::vector<double> m_values;
			// Class.
			int m_class;
			// True if the list is sorted.
			bool m_sorted;

		public:

			/**
			 * Construct a Stat object.
			 */
			Stat() :
					m_class(0), m_sorted(false) {
			}

			/**
			 * Sorts the values.
			 */
			void sort() {
				if (!m_sorted) {
					std::sort(m_values.begin(), m_values.end());
					m_sorted = true;
				}
			}

			/**
			 * Add a point value to the list.
			 * The list ceases to be sorted.
			 */
			void add(double value) {
				m_values.push_back(value);
				m_sorted = false;
			}

			/**
			 * Return the number of values.
			 */
			int count() {
				return m_values.size();
			}

			/**
			 * Return the minimum value.
			 */
			double min() {
				if (count() > 0) {
					sort();
					return m_values[0];
				}
				return nan("");
			}

			/**
			 * Return the maximum value.
			 */
			double max() {
				if (count() > 0) {
					sort();
					return m_values[m_values.size() - 1];
				}
				return nan("");
			}

			/**
			 * Return the sum of values.
			 */
			double sum() {
				double sum = 0;
				for (size_t i = 0; i < m_values.size(); ++i)
					sum += m_values[i];
				return sum;
			}

			/**
			 * Return the mean of values.
			 */
			double mean() {
				if (count() > 0)
					return sum() / count();
				return nan("");
			}

			/**
			 * Return the median of values.
			 */
			double median() {
				int num;
				if ((num = count()) > 0) {
					sort();
					if (num % 2 == 0) {
						return (m_values[num / 2] + m_values[num / 2 - 1]) / 2.0;
					} else {
						return m_values[num / 2];
					}
				}
				return nan("");
			}

			/**
			 * Return the variance of values.
			 */
			double variance() {
				if (count() > 0) {
					double m = mean();
					double s = 0.0;
					for (size_t i = 0; i < m_values.size(); ++i)
						s += g_sq(m_values[i] - m);
					return s / (count() - 1.0);
				}
				return nan("");
			}

			/**
			 * Return the standard deviation of values.
			 */
			double stddev() {
				if (count() > 0)
					return sqrt(variance());
				return nan("");
			}

			~Stat() {
			}

		};

		} // Util

		class RasterStatsConfig {
		public:
			std::vector<std::string> sourceFiles;
			std::string destFile;
			std::string classFile;
			uint16_t kernelSize;
			uint8_t method;

			RasterStatsConfig() :
				kernelSize(3),
				method(0) {
			}

			void check() {
				if(method == 3 && kernelSize != 3) {
					g_warn("Kernel size for aspect is 3. Adjusting.");
					kernelSize = 3;
				}
			}
		};

		/**
		 * Compute statistics on the pairs of files in the files vector, for
		 * each class in the classification file. Print out a CSV table.
		 */
		void rasterclassstats(const RasterStatsConfig &config) {

			const std::string &clsfile = config.classFile;
			const std::vector<std::string> &files = config.sourceFiles;

			using namespace geotools::raster::util;

			if (clsfile.empty())
				throw std::invalid_argument("A classification raster is required.");

			if (files.size() < 2)
				throw std::invalid_argument("At least two data rasters are required.");

			std::map<std::string, std::map<std::string, std::map<int, Stat> > > stats;
			std::set<int> classes;

			Raster<unsigned char> clsrast(clsfile);

			// Loop over every pair of files.
			for (size_t f0 = 0; f0 < files.size(); ++f0) {
				for (size_t f1 = f0 + 1; f1 < files.size(); ++f1) {

					Raster<float> frast0(files[f0]);
					Raster<float> frast1(files[f1]);

					if (frast0.resolutionX() != frast1.resolutionX()
							|| frast0.resolutionY() != frast1.resolutionY())
						g_runerr("Rasters must have the same resolution.");

					int f0nodata = frast0.nodata();
					int f1nodata = frast1.nodata();

					double minx = g_max(frast0.minx(), frast1.minx());
					double miny = g_max(frast0.miny(), frast1.miny());
					double maxx = g_min(frast0.maxx(), frast1.maxx());
					double maxy = g_min(frast0.maxy(), frast1.maxy());

					for (double y = miny; y < maxy;
							y += std::abs(frast0.resolutionY())) {
						for (double x = minx; x < maxx;
								x += std::abs(frast0.resolutionX())) {

							int c0 = frast0.toCol(x);
							int r0 = frast0.toRow(y);
							int c1 = frast1.toCol(x);
							int r1 = frast1.toRow(y);
							int cc = clsrast.toCol(x);
							int cr = clsrast.toRow(y);

							// If any of the indices are out of bands, skip.
							if (c0 < 0 || c1 < 0 || cc < 0 || r0 < 0 || r1 < 0 || cr < 0
									|| c0 >= frast0.cols() || c1 >= frast1.cols()
									|| cc >= clsrast.cols() || r0 >= frast0.rows()
									|| r1 >= frast1.rows() || cr >= clsrast.rows())
								continue;

							float v0 = frast0.get(c0, r0);
							float v1 = frast1.get(c1, r1);

							// If either pixel is nodata, skip.
							if (v0 == f0nodata || v1 == f1nodata)
								continue;

							// Add the difference to the Stat object for the pair/class
							int cls = clsrast.get(cc, cr);
							classes.insert(cls);
							stats[files[f0]][files[f1]][cls].add(v0 - v1);
						}
					}
				}

			}

			// Print the header.
			std::cout << "file1,file2,class,count,sum,min,max,mean,variance,stddev"
					<< std::endl;

			// Print the data rows.
			for (auto it = stats.begin(); it != stats.end(); ++it) {
				for (auto it0 = it->second.begin(); it0 != it->second.end(); ++it0) {
					for (int cls : classes) {
						std::cout << it->first << "," << it0->first << "," << cls;
						if (it0->second.find(cls) == it0->second.end()) {
							std::cout << ",,,,,,,";
						} else {
							std::cout << "," << it0->second[cls].sum() << ","
									<< it0->second[cls].count() << ","
									<< it0->second[cls].min() << ","
									<< it0->second[cls].max() << ","
									<< it0->second[cls].mean() << ","
									<< it0->second[cls].variance() << ","
									<< it0->second[cls].stddev();
						}
						std::cout << std::endl;
					}
				}
			}
		}

		double median(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			int size = values.size();
			if(!size) return rast.nodata(band);
			std::sort(values.begin(), values.end());
			return size % 2 == 0 ? (values[size / 2 - 1] + values[size / 2]) / 2.0 : values[size / 2];
		}

		double mean(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			double sum = 0.0;
			double nodata = rast.nodata(band);
			uint16_t count = 0;
			for(const double &v : values) {
				if(v != nodata) {
					sum += v;
					++count;
				}
			}
			return values.size() ? sum / count : nodata;
		}

		double variance(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			double sum = 0.0;
			double nodata = rast.nodata(band);
			uint16_t count = 0;
			for(const double &v : values) {
				if(v != nodata) {
					sum += v;
					++count;
				}
			}
			if(count == 0) return nodata;
			sum = count > 0 ? sum / count : nodata;
			double sum2 = 0;
			for(const double &v : values) {
				if(v != nodata)
					sum2 += g_sq(v - sum);
			}
			return sum2;
		}

		double stddev(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			double var = variance(values, band, rast);
			double nodata = rast.nodata(band);
			return var == nodata ? nodata : std::sqrt(var);
		}

		// Shamelessly stolen from http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/how-aspect-works.htm#ESRI_SECTION1_4198691F8852475A9F4BC71246579FAA
		double aspect(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			if(values.size() != 9)
				g_argerr("Kernel must have 9 elements for aspect.");
			double a = ((values[2] + 2.0 * values[5] + values[8]) - (values[0] + 2.0 * values[3] + values[6])) / 8.0;
			double b = ((values[6] + 2.0 * values[7] + values[8]) - (values[0] + 2.0 * values[1] + values[3])) / 8.0;
			double aspect = 57.29578 * std::atan2(a, -b);
			if(aspect < 0.0) {
				return 90.0 - aspect;
			} else if(aspect > 90.0) {
			    return 360.0 - aspect + 90.0;
			} else {
			    return 90.0 - aspect;
			}
		}

		// Shamelessly stolen from http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/how-aspect-works.htm#ESRI_SECTION1_4198691F8852475A9F4BC71246579FAA
		double slope(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			if(values.size() != 9)
				g_argerr("Kernel must have 9 elements for slope.");
			double centre = values[4];
			double dif = 0.0;
			double nodata = rast.nodata(band);
			for(uint16_t i = 0; i < 9; ++i) {
				double v = values[i];
				if(v != nodata)
					dif = g_max(dif, g_abs(v - centre));
			}
			std::cerr << dif << ", " << std::sqrt(g_sq(rast.resolutionX()) + g_sq(rast.resolutionY())) << "\n";
			return dif / std::sqrt(g_sq(rast.resolutionX()) + g_sq(rast.resolutionY()));
		}

		void normalize(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			double sd = stddev(values, band, rast);
			double mn = mean(values, band, rast);
			for(int i = 0; i < values.size(); ++i)
				values[i] = sd == 0 ? rast.nodata(band) : (values[i] - mn) / sd;
		}

		double linearity(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			if(values.size() != 9)
				g_argerr("Kernel must have 9 elements for aspect.");
			normalize(values, band, rast);
			double maxVar = 0;
			int maxIdx = 0;
			double d[3];
			double x[4][3] = {{0, 4, 8}, {1, 4, 7}, {2, 4, 6}, {3, 4, 5}};
			for(int i = 0; i < 4; ++i) {
				double sum = 0;
				for(int j = 0; j < 3; ++j) {
					d[j] = values[x[i][j]];
					sum += values[x[i][j]];
				}
				double mean = sum / 3;
				double var = 0;
				for(int j = 0; j < 3; ++j)
					var += g_sq(values[x[i][j]] - mean);
				if(var > maxVar) {
					maxVar = var;
					maxIdx = i;
				}
			}
			return maxIdx;
		}

		double linearity2(std::vector<double> &values, uint16_t band, const Raster<float> &rast) {
			if(values.size() != 9)
				g_argerr("Kernel must have 9 elements for aspect.");
			normalize(values, band, rast);
			double maxVar = 0;
			int maxIdx = 0;
			double d[3];
			double x[4][3] = {{0, 4, 8}, {1, 4, 7}, {2, 4, 6}, {3, 4, 5}};
			double vars[4];
			for(int i = 0; i < 4; ++i) {
				double sum = 0;
				for(int j = 0; j < 3; ++j) {
					d[j] = values[x[i][j]];
					sum += values[x[i][j]];
				}
				double mean = sum / 3;
				double var = 0;
				for(int j = 0; j < 3; ++j)
					var += g_sq(values[x[i][j]] - mean);
				vars[i] = var;
			}
			double sum = 0;
			for(int i = 0; i < 3; ++i) 
				sum += vars[i];
			double mean = sum / 3;
			double var = 0;
			for(int i = 0; i < 3; ++i)
				var += g_sq(vars[i] - mean);
			return var;
		}

		typedef double (*statFunc)(std::vector<double>&, uint16_t, const Raster<float>&);

		void process(MemRaster<float> &inrast, MemRaster<float> &outrast,
				const Raster<float> &writerast,
				uint16_t kernelSize, uint16_t band, statFunc fn) {

			if(kernelSize % 2 == 0) {
				g_warn("Kernel size can't be even. Bumping up by one.");
				++kernelSize;
			}
			if(kernelSize < 3)
				g_argerr("Kernel size must be three or greater.");

			MemRaster<float> kernel(kernelSize, kernelSize);
			for(uint16_t row = 0; row < inrast.rows() - kernelSize; ++row) {
				for(uint16_t col = 0; col < inrast.cols() - kernelSize; ++col) {
					inrast.readBlock(col, row, kernel);
					std::vector<double> values;
					for(uint16_t r = 0; r < kernelSize; ++r) {
						for(uint16_t c = 0; c < kernelSize; ++c)
							values.push_back(kernel.get(c, r));
					}
					double value = fn(values, band, writerast);
					outrast.set(col + kernelSize / 2, row + kernelSize / 2, value);
				}
			}
		}

		statFunc getMethod(const RasterStatsConfig &config) {
			switch(config.method) {
			case 1:
				return mean;
			case 2:
				return median;
			case 3:
				return aspect;
			case 4:
				return slope;
			case 7:
				return variance;
			case 8:
				return stddev;
			case 9:
				return linearity;
			default:
				g_argerr("Unknown method: " << config.method);
			}
		}

		void rasterstats(RasterStatsConfig &config, Callbacks *callbacks = nullptr, bool *cancel = nullptr) {

			config.check();

			statFunc fn = getMethod(config);
			double resX = 0.0, resY = 0.0;
			uint16_t bands = 0;
			std::vector<double> nodata;
			Bounds bounds;
			bounds.collapse();
			for(const std::string &file : config.sourceFiles) {
				Raster<float> rast(file);
				bounds.extend(rast.bounds());
				if((resX != 0.0 && resX != rast.resolutionX()) || (resY != 0.0 && rast.resolutionY()))
					g_argerr("All rasters must have the same resolution.");
				if(bands > 0 && bands != rast.bandCount())
					g_argerr("All rasters must have the same number of bands.");
				bands = rast.bandCount();
				resX = rast.resolutionX();
				resY = rast.resolutionY();
				if(nodata.empty()) {
					for(uint16_t band = 1; band <= bands; ++band)
						nodata.push_back(rast.nodata(band));
				}
			}

			MemRaster<float> inrast(bounds.maxCol(resX) + 1, bounds.maxRow(resY) + 1, true);
			MemRaster<float> outrast(bounds.maxCol(resX) + 1, bounds.maxRow(resY) + 1, true);

			Raster<float> writerast(config.destFile, bands, bounds.minx(), 
				bounds.miny(), bounds.maxx(), bounds.maxy(), resX, resY, 0);

			for(uint32_t band = 1; band <= bands; ++band) {
				
				inrast.fill(nodata[band - 1]);
				writerast.setNodata(nodata[band - 1]);
				
				for(const std::string &file : config.sourceFiles) {
					Raster<float> rast(file);
					rast.readBlock(band, 0, 0, inrast, writerast.toCol(rast.leftx()), writerast.toRow(rast.topy()));
				}

				process(inrast, outrast, writerast, config.kernelSize, band, fn);

				writerast.writeBlock(band, outrast);
			}

		}

		void contributingarea(RasterStatsConfig &config, Callbacks *callbacks = nullptr, bool *cancel = nullptr) {

			config.check();

			Raster<float> source(config.sourceFiles[0]);

			MemRaster<float> inrast(source.cols(), source.rows(), false);
			inrast.writeBlock(source);

			MemRaster<uint32_t> outrast(source.cols(), source.rows(), false);
			outrast.fill(0);

			std::queue<std::pair<uint16_t, uint16_t> > q;

			double nodata = source.nodata();
			uint64_t total = inrast.size();
			uint64_t count = 0;
			uint32_t cc = 0;

			std::vector<bool> visited(inrast.cols() * inrast.rows());

			for(uint16_t row = 0; row < inrast.rows(); ++row) {
				for(uint16_t col = 0; col < inrast.cols(); ++col) {
					if(++count % 1000 == 0)
						std::cerr << " - " << (int) ((float) count / total * 100) << "\n";
					std::fill(visited.begin(), visited.end(), false);
					q.push(std::make_pair(col, row));
					while(!q.empty()) {
						auto p = q.front();
						q.pop();
						uint16_t c = p.first;
						uint16_t r = p.second;
						if(visited[r * inrast.cols() + c]) continue;
						visited[r * inrast.cols() + c] = true;
						double v = inrast.get(c, r);
						if(v == nodata) continue;
						if(c > 0 && inrast.get(c - 1, r) < v) {
							outrast.set(c - 1, r, outrast.get(c - 1, r) + 1);
							q.push(std::make_pair(c - 1, r));
						}
						if(c <= inrast.cols() - 1 && inrast.get(c + 1, r) < v) {
							outrast.set(c + 1, r, outrast.get(c + 1, r) + 1);
							q.push(std::make_pair(c + 1, r));
						}
						if(r > 0 && inrast.get(c, r - 1) < v) {
							outrast.set(c, r - 1, outrast.get(c, r - 1) + 1);
							q.push(std::make_pair(c, r - 1));
						}
						if(r <= inrast.rows() - 1 && inrast.get(c, r + 1) < v) {
							outrast.set(c, r + 1, outrast.get(c, r + 1) + 1);
							q.push(std::make_pair(c, r + 1));
						}
					}
				}
			}

			Raster<uint32_t> dest(config.destFile, 1, source);
			dest.writeBlock(outrast);
		}

	} // raster

} // geotools

void usage() {
	std::cerr
			<< "Usage: rasterstats [options] <raster [raster [...]]>\n"
			<< "	Produces a table containing statistics for each class in the class raster\n"
			<< "	for the difference in every pixel for each pair of rasters.\n"
			<< "    Or, produces a raster containing statistics derived from a window\n"
			<< "    on an input raster.\n"
			<< " -m <method>     An action to perform. Currently available:\n"
			<< "                 mean, median, slope, aspect, contributing area (ca), class.\n"
			<< " -c <filename>   The classification raster used to generate classes (use with class method.)\n"
			<< "                 Implies class statistics option.\n"
			<< " -o <filename>   Output raster. Input rasters will be merged into output.\n"
			<< " -k <int>        Kernel size.\n";
}

int main(int argc, char ** argv) {

	if(argc == 1) {
		usage();
		return 1;
	}

	g_loglevel(G_LOG_DEBUG);

	using namespace geotools::raster;

	RasterStatsConfig config;

	std::map<std::string, uint8_t> statTypes;
	statTypes["mean"] = 1;
	statTypes["median"] = 2;
	statTypes["aspect"] = 3;
	statTypes["slope"] = 4;
	statTypes["ca"] = 5;
	statTypes["class"] = 6;
	statTypes["variance"] = 7;
	statTypes["stddev"] = 8;
	statTypes["linearity"] = 9;

	try {

		for (int i = 1; i < argc; ++i) {
			std::string arg(argv[i]);
			if(arg == "-c") {
				config.classFile = argv[++i];
			} else if(arg == "-o") {
				config.destFile = argv[++i];
			} else if(arg == "-k") {
				config.kernelSize = atoi(argv[++i]);
			} else if(arg == "-m") {
				std::string t = argv[++i];
				if(statTypes.find(t) == statTypes.end())
					g_argerr("Unknown method: " << t);
				config.method = statTypes[t];
			} else {
				config.sourceFiles.push_back(argv[i]);
			}
		}

		switch(config.method) {
		case 6:
			rasterclassstats(config);
			break;
		case 5:
			contributingarea(config);
			break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 7:
		case 8:
		case 9:
			rasterstats(config);
			break;
		default:
			g_argerr("Unknown method: " << config.method);
		}

	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		usage();
		return 1;
	}

	return 0;
}
