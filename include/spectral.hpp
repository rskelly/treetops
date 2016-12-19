#ifndef __SPECTRAL_HPP__
#define __SPECTRAL_HPP__

#include <vector>
#include <string>
#include <set>

#include "geotools.hpp"

namespace geotools {

    namespace spectral {

        namespace config {

            class SpectralConfig {
            public:
            	// Filenames of hyperspectral files.
                std::vector<std::string> spectralFilenames;
                // Filename of the raster that contains ROIs with IDs.
                std::string indexFilename;
                // Filename of output database.
                std::string outputFilename;
                // A list of bands to process.
                std::set<int> bands;
                // Nodata value for spectral files.
                unsigned short specNodata;
                // Nodata value for index file.
                unsigned int idxNodata;
                // True if a spectral nodata has been configured.
                bool hasSpecNodata;
                // True if an index nodata has been configured.
                bool hasIdxNodata;
                // The SRID of the output.
                int srid;

                SpectralConfig() :
					specNodata(0),
					idxNodata(0),
					hasSpecNodata(false),
					hasIdxNodata(false),
					srid(0) {
                }

                void check() const {
                    if (spectralFilenames.size() == 0)
                        g_argerr("There must be at least one spectral file.");
                    if (indexFilename.empty())
                        g_argerr("The index filename must be given.");
                    if (outputFilename.empty())
                        g_argerr("The output filename must be given.");
                    if (srid == 0)
                        g_warn("The SRID is zero.");
                    if(hasIdxNodata)
                    	g_warn("The index file null value is not set.");
                    if(hasSpecNodata)
                    	g_warn("The spectral file null value is not set.");
                }

            };


        } // config

        class Spectral {
        private:
        	geotools::util::Callbacks *m_callbacks;
        	bool *m_cancel;

        public:
        	Spectral();
            void extractSpectra(const geotools::spectral::config::SpectralConfig &config,
            		geotools::util::Callbacks *callbacks = nullptr, bool *cancel = nullptr);
        };

    } // spectra

} // geotools


#endif                                                                                                                               
