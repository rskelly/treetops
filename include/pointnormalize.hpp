#ifndef __POINTNORMALIZE_HPP__
#define __POINTNORMALIZE_HPP__

#include <string>
#include <list>

#include "util.hpp"

namespace geotools {

    namespace point {

        class PointNormalizeConfig {
        public:
            std::string outputDir;
            std::vector<std::string> sourceFiles;
            bool dropNegative;
            bool dropGround;
            uint32_t threads;
            bool overwrite;
            
            PointNormalizeConfig() :
                dropNegative(true),
                dropGround(true),
                threads(1),
                overwrite(true) {
                
            }
        };

        class PointNormalize {
        public:
            void normalize(const PointNormalizeConfig &config, 
                    const geotools::util::Callbacks *callbacks = nullptr, 
                    bool *cancel = nullptr);
        };

    } // point

} // geotools

#endif
