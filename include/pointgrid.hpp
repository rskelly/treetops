/*
 * pointgrid.hpp
 *
 *  Created on: Dec 19, 2016
 *      Author: rob
 */

#ifndef INCLUDE_POINTGRID_HPP_
#define INCLUDE_POINTGRID_HPP_

#include <string>
#include <vector>

namespace geotools {

	namespace point {

		class PointgridConfig {
		public:
			double resolutionX;
			double resolutionY;
			double alignX;
			double alignY;
			int srid;
			std::vector<std::string> sourceFiles;
			std::string destFile;
			uint8_t method;
			std::set<uint8_t> classes;
		};

		class PointGrid {
		public:

			void process(const PointgridConfig &config,
					geotools::util::Callbacks *callbacks = nullptr, bool *cancel = nullptr);
};
	}
}




#endif /* INCLUDE_POINTGRID_HPP_ */
