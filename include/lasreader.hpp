#ifndef __LASREADER_HPP__
#define __LASREADER_HPP__

#include <cstdio>
#include <vector>
#include <memory>
#include <string>
#include <functional>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "util.hpp"
#include "raster.hpp"
#include "laspoint.hpp"

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

	namespace las {

		class LASFilter {
		private:
			double m_radius;			  // If zero, whole cell is used. Otherwise must be squared.
			uint8_t m_classes[255];       // Class ID is index into the array.

		public:
			LASFilter();
			void setClasses(const std::set<uint8_t> &classes);
			void setRadius(double radius);
			bool keep(const LASPoint &pt) const;
		};

		class LASReaderCallback {
		public:
			virtual void status(float) const = 0;
			virtual ~LASReaderCallback() = 0;
		};

		// Represents a single las file and enables reading.
		class LASReader {
		private:
			const static uint64_t BATCH_SIZE = 1000000;
			std::FILE *m_f;
			std::string m_file;
			uint16_t m_sourceId;
			std::string m_version;
			uint16_t m_headerSize;
			uint32_t m_offset;
			uint8_t m_pointFormat;
			uint16_t m_pointLength;
			uint64_t m_pointCount;
			uint64_t m_pointCountByReturn[5];
			double m_xScale;
			double m_yScale;
			double m_zScale;
			double m_xOffset;
			double m_yOffset;
			double m_zOffset;
			double m_xMin;
			double m_yMin;
			double m_zMin;
			double m_xMax;
			double m_yMax;
			double m_zMax;

			uint64_t m_curPoint;
			uint64_t m_batchSize;

			std::unique_ptr<Buffer> m_buf;
			std::queue<LASPoint> m_pts;

			const LASFilter *m_filter;

			// Prepares the reader for streaming points
			// by loading the header, etc.
			void init();

			// Load a batch of points.
			bool loadBatch();

		public:

			LASReader(const std::string &file);

			~LASReader();

			// Set the filter.
			void setFilter(const LASFilter *filter);

			// Reset the reader to enable reading from the start.
			void reset();

			// Populate the point with data for the next point in the file.
			// Return true if there is another point, false otherwise.
			bool next(LASPoint &pt);

			// Return the bounds of the file.
			Bounds bounds() const;

			// Return the number of points in the file.
			uint64_t pointCount() const;

		};

		// Enables reading from multiple LAS files as a single
		// stream of points.
		class LASMultiReader {
		private:
			bool *m_cancel;
			uint32_t m_idx;
			uint32_t m_cols;
			uint64_t m_cellCount;
			uint64_t m_pointCount;
			double m_resolutionX;
			double m_resolutionY;
			LASReader *m_reader;

			Bounds m_bounds;
			std::vector<std::string> m_files;
			std::vector<uint32_t> m_finalizer;
			std::vector<Bounds> m_blockBounds;

			const LASFilter *m_filter;

			void init(const std::vector<std::string> &files);

		public:
			LASMultiReader(const std::vector<std::string> &files, double resolutionX, double resolutionY, bool *cancel = nullptr);

			~LASMultiReader();

			// Set the filter.
			void setFilter(const LASFilter *filter);

			// Initializes the finalization grid.
			// Pass a functor to capture status updates from 0 to 1.
			void buildFinalizer(const LASReaderCallback *callback = nullptr);

			// Reset the reader to enable reading from the start.
			void reset();

			// Return the total number of cells used in the finalizer.
			uint64_t cellCount() const;

			// Return the total number of points.
			uint64_t pointCount() const;

			// Set the bounds of the las file set. This is used by the finalizer.
			// This is necessary if the output depends on a discretization of the
			// space with different origins than the point cloud itself.
			void setBounds(const Bounds &bounds);

			// Return the bounds of the point cloud.
			Bounds bounds() const;

			// Returns a vector containing the bounds of each file
			// comprising the point cloud.
			std::vector<Bounds> blockBounds() const;

			// Returns the next point.
			// Returns true if a point was successfully written to the output,
			// false otherwise.
			// If a pointer to final is given, it is set when a cell in the finalize
			// is completed.
			// If a pointer to finalIdx is given, it is set to the index that
			// was finalized.
			// If a pointer to fileChanged is given, it is set to true
			// when the reader switches to the next available LAS file.
			bool next(LASPoint &pt, bool *final = nullptr,
					uint64_t *finalIdx = nullptr, bool *fileChanged = nullptr);

		};

	}
}

#endif
