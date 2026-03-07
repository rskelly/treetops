#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>
#include <vector>
#include <memory>
#include <cstdio>
#include <mutex>

#include <geos_c.h>

#include "settings.hpp"
#include "util.hpp"
#include "db.hpp"
#include "grid.hpp"
//#include "ds/mqtree.hpp"
//#include "ds/simple_interval_tree.hpp"

namespace tt {

	namespace data {

		class Treetop {
		public:
			int id;
			int col;
			int row;
			int window;
			float height;

			float ox;
			float oy;
			float oz;
			float groundZ;
			float sx;
			float sy;
			float sz;
			int parentId;

			Treetop() : Treetop(0, 0, 0, 0, 0) {}

			Treetop(int id, int col, int row, int window, float height) :
				id(id),
				col(col),
				row(row),
				window(window),
				height(height),
				ox(0), oy(0), oz(0),
				sx(0), sy(0), sz(0) {}

		};

	} // data

	/**
	 * The main class of the Treetops application. Contains the high-level
	 * API for performing the process.
	 */
	class Treetops {
	private:
		tt::config::Settings* m_settings;

	public:

		/**
		 * Construct the Treetops object.
		 */
		Treetops(tt::config::Settings& settings);

		/**
		 * A convenience method for smoothing the input raster before using it to
		 * generate crowns or treetops.
		 */
		void smooth();

		/**
		 * Locates tree top points on a canopy height model and saves them to a database.
		 */
		void treetops();

		/**
		 * Performs tree crown delineation using a (preferrably smoothed) input raster and a
		 * vector file (sqlite) containing tree tops as seeds. Output is an integer raster with
		 * cell values representing tree top IDs, and an optional vector which is the polygonized
		 * version of the raster. The table should have been generated using the treetops() method
		 * to ensure that its structure is correct.
		 */
		void treecrowns();

		~Treetops();
	};

	// 	/**
	// 	 * Represents a single pixel.
	// 	 */
	// 	class Px {
	// 	public:
	// 		int col;	///<! The column.
	// 		int row;	///<! The row.

	// 		/**
	// 		 * Construct a pixel using the given column and row.
	// 		 *
	// 		 * \param col The column.
	// 		 * \param row The row.
	// 		 */
	// 		Px(int col, int row) : col(col), row(row) {}

	// 		/**
	// 		 * Construct a default pixel at (0, 0).
	// 		 */
	// 		Px() : Px(0, 0) {}
	// 	};

	// 	/**
	// 	 * Represents a tree top.
	// 	 */
	// 	class Treetop {
	// 	public:
	// 		int col;	///<! The grid column of the top.
	// 		int row;	///<! The grid row of the top.
	// 		int window;	///<! The size of the window used to locate the top.

	// 		/**
	// 		 * Construct a Treetop with the given row, column and window size.
	// 		 *
	// 		 * \param col The grid column of the top.
	// 		 * \param row The grid row of the top.
	// 		 * \param window The size of the window used to locate the top.
	// 		 */
	// 		Treetop(int col, int row, int window) : col(col), row(row), window(window) {}

	// 		/**
	// 		 * Construct a default Treetop with default values (0, 0, 0).
	// 		 */
	// 		Treetop() : Treetop(0, 0, 0) {}
	// 	};

	// 	/**
	// 	 * Represents a delineated tree crown.
	// 	 */
	// 	class Crown {
	// 	public:
	// 		double x;	///<! The projected x coordinate of the crown.
	// 		double y;	///<! The projected y coordinate of the crown.
	// 		double z;	///<! The height of the crown.

	// 		/**
	// 		 * Construct a Crown with the given x, y and z coordinates.
	// 		 *
	// 		 * \param x The projected x coordinate of the crown.
	// 		 * \param y The projected y coordinate of the crown.
	// 		 * \param z The height of the crown.
	// 		 */
	// 		Crown(double x, double y, double z) :
	// 			x(x), y(y), z(z) {}

	// 		/**
	// 		 * Construct a default Crown with default values (0, 0, 0).
	// 		 */
	// 		Crown() : Crown(0, 0, std::numeric_limits<double>::lowest()) {}
	// 	};

	// 	class TreetopsConfig;

	// 	/**
	// 	 * Enum of mutable fields in the TreetopsConfig object.
	// 	 */
	// 	enum TreetopsConfigField : uint64_t {
	// 		BoundsField = 1L << 0,
	// 		BuildIndex = 1L << 1,
	// 		TableCacheSize = 1L << 2,
	// 		RowCacheSize = 1L << 3,
	// 		Threads = 1L << 4,
	// 		DoSmoothing = 1L << 5,
	// 		SmoothWindowSize = 1L << 6,
	// 		SmoothSigma = 1L << 7,
	// 		OriginalCHM = 1L << 8,
	// 		OriginalCHMBand = 1L << 9,
	// 		SmoothedCHM = 1L << 10,
	// 		SmoothedCHMDriver = 1L << 11,
	// 		DoTops = 1L << 12,
	// 		TopsThresholds = 1L << 13,
	// 		TreetopsDatabase = 1L << 14,
	// 		TreetopsDatabaseDriver = 1L << 15,
	// 		TopsMaxNulls = 1L << 16,
	// 		DoCrowns = 1L << 17,
	// 		CrownsThresholds = 1L << 18,
	// 		CrownsUpdateHeights = 1L << 19,
	// 		CrownsRaster = 1L << 20,
	// 		CrownsRasterDriver = 1L << 21,
	// 		CrownsDoDatabase = 1L << 22,
	// 		CrownsDatabase = 1L << 23,
	// 		CrownsDatabaseDriver = 1L << 24,
	// 		CrownsRemoveHoles = 1L << 25,
	// 		CrownsRemoveDangles = 1L << 26,
	// 		SettingsFile = 1L << 27,
	// 		CrownsKeepSmoothed = 1L << 28,
	// 		TopsDBFormatChanged = 1L << 29,		// Used when the program forces a change in format e.g., to accommodate filesize.
	// 		CrownsDBFormatChanged = 1L << 30	// Used when the program forces a change in format e.g., to accommodate filesize.
	// 	};

	// 	/**
	// 	 * \brief An exception to identify failures in saving tops or crowns.
	// 	 */
	// 	class CrownDBSaveException : public std::runtime_error {
	// 	public:
	// 		CrownDBSaveException(const std::string& msg) :
	// 			std::runtime_error(msg) {}
	// 	};

	// 	/**
	// 	 * A fast, sequential database for storing Crown and Treetop instances.
	// 	 *
	// 	 * Each top record is constant length, stored as binary. The poly_position field gives
	// 	 * the location of the polygon in the crowns file. If poly_len is zero, there is no crown,
	// 	 * if it is >0, the crown is read from the file at poly_position.
	// 	 *
	// 	 * id (size_t) | top (Top) | poly_len (size_t) | poly_position (size_t)
	// 	 *
	// 	 * The crowns file contains crown polygons. This file has no structure. Crowns
	// 	 * must be located and read using the poly_len and poly_position values in the
	// 	 * tops file. The crowns are WKB.
	// 	 */
	// 	class CrownDB {
	// 	private:
	// 		int m_tfd; 								///<! File descriptor of open tops database.
	// 		int m_cfd; 								///<! File descriptor of open crowns database.
	// 		size_t m_toffset;						///<! Current offset in the tops file.
	// 		size_t m_coffset;						///<! Current offset in the crowns file.
	// 		size_t m_tidx;							///<! The iterator index.
	// 		int m_count;
	// 		GEOSContextHandle_t m_gctx;				///<! GEOS context.
	// 		GEOSWKBReader* m_rdr;					///<! WKB Reader.
	// 		GEOSWKBWriter* m_wtr;					///<! WKB Writer.
	// 		std::string m_tfile;					///<! The tops database filename.
	// 		std::string m_cfile;					///<! The crowns database filename.
	// 		std::unordered_map<int, size_t> m_idx; 	///<! Mapping from top ID to offset in the file.

	// 		std::mutex m_imtx;
	// 		std::mutex m_tmtx;
	// 		std::mutex m_cmtx;

	// 		/**
	// 		 * \brief Write the crown to the crowns file and return its position.
	// 		 *
	// 		 * If the length of the geometry is greater than the original
	// 		 * geometry (if it exists) it'll be written at the end and the original
	// 		 * will be orphaned.
	// 		 *
	// 		 * The new length and position are set into the plen and ppos values.
	// 		 *
	// 		 * True is returned if the process succeeds.
	// 		 */
	// 		bool writeCrown(size_t& plen, size_t& ppos, const GEOSGeometry* crown, GEOSContextHandle_t gctx);

	// 	public:

	// 		/**
	// 		 * \brief Create the file DB using temporary files.
	// 		 */
	// 		CrownDB();

	// 		/**
	// 		 * \brief Return the internal GEOS context handle.
	// 		 */
	// 		GEOSContextHandle_t gctx();

	// 		/**
	// 		 * \brief Attempt to load the Top and crown for the given ID.
	// 		 *
	// 		 * \param id The top ID.
	// 		 * \param[out] top A reference to a writable Top.
	// 		 * \param[out] crown A reference to a writable pointer to the crown geometry.
	// 		 * \param gctx A GEOSContextHandle_t.
	// 		 * \param withCrown If true, loads the crown geometry. Otherwise not.
	// 		 * \return True if the record was found.
	// 		 */
	// 		bool find(int id, tt::util::Top& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx = nullptr, bool withCrown = true);

	// 		/**
	// 		 * \brief Read the next available Top and crown geometry.
	// 		 *
	// 		 * \param[inout] top A Top.
	// 		 * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
	// 		 * \param gctx An optional GEOS context.
	// 		 * \param loadCrown If false, the crown is not loaded.
	// 		 * \return True if a record was found.
	// 		 */
	// 		bool next(tt::util::Top& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx = nullptr, bool loadCrown = false);

	// 		/**
	// 		 * \brief Read the next available Top and crown geometry.
	// 		 *
	// 		 * \param[inout] top A Top.
	// 		 * \param[inout] crownwkb A char buffer representing the crown's WKB representation.
	// 		 * \param loadCrown If false, the crown is not loaded.
	// 		 * \return True if a record was found.
	// 		 */
	// 		bool next(tt::util::Top& top, std::vector<unsigned char>& crownwkb, bool loadCrown = false);

	// 		bool next(tt::util::Top& top);

	// 		int count();

	// 		/**
	// 		 * \brief Reset the iterator.
	// 		 */
	// 		void reset();

	// 		/**
	// 		 * \brief Updates or adds the Treetop and Crown.
	// 		 *
	// 		 * If the Top isn't already in the DB it will be inserted at the end.
	// 		 *
	// 		 * \param[inout] top A Treetop.
	// 		 * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
	// 		 * \return True on success.
	// 		 */
	// 		bool insert(const tt::util::Top& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx = nullptr);

	// 		/**
	// 		 * \brief Update the Treetop and Crown.
	// 		 *
	// 		 * \param[inout] top A Treetop.
	// 		 * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
	// 		 * return True on success.
	// 		 */
	// 		bool update(const tt::util::Top& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx = nullptr);

	// 		/**
	// 		 * \brief Save the Tops into the given database file and layer.
	// 		 */
	// 		void saveTops(const std::string& file, const std::string& driver,
	// 				const std::string& layer, const std::string& projection,
	// 				Monitor* monitor);

	// 		/**
	// 		 * \brief Save the crowns with geometry into the given database file and layer.
	// 		 */
	// 		void saveCrowns(const std::string& file, const std::string& driver,
	// 				const std::string& layer, const std::string& projection,
	// 				Monitor* monitor);

	// 		/**
	// 		 * \brief Return the filesize of the crowns file in bytes.
	// 		 */
	// 		size_t crownsDBSize();

	// 		/**
	// 		 * \brief Return the filesize of the tops file in bytes.
	// 		 */
	// 		size_t topsDBSize();

	// 		/**
	// 		 * \brief Return true if the tops database can be saved using the given driver.
	// 		 */
	// 		bool checkSaveTops(const std::string& driver);

	// 		/**
	// 		 * \brief Return true if the crowns database can be saved using the given driver.
	// 		 */
	// 		bool checkSaveCrowns(const std::string& driver);

	// 		~CrownDB();
	// 	};

	// 	/**
	// 	 * Gives an implementor the ability to receive updates from a TreetopsConfig object.
	// 	 */
	// 	class TreetopsConfigListener {
	// 	public:

	// 		/**
	// 		 * Called when the TreetopsConfig object generates an update.
	// 		 *
	// 		 * \param config The originating TreetopsConfig object.
	// 		 * \param field An value containing ORed field IDs given by the TreetopsConfigField enum.
	// 		 */
	// 		virtual void configUpdate(TreetopsConfig& config, long field) = 0;

	// 		virtual ~TreetopsConfigListener() {}
	// 	};

	// 	/**
	// 	 * Contains configuration information for performing tree top extraction.
	// 	 */
	// 	class TreetopsConfig {
	// 	private:
	// 		tt::util::Bounds<float> m_bounds;				///<! Defines the boundaries of work to be performed. Every step of the process will be confined, including smoothing and searching.
	// 		int m_srid;										///<! A spatial reference ID for the output files.
	// 		bool m_buildIndex;								///<! If true, build the index on the tops table. Can be slow.
	// 		int m_tableCacheSize;							///<! The cache size (B) for the sqlite database. A performance optimization.
	// 		int m_rowCacheSize;								///<! The cache size (B) for rows when reading the raster.
	// 		int m_threads;									///<! The number of threads to use in execution.
	// 		std::string m_settings;							///<! File for storing settings.

	// 		bool m_doSmoothing;								///<! Set to true to perform smoothing on the input raster. This will force a check that the smoothing params are valild.
	// 		bool m_doTops;									///<! If true, treetop location will be performed.
	// 		bool m_doCrowns;								///<! Set to true to delineate crowns.

	// 		std::string m_originalCHM;						///<! The path to the original CHM.
	// 		int m_originalCHMBand;							///<! The raster band to smooth.

	// 		std::string m_smoothedCHM;						///<! The path to the smoothed CHM.
	// 		std::string m_smoothedCHMDriver;				///<! The driver to use for creating the smoothed CHM.
	// 		int m_smoothWindowSize;							///<! The size of the smoothing window >=3; an odd number.
	// 		double m_smoothSigma;							///<! The std. deviation used for generating the Gaussian kernel. 0 < n <= 1.

	// 		std::string m_treetopsDatabase;					///<! The path to the treetops database.
	// 		std::string m_treetopsDatabaseDriver;			///<! The river to use for creating the database.
	// 		std::vector<TopThreshold> m_topsThresholds;		///<! For pixels equal or above each height (double) use the given
	// 														///<! window size to detect maxima. Previously-detected maxima will be
	// 														///<! obliterated if a new maximum is found whose window encompases
	// 														///<! the previous one.
	// 		double m_topsMaxNulls;							///<! The max proportion of pixels in a given kernel that are allowed
	// 														///<! to be null. Any kernel with this many or greater is ignored.

	// 		std::string m_crownsRaster;						///<! The path to the crowns raster.
	// 		std::string m_crownsRasterDriver;				///<! The driver to use for the raster.
	// 		std::string m_crownsDatabase;					///<! The path to the crowns database.
	// 		std::string m_crownsDatabaseDriver;				///<! The driver to use for the database.
	// 		std::vector<CrownThreshold> m_crownsThresholds;	///<! The crown delineation thresholds: min height, height percentage and crown radius.
	// 														///<! If a pixel is above the minimum height, is within a given percentage
	// 														///<! of the top height and is within a given radius of the top, it may
	// 														///<! be included in a crown.
	// 		bool m_crownsUpdateHeights;						///<! If true, the heights of the treetops will be updated
	// 														///<! using values from the original CHM within the bounds
	// 														///<! of the crowns
	// 		bool m_crownsDoDatabase;						///<! If true, a crowns database will be produced
	// 		bool m_crownsRemoveHoles;						///<! True to remove holes in vector polygons.
	// 		bool m_crownsRemoveDangles;						///<! True to remove dangles (diagonally-connected sub-polys) in vector polygons.
	// 		bool m_crownsKeepSmoothed;						///<! If true, keep the smoothed treetop heights.

	// 		TreetopsConfigListener* m_listener;				///<! Listens for updates from the config object.
	// 		Monitor* m_monitor;
	// 		bool m_active;									///<! True if the object is dispatching events. Otherwise silent.

	// 		bool m_locked;									///<! True if the object will not accept changes.

	// 		std::unique_ptr<tt::grid::Grid<float>> m_chm;
	// 		std::unique_ptr<tt::grid::Grid<float>> m_smoothed;
	// 		bool m_smoothExisted;
	// 		std::unique_ptr<tt::grid::Grid<uint32_t>> m_topsWindow;
	// 		std::unique_ptr<tt::grid::Grid<uint32_t>> m_topsID;
	// 		std::unique_ptr<tt::grid::Grid<uint32_t>> m_topsParent;
	// 		std::unique_ptr<tt::grid::Grid<uint32_t>> m_crowns;

	// 		std::unique_ptr<tt::ds::mqtree<tt::util::Top>> m_topsTree; 	///<! Used to maintain a mapped tops database if both treetops
	// 																				///<! and treecrowns are expected to be invoked in sequence.

	// 		std::unique_ptr<CrownDB> m_crownDb;										///<! Database file for storing crowns and tops; to convert to final DB format.

	// 	public:

	// 		tt::grid::Grid<float>& chm();
	// 		tt::grid::Grid<float>& smoothed();
	// 		bool smoothExisted() const;
	// 		tt::grid::Grid<uint32_t>& crowns();
	// 		tt::grid::Grid<uint32_t>& topsWindow();
	// 		tt::grid::Grid<uint32_t>& topsID();
	// 		tt::grid::Grid<uint32_t>& topsParent();
	// 		tt::ds::mqtree<tt::util::Top>& topsTree();

	// 		CrownDB& crownDB();

	// 		/**
	// 		 * \brief Reset anything that needs to be before starting a run.
	// 		 */
	// 		void reset();

	// 		void flush();

	// 		void destroy();

	// 		void setMonitor(Monitor* monitor);

	// 		Monitor* monitor();

	// 		// Build a TreetopsConfig with defaults.
	// 		TreetopsConfig();

	// 		/**
	// 		 * Locks the object to prevent changes.
	// 		 */
	// 		void lock();

	// 		/**
	// 		 * Unlocks the object to allow changes.
	// 		 */
	// 		void unlock();

	// 		/**
	// 		 * If there is a listener update it with the given field(s).
	// 		 *
	// 		 * \param field The updated field, or an integer containing the ored fields.
	// 		 */
	// 		void update(long field);

	// 		/**
	// 		 * Set to true to let the object dispatch updates.
	// 		 *
	// 		 * \param active True if the object is distpatching updates.
	// 		 * \return The original state.
	// 		 */
	// 		bool setActive(bool active);

	// 		/**
	// 		 * Returns true if the object is distpatching updates.
	// 		 *
	// 		 * \return True if the object is distpatching updates.
	// 		 */
	// 		bool active() const;

	// 		/**
	// 		 * Set the settings file path.
	// 		 *
	// 		 * \param filename The settings file path.
	// 		 */
	// 		void setSettings(const std::string& filename);

	// 		/**
	// 		 * Get the settings file path.
	// 		 *
	// 		 * \return The settings file path.
	// 		 */
	// 		const std::string& settings() const;

	// 		// Check that the settings are appropriate for a smoothing
	// 		// job, throw an exception otherwise.
	// 		void checkSmoothing() const;

	// 		// Check that the settings are appropriate for a treetops
	// 		// job, throw an exception otherwise.
	// 		void checkTops() const;

	// 		// Check that the settings are appropriate for a crowns
	// 		// job, throw an exception otherwise.
	// 		void checkCrowns() const;

	// 		// Check that the settings are appropriate for a merge
	// 		// job, throw an exception otherwise.
	// 		void checkMerge() const;

	// 		// Check the validity of the configuration.
	// 		void check() const;

	// 		// Returns true if any function can be successfully run.
	// 		bool canRun() const;

	// 		// Returns the tops thresholds as a comma-delimited list.
	// 		std::string topsThresholdsList() const;

	// 		// Parses a comma-delimited list of tops thresholds into an internal list.
	// 		void parseTopsThresholds(const std::string& str);

	// 		// Returns the crowns thresholds as a comma-delimited list.
	// 		std::string crownsThresholdsList() const;

	// 		// Parses a comma-delimited list of crowns thresholds into an internal list.
	// 		void parseCrownsThresholds(const std::string& str);

	// 		void setListener(TreetopsConfigListener* listener);
	// 		TreetopsConfigListener* listener();

	// 		void setBounds(const Bounds<float>& bounds);
	// 		const Bounds<float>& bounds() const;

	// 		// Attempts to discover the projection from one of the configured rasters.
	// 		// In order: the unsmoothed raster, smoothed, crowns. If none is found,
	// 		// returns zero.
	// 		std::string projection();

	// 		void setBuildIndex(bool build);
	// 		bool buildIndex() const;

	// 		void setTableCacheSize(int size);
	// 		int tableCacheSize() const;

	// 		void setRowCacheSize(int size);
	// 		int rowCacheSize() const;

	// 		void setThreads(int threads);
	// 		int threads() const;

	// 		void setOriginalCHM(const std::string& filename, bool formatAll = false);
	// 		const std::string& originalCHM() const;

	// 		void setOriginalCHMBand(int band);
	// 		int originalCHMBand() const;

	// 		void setSmoothedCHM(const std::string& filename);
	// 		const std::string& smoothedCHM() const;

	// 		void setSmoothedCHMDriver(const std::string& driver);
	// 		const std::string& smoothedCHMDriver() const;

	// 		void setTreetopsDatabase(const std::string& filename);
	// 		const std::string& treetopsDatabase() const;

	// 		/**
	// 		 * \brief Set the treetops database driver.
	// 		 *
	// 		 * \param driver The driver.
	// 		 * \param forced True if the program has forced the driver to change due to file size or errors.
	// 		 */
	// 		void setTreetopsDatabaseDriver(const std::string& driver, bool forced = false);
	// 		const std::string& treetopsDatabaseDriver() const;

	// 		void setCrownsRaster(const std::string& filename);
	// 		const std::string& crownsRaster() const;

	// 		void setCrownsRasterDriver(const std::string& filename);
	// 		const std::string& crownsRasterDriver() const;

	// 		void setCrownsDatabase(const std::string& filename);
	// 		const std::string& crownsDatabase() const;

	// 		/**
	// 		 * \brief Set the crowns database driver.
	// 		 *
	// 		 * \param driver The driver.
	// 		 * \param forced True if the program has forced the driver to change due to file size or errors.
	// 		 */
	// 		void setCrownsDatabaseDriver(const std::string& driver, bool forced = false);
	// 		const std::string& crownsDatabaseDriver() const;

	// 		void setCrownsKeepSmoothed(bool keep);
	// 		bool crownsKeepSmoothed() const;

	// 		void setDoSmoothing(bool smoothing);
	// 		bool doSmoothing() const;

	// 		void setSmoothWindowSize(int size);
	// 		int smoothWindowSize() const;

	// 		void setSmoothSigma(double sigma);
	// 		double smoothSigma() const;

	// 		void setDoTops(bool tops);
	// 		bool doTops() const;

	// 		void setTopsThresholds(const std::vector<TopThreshold>& thresholds);
	// 		const std::vector<TopThreshold>& topsThresholds() const;

	// 		void setTopsMaxNulls(double nulls);
	// 		double topsMaxNulls() const;

	// 		void setDoCrowns(bool crowns);
	// 		bool doCrowns() const;

	// 		void setCrownsThresholds(const std::vector<CrownThreshold>& thresholds);
	// 		const std::vector<CrownThreshold> crownsThresholds() const;

	// 		void setCrownsUpdateHeights(bool updateHeights);
	// 		bool crownsUpdateHeights() const;

	// 		void setCrownsDoDatabase(bool doDatabase);
	// 		bool crownsDoDatabase() const;

	// 		void setCrownsRemoveHoles(bool removeHoles);
	// 		bool crownsRemoveHoles() const;

	// 		void setCrownsRemoveDangles(bool removeDangles);
	// 		bool crownsRemoveDangles() const;

	// 		~TreetopsConfig();

	// 	};

	// } // config

	// namespace util {

	// 	using namespace tt::grid;
	// 	using namespace tt::ds;
	// 	using namespace tt::config;

	// 	/**
	// 	 * A simple class for maintaining information about a tree top.
	// 	 */
	// 	class Top {
	// 	public:
	// 		size_t pos;			///<! The position in the qtree's memory.
	// 		size_t id;			///<! The ID of this top.
	// 		size_t parentID;	///<! The ID of this top's parent.
	// 		double ox, oy, oz; 	///<! Original x, y, z value
	// 		double sx, sy, sz; 	///<! Smoothed x, y, z value.
	// 		double groundZ;		///<! The height of the ground extracted from a DTM.
	// 		int sc, sr;    		///<! Smoothed col, row.

	// 		/**
	// 		 * Construct a top.
	// 		 */
	// 		Top();

	// 		/**
	// 		 * Construct a top with the given ID, parent and properties.
	// 		 *
	// 		 * \param id The unique ID of this top.
	// 		 * \param parentId The ID of this top's parent. The parent is
	// 		 *                 a top which is higher and whose window encompasses
	// 		 *                 the current one.
	// 		 * \param ox The original x coordinate.
	// 		 * \param oy The original y coordinate.
	// 		 * \param oz The original z coordinate.
	// 		 * \param sx The smoothed x coordinate.
	// 		 * \param sy The smoothed y coordinate.
	// 		 * \param sz The smoothed z coordinate.
	// 		 * \param sc The smoothed column.
	// 		 * \param sr The smoothed row.
	// 		 */
	// 		Top(size_t id, size_t parentId,
	// 						double ox, double oy, double oz,
	// 						double sx, double sy, double sz,
	// 						double groundZ,
	// 						int sc, int sr);

	// 		/**
	// 		 * Update a top with the given ID, parent and properties.
	// 		 *
	// 		 * \param id The unique ID of this top.
	// 		 * \param parentId The ID of this top's parent. The parent is
	// 		 *                 a top which is higher and whose window encompasses
	// 		 *                 the current one.
	// 		 * \param ox The original x coordinate.
	// 		 * \param oy The original y coordinate.
	// 		 * \param oz The original z coordinate.
	// 		 * \param sx The smoothed x coordinate.
	// 		 * \param sy The smoothed y coordinate.
	// 		 * \param sz The smoothed z coordinate.
	// 		 * \param sc The smoothed column.
	// 		 * \param sr The smoothed row.
	// 		 */
	// 		void update(size_t id, size_t parentId,
	// 							double ox, double oy, double oz,
	// 							double sx, double sy, double sz,
	// 							int sc, int sr);

	// 		/**
	// 		 * The smoothed x coordinate.
	// 		 *
	// 		 * \return The smoothed x coordinate.
	// 		 */
	// 		double x() const;

	// 		/**
	// 		 * The smoothed y coordinate.
	// 		 *
	// 		 * \return The smoothed y coordinate.
	// 		 */
	// 		double y() const;
	// 	};



	// 	/**
	// 	 * Save treetops stored in mapped memory into a database.
	// 	 *
	// 	 * \param topsDatabase The filename of the database.
	// 	 * \param topsLayer The layer name.
	// 	 * \param driver The database driver.
	// 	 * \param projection The WKT projection of the database.
	// 	 * \param qt A QTree containing Tops.
	// 	 * \param monitor Pointer to a Monitor object.
	// 	 */
	// 	int saveTops(const std::string& topsDatabase, const std::string& topsLayer,
	// 			const std::string& driver, const std::string& projection, mqtree<Top>& qt,
	// 			Monitor* monitor);

	// 	/**
	// 	 * Load tops from a database into mapped memory
	// 	 *
	// 	 * \param topsDatabase The filename of the database.
	// 	 * \param topsLayer The layer name.
	// 	 * \param driver The database driver.
	// 	 * \param qt A QTree containing Tops.
	// 	 * \param monitor Pointer to a Monitor object.
	// 	 */
	// 	void loadTops(const std::string& topsDatabase, const std::string& topsLayer,
	// 			mqtree<Top>& qt, Monitor* monitor);

	// 	/**
	// 	 * Returns true if the pixel at the center of the given circular window is
	// 	 * the maximum value in the window. Assigns the max pixel value to max,
	// 	 * and the proportion of nulls [0-1] to nulls.
	// 	 *
	// 	 * \param raster The source raster as vector of floats.
	// 	 * \param col The column of interest.
	// 	 * \param row The row of interest.
	// 	 * \param cols The number of columns in the raster.
	// 	 * \param rows The number of rows in the raster.
	// 	 * \param window The size of the window.
	// 	 * \param maxOut The maximum pixel value in the window.
	// 	 * \param nullsOut The proportion of pixels that are null.
	// 	 * \return True if the center pixel is the maximum.
	// 	 */
	// 	bool isMaxCenter(std::vector<float>& raster,
	// 			int col, int row, int cols, int rows, int window, float nodata,
	// 			float& maxOut, float& nullsOut);

	// 	/**
	// 	 * Returns the max value of pixels in the kernel. Sets the
	// 	 * col and row of the pixel to mr and mc, and the max to max.
	// 	 *
	// 	 * \param raster The source raster.
	// 	 * \param col The column of interest.
	// 	 * \param row The row of interest.
	// 	 * \param window The size of the window.
	// 	 * \param max The maximum value.
	// 	 * \param mc The column of the max value.
	// 	 * \param mr The row of the max value.
	// 	 */
	// 	void getKernelMax(Grid<float>& raster, int col, int row, int window,
	// 		double& max, int& mc, int& mr);

	// 	/**
	// 	 * Set all cells in the circular kernel to zero except the center.
	// 	 *
	// 	 * \param raster The source raster.
	// 	 * \param col The column of interest.
	// 	 * \param row The row of interest.
	// 	 * \param window The size of the window.
	// 	 */
	// 	void zeroKernel(Grid<int>& raster, int col, int row, int window);

	// 	/**
	// 	 * For each ID represented in the crowns raster, finds the highest pixel value
	// 	 * in the CHM within the pixels corresponding to that ID. Produces
	// 	 * a map relating the ID to a threshold containing the 3D coordinate of the
	// 	 * highest pixel.
	// 	 *
	// 	 * \param chm The CHM raster.
	// 	 * \param crowns The crowns raster.
	// 	 * \param heights Them map of crowns by Top ID.
	// 	 * \param monitor Pointer to a Monitor object.
	// 	 */
	// 	void findCrownMax(tt::grid::Grid<float>& chm, tt::grid::Grid<uint32_t>& crowns,
	// 			std::unordered_map<size_t, Crown>& heights, Monitor* monitor);

	// 	/**
	// 	 * Returns true if a pixel, represented by c, r, z is a valid crown
	// 	 * pixel, given the thresholds and the location of the treetop, represented
	// 	 * by n.
	// 	 *
	// 	 * \param c The column.
	// 	 * \param r The row.
	// 	 * \param z The elevation.
	// 	 * \param nodata The nodata value.
	// 	 * \param res The resolution.
	// 	 * \param n The Node instance pointer.
	// 	 * \param config The configuration object.
	// 	 * \param st A SimpleIntervalTree containing the IDs within height ranges (TODO: Is this correct?)
	// 	 */
	// 	bool verifyCrownPixel(int c, int r, double z, double nodata, double res,
	// 			const std::unique_ptr<Node>& n, TreetopsConfig* config,
	// 			const tt::ds::SimpleIntervalTree<double, size_t>& st);

	// 	/**
	// 	 * Returns the largest radius threshold.
	// 	 *
	// 	 * \param config The configuration to search.
	// 	 */
	// 	double maxCrownRadius(TreetopsConfig* config);

	// 	/**
	// 	 * Locates the treetop height from the original CHM. Only searches
	// 	 * within the delineated crown for the highest pixel.
	// 	 *
	// 	 * \param config The configuration object.
	// 	 */
	// 	void updateOriginalCHMHeights(TreetopsConfig* config);

	// 	/**
	// 	 * Delineate the crowns using the given treetops.
	// 	 *
	// 	 * \param config The configuration object.
	// 	 */
	// 	void delineateCrowns(TreetopsConfig* config);

	// 	/**
	// 	 * Create a vector database of the crowns.
	// 	 *
	// 	 * \param tmpFile The temporary database used to store both treetops and crowns.
	// 	 * \param config The configuration object.
	// 	 */
	// 	void polygonizeCrowns(TreetopsConfig* config);



	// 	/**
	// 	 * Locates treetops and updates their properties in
	// 	 * the treetop-tracking rasters: the window grid and the ID grid.
	// 	 *
	// 	 * \param config The configuration object.
	// 	 */
	// 	void treetopsA(TreetopsConfig* config);

	// 	/**
	// 	 * Locates the parent treetops and assigns the IDs.
	// 	 *
	// 	 * \param config The configuration object.
	// 	 */
	// 	void treetopsB(TreetopsConfig* config);

	// 	/**
	// 	 * Creates the treetop objects using the rasters as inputs.
	// 	 *
	// 	 * \param config The configuration object.
	// 	 */
	// 	void treetopsC(TreetopsConfig* config);

	// 	/**
	// 	 * Describes the possible DB conversion states.
	// 	 */
	// 	enum TTDBState {
	// 		ConvertSuccess,	///<! Copied the original SQLite database to output format.
	// 		ConvertMoved,	///<! Moved the internal DB to use as output.
	// 		ConvertNone,
	// 		Cancelled
	// 	};

	// } // util



} // tt

#endif
