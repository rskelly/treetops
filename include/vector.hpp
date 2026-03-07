#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__

#include <geos_c.h>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>

#include "settings.hpp"
#include "treetops.hpp"
#include "util.hpp"

using namespace tt::data;

namespace tt {
namespace vec {

    /**
     * A fast, sequential database for storing Crown and Treetop instances.
     *
     * Each top record is constant length, stored as binary. The poly_position field gives
     * the location of the polygon in the crowns file. If poly_len is zero, there is no crown,
     * if it is >0, the crown is read from the file at poly_position.
     *
     * id (int) | top (Top) | poly_len (int) | poly_position (int)
     *
     * The crowns file contains crown polygons. This file has no structure. Crowns
     * must be located and read using the poly_len and poly_position values in the
     * tops file. The crowns are WKB.
     */
    class CrownDB {
    private:
        int m_tfd; 								///<! File descriptor of open tops database.
        int m_cfd; 								///<! File descriptor of open crowns database.
        int m_toffset;						///<! Current offset in the tops file.
        int m_coffset;						///<! Current offset in the crowns file.
        int m_tidx;							///<! The iterator index.
        int m_count;
        GEOSContextHandle_t m_gctx;				///<! GEOS context.
        GEOSWKBReader* m_rdr;					///<! WKB Reader.
        GEOSWKBWriter* m_wtr;					///<! WKB Writer.
        std::string m_tfile;					///<! The tops database filename.
        std::string m_cfile;					///<! The crowns database filename.
        std::unordered_map<int, int> m_idx; 	///<! Mapping from top ID to offset in the file.

        std::mutex m_imtx;
        std::mutex m_tmtx;
        std::mutex m_cmtx;

        /**
         * \brief Write the crown to the crowns file and return its position.
         *
         * If the length of the geometry is greater than the original
         * geometry (if it exists) it'll be written at the end and the original
         * will be orphaned.
         *
         * The new length and position are set into the plen and ppos values.
         *
         * True is returned if the process succeeds.
         */
        bool writeCrown(int& plen, int& ppos, const GEOSGeometry* crown, GEOSContextHandle_t gctx);

    public:

        /**
         * \brief Create the file DB using temporary files.
         */
        CrownDB();

        /**
         * \brief Return the internal GEOS context handle.
         */
        GEOSContextHandle_t gctx();

        /**
         * \brief Attempt to load the Top and crown for the given ID.
         *
         * \param id The top ID.
         * \param[out] top A reference to a writable Top.
         * \param[out] crown A reference to a writable pointer to the crown geometry.
         * \param gctx A GEOSContextHandle_t.
         * \param withCrown If true, loads the crown geometry. Otherwise not.
         * \return True if the record was found.
         */
        bool find(int id, tt::data::Treetop& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx = nullptr, bool withCrown = true);

        /**
         * \brief Read the next available Top and crown geometry.
         *
         * \param[inout] top A Top.
         * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
         * \param gctx An optional GEOS context.
         * \param loadCrown If false, the crown is not loaded.
         * \return True if a record was found.
         */
        bool next(Treetop& top, GEOSGeometry*& crown, GEOSContextHandle_t gctx = nullptr, bool loadCrown = false);

        /**
         * \brief Read the next available Top and crown geometry.
         *
         * \param[inout] top A Top.
         * \param[inout] crownwkb A char buffer representing the crown's WKB representation.
         * \param loadCrown If false, the crown is not loaded.
         * \return True if a record was found.
         */
        bool next(Treetop& top, std::vector<unsigned char>& crownwkb, bool loadCrown = false);

        bool next(Treetop& top);

        int count();

        /**
         * \brief Reset the iterator.
         */
        void reset();

        /**
         * \brief Updates or adds the Treetop and Crown.
         *
         * If the Top isn't already in the DB it will be inserted at the end.
         *
         * \param[inout] top A Treetop.
         * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
         * \return True on success.
         */
        bool insert(const Treetop& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx = nullptr);

        /**
         * \brief Update the Treetop and Crown.
         *
         * \param[inout] top A Treetop.
         * \param[inout] crown A GEOSGeometry representing a crown. The caller is responsible for freeing the geometry.
         * return True on success.
         */
        bool update(const Treetop& top, const GEOSGeometry* crown, GEOSContextHandle_t gctx = nullptr);

        /**
         * \brief Save the Tops into the given database file and layer.
         */
        void saveTops(const std::string& file, const std::string& driver,
                const std::string& layer, const std::string& projection);

        /**
         * \brief Save the crowns with geometry into the given database file and layer.
         */
        void saveCrowns(const std::string& file, const std::string& driver,
                const std::string& layer, const std::string& projection);

        /**
         * \brief Return the filesize of the crowns file in bytes.
         */
        int crownsDBSize();

        /**
         * \brief Return the filesize of the tops file in bytes.
         */
        int topsDBSize();

        /**
         * \brief Return true if the tops database can be saved using the given driver.
         */
        bool checkSaveTops(const std::string& driver);

        /**
         * \brief Return true if the crowns database can be saved using the given driver.
         */
        bool checkSaveCrowns(const std::string& driver);

        ~CrownDB();
    };

}
}

#endif // __VECTOR_HPP__