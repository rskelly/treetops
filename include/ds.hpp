#ifndef __DS_HPP__
#define __DS_HPP__

#include <string>
#include <queue>
#include <map>
#include <cstdint>

#include "treetops.hpp"

namespace tt {
namespace ds {
    
	/**
	 * Represents a grid cell and maintains some properties of the seed that originated it.
	 */
	class Node {
	public:
		int id;			///<! The seed (top) ID.
		int c, r;		///<! The column and row of the cell.
		double z;		///<! The z coordinate of the cell.
		int tc, tr;		///<! The column and row of the top.
		double tz;		///<! The height of the top.

		/**
		 * Construct a node using the ID and position of a Top.
		 *
		 * \param id The ID of the Top.
		 * \param c The column.
		 * \param r The row.
		 * \param z The elevation of the cell.
		 * \param tc The column of the Top.
		 * \param tr The row of the Top.
		 * \param tz The z elevation of the Top.
		 */
		Node(int id, int c, int r, double z, int tc, int tr, double tz);

		/**
		 * Construct a Node using the given Top as a source of data.
		 *
		 * \param top A Top.
		 */
		Node(tt::data::Treetop& top);
	};
        
    /**
     * Stores a single pixel column and row.
     */
    // class Px {
    // public:
    //     Px(uint32_t col, uint32_t row);
    //     uint32_t m_col;   // The column coordinate.
    //     uint32_t m_row;   // Thr row coordinate.
    // };

    /**
     * Represents a single tree top, with a unique ID, column,
     * row and pixel queue. A count
     */
    // class Tt {
    // public:
    //     /**
    //      * Construct a treetop with the given column and row. The ID is generated
    //      * automatically.
    //      */
    //     Tt(uint32_t col, uint32_t row);
    //     uint32_t m_id;                // The unique ID of the tree top.
    //     uint32_t m_col;               // The column coordinate.
    //     uint32_t m_row;               // Thr row coordinate.
    //     uint32_t m_px_count;          // The number of pixels in the queue remaining from the previous search.
    //     std::queue<Px> m_px_queue;    // The pixel queue.
    // };

    /**
     * The configuration. 
     */
    // class Config {
    // public:
    //     /**
    //      * Construct the configuration and set default values.
    //      */
    //     Config();
    //     /**
    //      * Construct a config object using properties from the given map.
    //      * Map values will be coerced to the appropriate types. Default
    //      * properties will maintain their original values if not set.
    //      */
    //     Config(const std::map<std::string, std::string>&);
    //     static void saveToFile(const Config&, const std::string&);
    //     std::string m_input_file;         // The path to the canopy DSM raster.
    //     std::string m_crowns_raster;      // The path to the tree crowns raster.
    //     std::string m_smoothed_raster;    // The path to the Gaussian smoothed raster,
    //     std::string m_crowns_db;          // The path to the crowns database (vector).
    //     std::string m_tops_db;            // The path to the tops database (vector).
    //     std::string m_config;             // The path to the JSON configuration file.
    //     float m_gauss_std;                // The standard deviation for the Gaussian smoother.
    //     float m_gauss_radius;             // The radius for the Gaussian smoother.
    //     float m_tt_search_radius;         // The tree top search radius.
    //     float m_crown_radius;             // The crown radius limit.
    //     float m_crown_fraction;           // The crown height fraction.
    //     bool m_trim_dangles;              // Trim dangles from polygonized crowns.
    //     bool m_remove_holes;              // Remove holes from polygonized crowns.
    // };

}
}

#endif