#ifndef __LASPOINT_HPP__
#define __LASPOINT_HPP__

#include <fstream>
#include <cstdio>

#include "liblas/liblas.hpp"

#include "util.hpp"

#define LASPOINT_SIZE 79

using namespace geotools::util;

namespace geotools {

    namespace las {

        class LASPoint {
        private:

            void readLAS0(char *buf);
            void readLAS1(char *buf);
            void readLAS2(char *buf);
            void readLAS6(char *buf);
            void readLAS8(char *buf);

        public:

            double x, y, z;
            uint16_t intensity;
            uint8_t returnNum, numReturns;
            uint8_t scanDirection;
            bool isEdge;
            uint8_t cls;
            uint8_t clsFlags;
            int8_t scanAngle;
            unsigned char userData;
            uint16_t sourceId;
            double gpsTime;
            uint8_t red, green, blue, nir;
            unsigned char wavePacketDesc;
            uint64_t waveOffset;
            uint32_t wavePacketSize;
            float waveLocation;
            float xt;
            float yt;
            float zt;
            uint8_t channel;
            uint8_t format;
            
            LASPoint();

            LASPoint(const liblas::Point &pt);

            ~LASPoint();

            // Sets the scale values used by all points.
            static void setScale(double x, double y, double z);

            // Configure a liblas point from this point's values.
            void toLibLAS(liblas::Point &pt) const;

            // Configure this point from a liblas point's values.
            void fromLibLAS(const liblas::Point &pt);

            bool operator<(const LASPoint&) const;
            bool operator>(const LASPoint&) const;
            bool operator==(const LASPoint&) const;
            bool operator!=(const LASPoint&) const;
            bool operator<=(const LASPoint&) const;
            bool operator>=(const LASPoint&) const;

            void write(std::ostream &str) const;

            void read(std::istream &str);

            void write(std::FILE *str);

            void read(std::FILE *str);

            void readLAS(char *buf, uint8_t format);

            void write(void *str) const;

            void read(void *str);

            // Return true if this is a last return.
            bool last() const;

            // Return true if this is a first return.
            bool first() const;

            // Return true if this is an intermediate return.
            bool intermediate() const;

            // Return true if this is a ground point (class 2).
            bool ground() const;

            // Return true if this is a single return.
            bool single() const;

            static double scaleX();
            static double scaleY();
            static double scaleZ();

            static uint64_t dataSize();

        };
    }
}


#endif
