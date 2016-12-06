#include <fstream>
#include <cstdio>

#include "liblas/liblas.hpp"

#include "geotools.hpp"
#include "laspoint.hpp"
#include "util.hpp"

// The number of points to store in a row before writing.
// TODO: Configure/optimize.
#define CACHE_LEN 2048

// The size in bytes of the header row.
#define HEADER_LEN 64

using namespace geotools::las;
using namespace geotools::util;

double las_scaleX = 0, las_scaleY = 0, las_scaleZ = 0;

LASPoint::LASPoint() :
		x(0), y(0), z(0),
		intensity(0),
		returnNum(0), numReturns(0),
		scanDirection(0),
		isEdge(false),
		cls(0), clsFlags(0),
		scanAngle(0),
		userData(0),
		sourceId(0),
		gpsTime(0),
		red(0), green(0), blue(0), nir(0),
		wavePacketDesc(0),
		waveOffset(0),
		wavePacketSize(0),
		waveLocation(0),
		xt(0), yt(0), zt(0),
		channel(0),
		format(0) {
}

LASPoint::LASPoint(const liblas::Point &pt) {
	fromLibLAS(pt);
}

void LASPoint::setScale(double x, double y, double z) {
	las_scaleX = x;
	las_scaleY = y;
	las_scaleZ = z;
}

double LASPoint::scaleX() {
	return las_scaleX;
}

double LASPoint::scaleY() {
	return las_scaleY;
}

double LASPoint::scaleZ() {
	return las_scaleZ;
}

uint64_t LASPoint::dataSize() {
	return LASPOINT_SIZE;
}

bool LASPoint::operator<(const LASPoint &p) const {
	//	g_debug(" -- laspoint < " << y << " < " << p.y << " -> " << (y < p.y));
	return y < p.y;
}

bool LASPoint::operator>(const LASPoint &p) const {
	return y > p.y;
}

bool LASPoint::operator==(const LASPoint &p) const {
	return y == p.y && x == p.x && z == p.z;
}

bool LASPoint::operator<=(const LASPoint &p) const {
	return y <= p.y;
}

bool LASPoint::operator>=(const LASPoint &p) const {
	return y >= p.y;
}

bool LASPoint::operator!=(const LASPoint &p) const {
	return y != p.y && x != p.x && z != p.z;
}

void LASPoint::fromLibLAS(const liblas::Point &pt) {
	x = pt.GetX();
	y = pt.GetY();
	z = pt.GetZ();
	intensity = pt.GetIntensity();
	returnNum = pt.GetReturnNumber();
	numReturns = pt.GetNumberOfReturns();
	cls = pt.GetClassification().GetClass();
	scanAngle = pt.GetScanAngleRank();
	isEdge = pt.GetFlightLineEdge();
	sourceId = pt.GetPointSourceID();
	scanDirection = pt.GetScanDirection();
	userData = pt.GetUserData();
	if (format == 1 || format > 2)
		gpsTime = pt.GetTime();
	if (format == 2 || format == 3 || format == 7 || format == 8
			|| format == 10) {
		const liblas::Color &c = pt.GetColor();
		red = c.GetRed();
		green = c.GetGreen();
		blue = c.GetBlue();
		// nir = c.?
	}
}

void LASPoint::toLibLAS(liblas::Point &pt) const {
	pt.SetX(x);
	pt.SetY(y);
	pt.SetZ(z);
	pt.SetIntensity(intensity);
	pt.SetReturnNumber(returnNum);
	pt.SetNumberOfReturns(numReturns);
	pt.SetScanDirection((uint16_t) scanDirection);
	pt.SetFlightLineEdge((uint16_t) isEdge);
	pt.SetClassification(liblas::Classification(cls, false, false, false));
	pt.SetScanAngleRank(scanAngle);
	pt.SetUserData(userData);
	pt.SetPointSourceID(sourceId);
	if (format == 1 || format > 2)
		pt.SetTime(gpsTime);
	if (format == 2 || format == 3 || format == 7 || format == 8
			|| format == 10)
		pt.SetColor(liblas::Color(red, green, blue)); // no nir
}

void LASPoint::read(std::istream &str) {
	int32_t xx, yy, zz;
	str.read((char *) &xx, sizeof(int32_t));
	str.read((char *) &yy, sizeof(int32_t));
	str.read((char *) &zz, sizeof(int32_t));
	str.read((char *) &intensity, sizeof(uint16_t));
	str.read((char *) &returnNum, sizeof(uint16_t));
	str.read((char *) &numReturns, sizeof(uint16_t));
	str.read((char *) &cls, sizeof(uint8_t));
	str.read((char *) &scanAngle, sizeof(int8_t));
	x = (double) (xx * las_scaleX);
	y = (double) (yy * las_scaleY);
	z = (double) (zz * las_scaleZ);
}

void LASPoint::write(std::ostream &str) const {
	int32_t xx = (int32_t) (x / las_scaleX);
	int32_t yy = (int32_t) (y / las_scaleY);
	int32_t zz = (int32_t) (z / las_scaleZ);
	str.write((char *) &xx, sizeof(int32_t));
	str.write((char *) &yy, sizeof(int32_t));
	str.write((char *) &zz, sizeof(int32_t));
	str.write((char *) &intensity, sizeof(uint16_t));
	str.write((char *) &returnNum, sizeof(uint16_t));
	str.write((char *) &numReturns, sizeof(uint16_t));
	str.write((char *) &cls, sizeof(uint8_t));
	str.write((char *) &scanAngle, sizeof(int8_t));
}

void LASPoint::read(std::FILE *str) {
	Buffer buffer(LASPoint::dataSize());
	if (!std::fread(buffer.buf, LASPoint::dataSize(), 1, str))
		g_runerr("Nothing read from file.");
	read(buffer.buf);
}

void LASPoint::readLAS0(char *buf) {
	int32_t xx, yy, zz;
	uint8_t props, cprops;
	xx = *((int32_t *) buf);
	buf += sizeof(int32_t);
	yy = *((int32_t *) buf);
	buf += sizeof(int32_t);
	zz = *((int32_t *) buf);
	buf += sizeof(int32_t);
	intensity = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	props = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	cprops = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	scanAngle = *((int8_t *) buf);
	buf += sizeof(int8_t);
	userData = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	sourceId = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	x = las_scaleX * xx;
	y = las_scaleY * yy;
	z = las_scaleZ * zz;
	returnNum = props & 6;
	numReturns = (props >> 3) & 6;
	scanDirection = ((props >> 6) & 1) == 1;
	isEdge = ((props >> 7) & 1) == 1;
	cls = cprops & 31;
	clsFlags = cprops >> 5;
}

void LASPoint::readLAS1(char *buf) {
	gpsTime = *((double *) buf);
	buf += sizeof(double);
}

void LASPoint::readLAS2(char *buf) {
	red = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	green = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	blue = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
}

void LASPoint::readLAS8(char *buf) {
	nir = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
}

void LASPoint::readLAS6(char *buf) {
	int32_t xx, yy, zz;
	uint8_t props, cprops;
	xx = *((int32_t *) buf);
	buf += sizeof(int32_t);
	yy = *((int32_t *) buf);
	buf += sizeof(int32_t);
	zz = *((int32_t *) buf);
	buf += sizeof(int32_t);
	intensity = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	props = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	cprops = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	cls = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	userData = *((uint8_t *) buf);
	buf += sizeof(uint8_t);
	scanAngle = *((int16_t *) buf);
	buf += sizeof(int16_t);
	sourceId = *((uint16_t *) buf);
	buf += sizeof(uint16_t);
	gpsTime = *((double *) buf);
	buf += sizeof(double);
	x = las_scaleX * xx;
	y = las_scaleY * yy;
	z = las_scaleZ * zz;
	returnNum = props & 15;
	numReturns = (props >> 4) & 15;
	clsFlags = cprops & 15;
	channel = (cprops >> 4) & 3;
	scanDirection = ((cprops >> 6) & 1) == 1;
	isEdge = ((cprops >> 7) & 1) == 1;

}

void LASPoint::readLAS(char *buf, uint8_t format) {
	this->format = format;
	switch (format) {
	case 0:
		readLAS0(buf);
		break;
	case 1:
		readLAS0(buf);
		readLAS1(buf);
		break;
	case 2:
		readLAS0(buf);
		readLAS2(buf);
		break;
	case 3:
		readLAS0(buf);
		readLAS1(buf);
		readLAS2(buf);
		break;
	case 6:
		readLAS6(buf);
		break;
	case 7:
		readLAS6(buf);
		readLAS2(buf);
		break;
	case 8:
		readLAS6(buf);
		readLAS2(buf);
		readLAS8(buf);
		break;
	case 4:
	case 5:
	case 9:
	case 10:
		g_runerr("Waveform not implemented.");
		break;
	default:
		g_runerr("Invalid format " << format);
		break;
	}
}

void LASPoint::write(std::FILE *str) {
	Buffer buffer(LASPoint::dataSize());
	write(buffer.buf);
	std::fwrite(buffer.buf, LASPoint::dataSize(), 1, str);
}

void LASPoint::read(void *str) {
	char *ptr = (char *) str;
	x = *((int32_t *) ptr) * las_scaleX;
	ptr += sizeof(int32_t);
	y = *((int32_t *) ptr) * las_scaleY;
	ptr += sizeof(int32_t);
	z = *((int32_t *) ptr) * las_scaleZ;
	ptr += sizeof(int32_t);
	intensity = *((uint16_t *) ptr);
	ptr += sizeof(uint16_t);
	returnNum = *((uint16_t *) ptr);
	ptr += sizeof(uint16_t);
	numReturns = *((uint16_t *) ptr);
	ptr += sizeof(uint16_t);
	cls = *((uint8_t *) ptr);
	ptr += sizeof(uint8_t);
	scanAngle = *((int8_t *) ptr);
}

void LASPoint::write(void *str) const {
	char *ptr = (char *) str;
	*((int32_t *) ptr) = (int32_t) (x / las_scaleX);
	ptr += sizeof(int32_t);
	*((int32_t *) ptr) = (int32_t) (y / las_scaleY);
	ptr += sizeof(int32_t);
	*((int32_t *) ptr) = (int32_t) (z / las_scaleZ);
	ptr += sizeof(int32_t);
	*((uint16_t *) ptr) = intensity;
	ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = returnNum;
	ptr += sizeof(uint16_t);
	*((uint16_t *) ptr) = numReturns;
	ptr += sizeof(uint16_t);
	*((uint8_t *) ptr) = cls;
	ptr += sizeof(uint8_t);
	*((int8_t *) ptr) = scanAngle;
}

bool LASPoint::last() const {
	return numReturns > 0 && returnNum == numReturns;
}

bool LASPoint::first() const {
	return numReturns > 0 && returnNum == 1;
}

bool LASPoint::intermediate() const {
	return numReturns > 2 && returnNum > 1 && returnNum < numReturns;
}

bool LASPoint::ground() const {
	return cls == 2;
}

bool LASPoint::single() const {
	return numReturns == 1;
}

LASPoint::~LASPoint() {
}

