#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

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
		groundZ(0),
		sx(0), sy(0), sz(0),
		parentId(0) {}
};

} // data
} // tt

#endif
