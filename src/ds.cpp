#include "treetops.hpp"
#include "ds.hpp"

using namespace tt::ds;
using namespace tt::data;

Node::Node(int id, int c, int r, double z, int tc, int tr, double tz) :
    id(id),
    c(c), r(r), z(z),
    tc(tc), tr(tr), tz(tz) {
}

Node::Node(Treetop& top) :
    Node(top.id, top.col, top.row, top.height, top.col, top.row, top.height) {

}