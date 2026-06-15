#ifndef __STATUS_HPP__
#define __STATUS_HPP__

#include <functional>
#include <string>

namespace tt {
namespace status {

struct Update {
	std::string phase;
	int progress = 0;
	std::string message;
};

using Callback = std::function<void(const Update&)>;

void setCallback(Callback callback);
void clearCallback();
void report(const std::string& phase, int progress, const std::string& message);
void report(const Update& update);

} // status
} // tt

#endif
