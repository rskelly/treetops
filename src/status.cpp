#include "status.hpp"

#include <iostream>
#include <mutex>

#include <nlohmann/json.hpp>

namespace tt {
namespace status {

namespace {

std::mutex g_mutex;
Callback g_callback;

void defaultCallback(const Update& update) {
	nlohmann::json payload{
		{"phase", update.phase},
		{"progress", update.progress},
		{"message", update.message}
	};
	std::cout << "@status:" << payload.dump() << std::endl;
}

} // namespace

void setCallback(Callback callback) {
	std::lock_guard<std::mutex> lock(g_mutex);
	g_callback = std::move(callback);
}

void clearCallback() {
	std::lock_guard<std::mutex> lock(g_mutex);
	g_callback = nullptr;
}

void report(const std::string& phase, int progress, const std::string& message) {
	report(Update{phase, progress, message});
}

void report(const Update& update) {
	Callback callback;
	{
		std::lock_guard<std::mutex> lock(g_mutex);
		callback = g_callback ? g_callback : Callback(defaultCallback);
	}
	callback(update);
}

} // status
} // tt
