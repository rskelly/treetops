#include "treetops.hpp"
#include "process.hpp"

using namespace tt;

Treetops::Treetops(tt::config::Config& settings) :
    m_settings(&settings) {}

void Treetops::smooth() {
    tt::proc::Processor processor(*m_settings);
    // Smoothing is run as part of the full pipeline when doSmoothing is enabled.
    processor.run();
}

void Treetops::treetops() {
    m_settings->set("doSmoothing", false);
    m_settings->set("doCrowns", false);
    tt::proc::Processor(*m_settings).run();
}

void Treetops::treecrowns() {
    m_settings->set("doSmoothing", false);
    m_settings->set("doTops", false);
    tt::proc::Processor(*m_settings).run();
}

Treetops::~Treetops() {}
