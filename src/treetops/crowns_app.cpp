#include <iostream>

#include "crowns.hpp"

using namespace dijital::crowns;
using namespace dijital::crowns::config;

int main(int argc, char** argv) {

    if(argc < 2) {
        std::cerr << "A configuration file is required.";
        return 1;
    }
    
    std::string config_file(argv[1]);


    Crowns crowns;
    CrownsAppConfig& config = crowns.config();
    config.loadFromJSON(config_file);
    crowns.run();

    return 0;
}