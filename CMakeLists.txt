cmake_minimum_required(VERSION 3.22)

project(treetops VERSION 1.0 LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable automatic handling of Qt MOC, UIC, RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Optionally set CMAKE_PREFIX_PATH if Qt is installed in a custom location
if(NOT DEFINED CMAKE_PREFIX_PATH)
    set(CMAKE_PREFIX_PATH "$ENV{QT_ROOT}" CACHE PATH "Qt install path")
endif()

# Find Qt 6 packages
find_package(Qt6 COMPONENTS Widgets Gui Core REQUIRED)

# Optionally find GIS/math libraries (GDAL, PROJ, GEOS)
find_package(GDAL CONFIG QUIET)
find_package(PROJ CONFIG QUIET)
find_package(GEOS 3 CONFIG QUIET)
find_package(PDAL CONFIG QUIET)

include_directories(BEFORE SYSTEM include)

# Collect all source and header files
file(GLOB_RECURSE SOURCES src/*.cpp)
file(GLOB_RECURSE HEADERS include/*.hpp)

# Add the executable
add_executable(treetops ${SOURCES} ${HEADERS})

# Link Qt libraries
target_link_libraries(treetops PRIVATE Qt6::Widgets Qt6::Gui Qt6::Core)

# Link GIS/math libraries if found
if(GDAL_FOUND)
    target_link_libraries(treetops PRIVATE GDAL::GDAL)
endif()

if(PROJ_FOUND)
    target_link_libraries(treetops PRIVATE PROJ::proj)
endif()

#if(GEOS_FOUND)
    target_link_libraries(treetops PRIVATE GEOS::geos GEOS::geos_c)
#endif()

if(PDAL_FOUND)
    target_link_libraries(treetops PRIVATE GEOS::geos)
endif()

# Optional: set output directory for the executable
set_target_properties(treetops PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Include directories (optional, e.g., for GIS headers)
if(GDAL_FOUND)
    target_include_directories(treetops PRIVATE ${GDAL_INCLUDE_DIRS})
endif()

if(PROJ_FOUND)
    target_include_directories(treetops PRIVATE ${PROJ_INCLUDE_DIRS})
endif()

if(GEOS_FOUND)
    target_include_directories(treetops PRIVATE ${GEOS_INCLUDE_DIRS})
endif()

# Optional: enable testing if you have CTest configured
enable_testing()

