cmake_minimum_required(VERSION 3.10)

project(treetops)

enable_language(CXX)

set(SQLITE_LIBRARY sqlite3)
set(SPATIALITE_LIBRARY spatialite)

find_package(GDAL 3 REQUIRED)
find_package(nlohmann_json 3 REQUIRED)
find_library(GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})

execute_process(COMMAND git rev-parse --short HEAD
	RESULT_VARIABLE GIT_RESULT
	OUTPUT_VARIABLE GIT_OUTPUT
	ERROR_VARIABLE GIT_ERROR)
string(STRIP "${GIT_OUTPUT}" GIT_REV)
message(STATUS "Finding Git rev: ${GIT_REV} (Result: ${GIT_RESULT}, Error: ${GIT_ERROR})")

include_directories(BEFORE SYSTEM
	include
	${GDAL_INCLUDE_DIR}
	${GEOS_INCLUDEDIR}
)

link_directories(
	${GEOS_LIBRARYDIR}
	${GIS_LIBRARYDIR}
)

set(TREETOPS_CORE_SOURCES
	src/config.cpp
	src/util.cpp
	src/ds.cpp
	src/process.cpp
	src/vector.cpp
)

add_library(treetops_core STATIC ${TREETOPS_CORE_SOURCES})
set_property(TARGET treetops_core PROPERTY CXX_STANDARD 20)
target_compile_options(treetops_core PRIVATE -fPIC)
target_link_libraries(treetops_core PUBLIC
	${GEOS_LIBRARY}
	${GDAL_LIBRARY}
	${SQLITE_LIBRARY}
	${SPATIALITE_LIBRARY}
	nlohmann_json::nlohmann_json
)

add_executable(treetops-cli src/treetops_cli.cpp)
set_property(TARGET treetops-cli PROPERTY CXX_STANDARD 20)
target_link_libraries(treetops-cli PRIVATE treetops_core)
