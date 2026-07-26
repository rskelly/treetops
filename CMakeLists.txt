cmake_minimum_required(VERSION 3.10)

project(treetops VERSION 0.1.0)

set(TREETOPS_VERSION "${PROJECT_VERSION}")

enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(GDAL 3 REQUIRED)
find_package(nlohmann_json 3 REQUIRED)

if(WIN32)
	find_library(GEOS_LIBRARY NAMES geos_c geos_c-static PATHS ${GEOS_LIBRARYDIR})
	find_library(SQLITE_LIBRARY NAMES sqlite3)
	find_library(SPATIALITE_LIBRARY NAMES spatialite libspatialite)
else()
	find_library(GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})
	set(SQLITE_LIBRARY sqlite3)
	set(SPATIALITE_LIBRARY spatialite)
endif()

if(NOT GEOS_LIBRARY)
	message(FATAL_ERROR "GEOS library (geos_c) not found")
endif()

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

if(GEOS_LIBRARYDIR OR GIS_LIBRARYDIR)
	link_directories(
		${GEOS_LIBRARYDIR}
		${GIS_LIBRARYDIR}
	)
endif()

set(TREETOPS_CORE_SOURCES
	src/config.cpp
	src/status.cpp
	src/util.cpp
	src/ds.cpp
	src/process.cpp
	src/geom.cpp
)

add_library(treetops_core STATIC ${TREETOPS_CORE_SOURCES})
set_property(TARGET treetops_core PROPERTY CXX_STANDARD 20)

if(UNIX AND NOT APPLE)
	target_compile_options(treetops_core PRIVATE -fPIC)
endif()

target_link_libraries(treetops_core PUBLIC
	${GEOS_LIBRARY}
	${GDAL_LIBRARY}
	${SQLITE_LIBRARY}
	${SPATIALITE_LIBRARY}
	nlohmann_json::nlohmann_json
)

add_executable(treetops-cli src/treetops_cli.cpp)
set_property(TARGET treetops-cli PROPERTY CXX_STANDARD 20)
target_compile_definitions(treetops-cli PRIVATE TREETOPS_VERSION="${TREETOPS_VERSION}")
target_link_libraries(treetops-cli PRIVATE treetops_core)

if(WIN32 AND VCPKG_TARGET_TRIPLET)
	set(X_VCPKG_APPLOCAL_DEPS_INSTALL ON CACHE BOOL "" FORCE)
endif()

install(TARGETS treetops-cli RUNTIME DESTINATION .)
