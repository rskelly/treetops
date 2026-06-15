cmake_minimum_required(VERSION 3.10)

project(treetops)

enable_language(C)
enable_language(CXX)

set(THREADS_PREFER_PTHREAD_FLAG ON)

set(SQLITE_LIBRARY sqlite3)
set(SPATIALITE_LIBRARY spatialite)

find_package(GDAL 3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(nlohmann_json 3 REQUIRED)
find_library(GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})

option(TREETOPS_BUILD_GUI "Build the Qt GUI application" OFF)

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
	src/fitpack/fitpack.cpp
	src/config.cpp
	src/util.cpp
	src/ds.cpp
	src/process.cpp
	src/vector.cpp
	src/treetops.cpp
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

if(TREETOPS_BUILD_GUI)
	find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

	set(CMAKE_AUTOUIC ON)
	set(CMAKE_AUTOMOC ON)
	set(CMAKE_AUTORCC ON)
	set(CMAKE_AUTOUIC_SEARCH_PATHS ${CMAKE_CURRENT_SOURCE_DIR}/src/ui)

	include_directories(BEFORE SYSTEM ${Qt5Widgets_INCLUDE_DIRS})

	add_executable(treetops-gui
		include/settings.hpp
		src/settings.cpp
		src/ui/ui_util.cpp
		src/ui/treetops.ui
		src/ui/tops_thresholds.ui
		src/ui/crowns_thresholds.ui
		include/ui/treetops_ui.hpp
		src/ui/treetops_ui.cpp
		include/ui/tops_threshold_item.hpp
		src/ui/tops_threshold_item.cpp
		include/ui/tops_thresholds_ui.hpp
		src/ui/tops_thresholds_ui.cpp
		include/ui/crowns_threshold_item.hpp
		src/ui/crowns_threshold_item.cpp
		include/ui/crowns_thresholds_ui.hpp
		src/ui/crowns_thresholds_ui.cpp
		src/treetops_app.cpp
	)

	set_property(TARGET treetops-gui PROPERTY CXX_STANDARD 20)
	target_compile_options(treetops-gui PRIVATE -fPIC)
	target_link_libraries(treetops-gui PRIVATE treetops_core Qt5::Widgets)
endif()
