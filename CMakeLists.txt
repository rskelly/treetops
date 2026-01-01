cmake_minimum_required(VERSION 3.10)

project (treetops)

enable_language(C)
enable_language(CXX)

set (THREADS_PREFER_PTHREAD_FLAG ON)

set (SQLITE_LIBRARY sqlite3)
set (SPATIALITE_LIBRARY spatialite)

find_package (GDAL 3 REQUIRED)
find_package (OpenMP REQUIRED)
find_package (Qt5 COMPONENTS Core Widgets REQUIRED)
find_package (nlohmann_json 3 REQUIRED)

find_library (GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})

# if(OPENMP_FOUND)
# 	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
# 	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
# endif()

# Get the git rev and create an environment variable to be used in compilation.

execute_process (COMMAND git rev-parse --short HEAD RESULT_VARIABLE GIT_RESULT 
	OUTPUT_VARIABLE GIT_OUTPUT ERROR_VARIABLE GIT_ERROR)
string (STRIP "${GIT_OUTPUT}" GIT_REV)
message (STATUS "Finding Git rev: ${GIT_REV} (Result: ${GIT_RESULT}, Error: ${GIT_ERROR})")

set(CMAKE_AUTOUIC_SEARCH_PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui
)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

include_directories(BEFORE SYSTEM
	include
	${Qt5Widgets_INCLUDE_DIRS}
	${EIGEN_INCLUDEDIR}	
	${GDAL_INCLUDE_DIR}
	${GEOS_INCLUDEDIR} 
)
	 
link_directories(
	${GEOS_LIBRARYDIR} 
	${GIS_LIBRARYDIR} 
)
	
add_executable(treetops 
	src/fitpack/fitpack.cpp
	src/settings.cpp
	src/util.cpp
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
#	src/grid.cpp
#	src/filelist.cpp
    src/treetops.cpp
)

target_compile_options(treetops PRIVATE
    -fPIC
)

target_link_libraries(treetops PRIVATE
 	#Qt5::Gui 
	#Qt5::Core 
	Qt5::Widgets 
	${GDAL_LIBRARY}
 	${SQLITE_LIBRARY} 
 	${SPATIALITE_LIBRARY} 
)
