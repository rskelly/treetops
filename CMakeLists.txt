cmake_minimum_required(VERSION 3.1)

set (WITH_GUI true CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")

set (Boost_DEBUG OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_STATIC_RUNTIME OFF)

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
elseif(${WIN32})
	set (Boost_COMPILER -vc141)
	set (Boost_NO_SYSTEM_PATHS  ON)
	set (DEV_DIR "C:\\dev")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\gisinternals")
	set (LIBGEO_LIBRARYDIR "${GIT_DIR}\\treetops\\libgeo\\build")
	set (LIBGEO_INCLUDEDIR "${GIT_DIR}\\treetops\\libgeo\\include")
	set (GEOS_LIBRARYDIR "${DEV_DIR}\\geos\\lib")
	set (GEOS_INCLUDEDIR "${DEV_DIR}\\geos\\include")
	set (BOOST_INCLUDEDIR "${DEV_DIR}\\boost_1_67_0")
	set (BOOST_LIBRARYDIR "${DEV_DIR}\\boost_1_67_0\\stage\\lib")
	set (EIGEN_INCLUDEDIR "${DEV_DIR}\\eigen3")
	set (Boost_INCLUDE_DIR "${BOOST_INCLUDEDIR}")
	set (GDAL_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_LIBRARYDIR "${GIS_DIR}\\lib")
	set (libLAS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (libLAS_LIBRARY_DIRS "${GIS_DIR}\\lib")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GIS_DIR}\\bin;${DEV_DIR}\\Qt\\5.11.0\\msvc2017_64\\lib\\cmake")
	set (SQLITE_LIBRARY sqlite3_i)
	set (SPATIALITE_LIBRARY spatialite_i)
	set (GDAL_LIBRARY gdal_i)
	include_directories(BEFORE SYSTEM .\\build\\ui_util_autogen\\include ${LIBGEO_INCLUDEDIR} ${EIGEN_INCLUDEDIR} ${GEOS_INCLUDEDIR})
	link_directories(${LIBGEO_LIBRARYDIR} ${GEOS_LIBRARYDIR} ${GEOS_LIBRARYDIR} ${BOOST_LIBRARYDIR} ${GDAL_LIBRARYDIR})
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
	include_directories (BEFORE SYSTEM /usr/include/eigen3)
	include_directories (BEFORE SYSTEM /opt/include)
	link_directories (BEFORE SYSTEM /opt/lib)
endif()

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)

find_library (GEOS_LIBRARY geos PATHS ${GEOS_LIBRARYDIR})
 
if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

if(${WITH_GUI})
	set (CMAKE_AUTOUIC ON)
	set (CMAKE_AUTOMOC ON)
	set (CMAKE_INCLUDE_CURRENT_DIR ON)
	set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")
endif()

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)

find_library(GEOS_LIBRARY geos PATHS ${GEOS_LIBRARYDIR})

if(${WITH_GUI})
	find_package(Qt5Widgets REQUIRED)
endif()

project (treetops)
add_subdirectory (libtreetops)
add_subdirectory (libgeo)

# Configure directories

link_directories(. ${Boost_LIBRARY_DIRS})

# Build support libraries ############################################################################################

if(${WITH_GUI})
	add_library (ui_util STATIC src/ui_util.cpp src/crs_selector_ui.cpp src/tops_thresholds_ui.cpp src/crowns_thresholds_ui.cpp)
	target_include_directories(ui_util PUBLIC libgeo/include ui_util_autogen/include include  ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} 
		${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
	target_link_libraries(ui_util Qt5::Widgets)
endif()


# Build programs #####################################################################################################

if(${WITH_GUI})
	add_executable (treetops-app src/treetops.cpp src/treetops_ui.cpp)
	target_include_directories(treetops-app PUBLIC libgeo/include include ${INC_PATH} ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} 
		${INC_PATH} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
	target_link_libraries (treetops-app PUBLIC Qt5::Gui Qt5::Core Qt5::Widgets treetops ui_util)
else()
	add_executable (treetops-app src/treetops.cpp src/treetops_ui.hpp)
	target_include_directories(treetops-app PUBLIC libgeo/include include ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} 
		${INC_PATH} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
	target_link_libraries (treetops-app PUBLIC treetops)
endif()

# Treetops has the same name as its library. Rename it as a last step
# so there's no conflict.
set_target_properties (treetops-app PROPERTIES OUTPUT_NAME treetops)


# Install ###########################################################################################################

install (TARGETS treetops-app DESTINATION bin)
