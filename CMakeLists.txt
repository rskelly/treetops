cmake_minimum_required(VERSION 3.1)

set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (WITH_GUI true CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

set (GEOS_LIBRARY geos)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -D_GLIBCXX_PARALLEL")

if (${APPLE})

	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH "/usr/local/Cellar/qt5/5.6.1-1/lib/cmake/Qt5Widgets")
	endif()

	# For GDAL2
	set (GDAL_LIBRARY /usr/local/opt/gdal2/lib/libgdal.dylib)
	set (GDAL_INCLUDE_DIR /usr/local/opt/gdal2/include)
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")

elseif(${WIN32})
	
	# Preparing to build on Windows
	# 1) download and install vs 2015
	# 2) use web installer to install qt5
	# 3) download prebuild boost for vs msvc-14 -- build using bootstrap, b2
	# 3) download gisinternals lib and headers zips for msvc-12 -- unzip to C:\opt
	# 4) download and install cmake
	# 5) download and install git
	# --- 6) checkout liblas from git
	# 7) download eigen -- extract to c:\opt\include -- rename folder to eigen3
	
	# Configure boost.
	set (Boost_COMPILER -vc140)
	set (Boost_NO_SYSTEM_PATHS  ON)
	set (BOOST_INCLUDEDIR "E:\\boost_1_63_0")
	set (BOOST_LIBRARYDIR "E:\\boost_1_63_0\\lib64-msvc-14.0")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:\\Qt\\5.7\\msvc2015_64\\lib\\cmake\\Qt5Widgets")
	set (LIB_PATH "E:\\opt\\lib")
	set (INC_PATH "E:\\opt\\include")
	set (GDAL_INCLUDE_DIR "E:\\opt\\include")
	set (GDAL_LIBRARY "E:\\opt\\lib\\gdal_i.lib")
else()
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:/usr/lib/x86_64-linux-gnu/cmake/Qt5")
		#set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:/home/rob/Qt5.7.0/5.7/gcc_64/lib/cmake/Qt5Widgets/")
	endif()
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude")

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)
find_package (libLAS 1.8 REQUIRED)

if(${WITH_GUI})
	find_package(Qt5Widgets REQUIRED)
endif()

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

project (geotools)

# Configure directories
link_directories(. ${Boost_LIBRARY_DIRS} ${LIB_PATH})
include_directories(include ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${INC_PATH})


# Build support libraries ############################################################################################

add_library (geotools STATIC src/geotools.cpp)

add_library (util STATIC src/util.cpp)
target_link_libraries (util ${Boost_LIBRARIES} geotools)

add_library (db STATIC src/db.cpp)
target_link_libraries (db ${GDAL_LIBRARY} geotools)

add_library (raster STATIC src/raster.cpp)
target_link_libraries (raster ${GDAL_LIBRARY} util geotools db)

add_library (cellstats STATIC src/cellstats.cpp src/laspoint.cpp)
target_link_libraries (cellstats gmp CGAL boost_system las geotools)

add_library (lasreader STATIC src/lasreader.cpp src/laspoint.cpp)
target_link_libraries (lasreader gmp CGAL boost_system las geotools raster)

if(${WITH_GUI})
	add_library (ui_util STATIC src/apps/ui_util.cpp)
	target_link_libraries(ui_util Qt5::Widgets)
endif()


# Build programs #####################################################################################################

#add_executable (flood src/flood.cpp)
#target_link_libraries (flood ${GDAL_LIBRARY} raster util geotools)


#add_executable (rasterstats src/rasterstats.cpp)
#target_link_libraries (rasterstats ${GDAL_LIBRARY} raster util geotools)


#add_library (mosaic SHARED src/mosaic.cpp)
#target_link_libraries (mosaic ${GDAL_LIBRARY} ${Boost_LIBRARIES} raster util geotools)

#if(${WITH_GUI})
#	add_executable (mosaic-app src/apps/mosaic.cpp src/apps/mosaic_ui.cpp)
#	target_link_libraries (mosaic-app PUBLIC mosaic Qt5::Gui Qt5::Core Qt5::Widgets)
#else()
#	add_executable (mosaic-app src/apps/mosaic.cpp)
#	target_link_libraries (mosaic-app PUBLIC mosaic)
#endif()


add_library (treetops SHARED src/treetops.cpp)
target_link_libraries (treetops ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${GEOS_LIBRARY} ${SQLITE3_LIBRARY} ${SPATIALITE_LIBRARY} raster util geotools db)

if(${WITH_GUI})
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.cpp src/apps/crs_selector_ui.cpp src/apps/tops_thresholds_ui.cpp)
	target_link_libraries (treetops-app PUBLIC treetops Qt5::Gui Qt5::Core Qt5::Widgets ui_util)
else()
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.hpp)
	target_link_libraries (treetops-app PUBLIC treetops)
endif()


#add_library (spectral SHARED src/spectral.cpp)
#target_link_libraries (spectral ${GDAL_LIBRARY} ${SQLITE3_LIBRARY} ${SPATIALITE_LIBRARY} raster util geotools db)
#add_executable (spectral-app src/apps/spectral.cpp)
#target_link_libraries (spectral-app PUBLIC spectral)


#add_library(pointstats SHARED src/pointstats.cpp)
#target_link_libraries (pointstats ${GDAL_LIBRARY} laszip las ${Boost_LIBRARIES} raster util geotools cellstats lasreader)
#if(${WITH_GUI})
#	add_executable (pointstats-app src/apps/pointstats.cpp src/apps/pointstats_ui.cpp src/apps/crs_selector_ui.cpp src/apps/filelist.cpp)
#	target_link_libraries (pointstats-app PUBLIC pointstats Qt5::Gui Qt5::Core Qt5::Widgets)
#else()
#	add_executable (pointstats-app src/apps/pointstats.cpp)
#	target_link_libraries (pointstats-app PUBLIC pointstats)
#endif()

#add_library(pointnormalize SHARED src/pointnormalize.cpp)
#target_link_libraries (pointnormalize ${GDAL_LIBRARY} CGAL gmp laszip las ${Boost_LIBRARIES} raster geotools lasreader)
#if(${WITH_GUI})
#	add_executable (pointnormalize-app src/apps/pointnormalize.cpp src/apps/pointnormalize_ui.cpp src/apps/filelist.cpp)
#	target_link_libraries (pointnormalize-app PUBLIC pointnormalize Qt5::Gui Qt5::Core Qt5::Widgets)
#else()
#	add_executable (pointnormalize-app src/apps/pointnormalize.cpp)
#	target_link_libraries (pointnormalize-app PUBLIC pointnormalize)
#endif()

#add_executable (lasboundary src/lasboundary.cpp src/lasutil.cpp)
#target_link_libraries (lasboundary PUBLIC ${GDAL_LIBRARY} lasreader ${GEOS_LIBRARY} ${Boost_LIBRARIES} raster util geotools)


#add_library(lasvalidate SHARED src/lasvalidate.cpp src/lasutil.cpp)
#target_link_libraries (lasvalidate las laszip ${GEOS_LIBRARY} raster util geotools)
#add_executable (lasvalidate-app src/apps/lasvalidate.cpp)
#target_link_libraries (lasvalidate-app PUBLIC lasvalidate)


#add_executable (pointzonalstats src/pointzonalstats.cpp src/lasutil.cpp)
#target_link_libraries (pointzonalstats ${GDAL_LIBRARY} las laszip ${GEOS_LIBRARY} ${Boost_LIBRARIES} raster util geotools)

#add_executable (rasternormalize src/rasternormalize.cpp)
#target_link_libraries (rasternormalize raster)

#add_executable (pointgrid src/pointgrid.cpp)
#target_link_libraries (pointgrid ${GDAL_LIBRARY} ${GEOS_LIBRARY} ${Boost_LIBRARIES} las laszip raster util geotools lasreader)

#add_executable (lasmerge src/lasmerge.cpp src/lasutil.cpp)
#target_link_libraries (lasmerge ${GDAL_LIBRARY} ${Boost_LIBRARIES} las laszip util geotools)


#add_executable (lasreclass src/lasreclass.cpp src/lasutil.cpp)
#target_link_libraries (lasreclass las laszip ${Boost_LIBRARIES} util geotools)


#add_executable (lasclip src/lasclip.cpp src/lasutil.cpp)
#target_link_libraries (lasclip ${GDAL_LIBRARY} las laszip ${GEOS_LIBRARY} ${Boost_LIBRARIES} util geotools)


#add_executable (las2csrs src/las2csrs.cpp src/lasutil.cpp)
#target_link_libraries (las2csrs ${GDAL_LIBRARY} las laszip proj ${Boost_LIBRARIES} util raster geotools)


# Install ###########################################################################################################

#install (TARGETS flood DESTINATION bin)
#install (TARGETS pointgrid DESTINATION bin)

install (TARGETS treetops DESTINATION lib)
install (TARGETS treetops-app DESTINATION bin)
install (FILES include/treetops.hpp DESTINATION include)

#install (TARGETS spectral DESTINATION lib)
#install (TARGETS spectral-app DESTINATION bin)
#install (FILES include/spectra.hpp DESTINATION include)

#install (TARGETS rasterstats DESTINATION bin)

#install (TARGETS mosaic-app DESTINATION bin)

#install (TARGETS pointstats DESTINATION lib COMPONENT libraries)
#install (TARGETS pointstats-app RUNTIME DESTINATION bin COMPONENT applications)
#install (FILES include/pointstats.hpp DESTINATION include COMPONENT headers)

#install (TARGETS lasvalidate DESTINATION lib)
#install (TARGETS lasvalidate-app DESTINATION bin)
#install (FILES include/lasvalidate.hpp DESTINATION include)

#install (TARGETS pointzonalstats DESTINATION bin)

#install (TARGETS lasmerge DESTINATION bin)

#install (TARGETS lasreclass DESTINATION bin)

#install (TARGETS lasclip DESTINATION bin)

#install (TARGETS las2csrs DESTINATION bin)
#install (DIRECTORY share/ DESTINATION share PATTERN "share/*")

