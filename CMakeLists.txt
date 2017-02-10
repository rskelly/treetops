cmake_minimum_required(VERSION 3.1)

set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (WITH_GUI true CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -D_GLIBCXX_PARALLEL")

if (${APPLE})
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH "/usr/local/Cellar/qt5/5.6.1-1/lib/cmake/Qt5Widgets")
	endif()
	set (GDAL_LIBRARY /usr/local/opt/gdal2/lib/libgdal.dylib)
	set (GDAL_INCLUDE_DIR /usr/local/opt/gdal2/include)
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
elseif(${WIN32})
	set (Boost_COMPILER -vc140)
	set (Boost_NO_SYSTEM_PATHS  ON)
	set (BOOST_INCLUDEDIR "E:\\boost_1_63_0")
	set (BOOST_LIBRARYDIR "E:\\boost_1_63_0\\lib64-msvc-14.0")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:\\Qt\\5.7\\msvc2015_64\\lib\\cmake\\Qt5Widgets")
	set (LIB_PATH "E:\\opt\\lib")
	set (INC_PATH "E:\\opt\\include")
	set (GDAL_INCLUDE_DIR "E:\\opt\\include")
	set (GDAL_LIBRARY "E:\\opt\\lib\\gdal_i.lib")
	set (libLAS_INCLUDE_DIR "E:\\opt\\include")
	set (libLAS_LIBRARY_DIRS "E:\\opt\\lib")
else()
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:/usr/lib/x86_64-linux-gnu/cmake/Qt5")
	endif()
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

if(${WITH_GUI})
	set (CMAKE_AUTOUIC ON)
	set (CMAKE_AUTOMOC ON)
	set (CMAKE_INCLUDE_CURRENT_DIR ON)
	set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc/apps")
endif()

if(${WITH_GUI})
	find_package(Qt5Widgets REQUIRED)
endif()

project (treetops)
add_subdirectory (libtreetops)

# Configure directories
link_directories(. ${Boost_LIBRARY_DIRS} ${LIB_PATH})
include_directories(libtreetops/libgeo/include libtreetops/include src/apps include 
	${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${INC_PATH})


# Build support libraries ############################################################################################

if(${WITH_GUI})
	add_library (ui_util STATIC src/apps/ui_util.cpp src/apps/crs_selector_ui.cpp src/apps/tops_thresholds_ui.cpp)
	target_link_libraries(ui_util Qt5::Widgets)
endif()


# Build programs #####################################################################################################

if(${WITH_GUI})
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.cpp )
	target_link_libraries (treetops-app PUBLIC Qt5::Gui Qt5::Core Qt5::Widgets treetops ui_util)
else()
	add_executable (treetops-app src/apps/treetops.cpp src/apps/treetops_ui.hpp)
	target_link_libraries (treetops-app PUBLIC treetops)
endif()


# Install ###########################################################################################################

install (TARGETS treetops-app DESTINATION bin)
