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
#set (CMAKE_POSITION_INDEPENDENT_CODE ON) # Qt complains if this is set and not -fPIC below...
set (CMAKE_CXX_FLAGS "-fPIC -lpthread -fopenmp -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/opt/lib;/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/opt/lib;/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
#set (CMAKE_VERBOSE_MAKEFILE ON)

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
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:\\Qt\\5.9.1\\msvc2015_64")
	set (LIB_PATH "C:\\dev\\geos\\makefiles\\lib\\Release;C:\\dev\\gisinternals\\lib")
	set (INC_PATH "C:\\dev\\geos\\include;C:\\dev\\gisinternals\\include;C:\\dev\\eigen3;C:\\dev\\treetops\\makefiles\\ui_util_autogen\\include")
	set (GDAL_INCLUDE_DIR "C:\\dev\\gisinternals\\include")
	set (GDAL_LIBRARY "C:\\dev\\gisinternals\\lib\\gdal_i.lib")
	set (libLAS_INCLUDE_DIR "C:\\dev\\gisinternals\\include")
	set (libLAS_LIBRARY_DIRS "C:\\dev\\gisinternals\\lib")
else()
endif()

if(${WITH_GUI})
	set (CMAKE_AUTOUIC ON)
	set (CMAKE_AUTOMOC ON)
	set (CMAKE_INCLUDE_CURRENT_DIR ON)
	set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)
find_package (OpenMP REQUIRED)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

if(${WITH_GUI})
	find_package(Qt5Widgets REQUIRED)
endif()

project (treetops)
add_subdirectory (libtreetops)
add_subdirectory (libgeo)

# Configure directories

include_directories(BEFORE SYSTEM include build/ui_util_autogen/include libgeo/include 
	libtreetops/include)
	
link_directories(build/lib /opt/lib)

# Build support libraries ############################################################################################

if(${WITH_GUI})
	add_library (ui_util STATIC src/ui_util.cpp src/crs_selector_ui.cpp src/tops_thresholds_ui.cpp src/crowns_thresholds_ui.cpp)
	target_link_libraries(ui_util Qt5::Widgets georaster treetops)
	target_include_directories(ui_util PUBLIC include ui_util_autogen/include  
		${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
endif()


# Build programs #####################################################################################################

if(${WITH_GUI})
	add_executable (treetops-app src/treetops.cpp src/treetops_ui.cpp src/settings.cpp)
	target_include_directories(treetops-app PUBLIC libgeo/include include ${INC_PATH} 
		${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
	target_link_libraries (treetops-app PUBLIC Qt5::Gui Qt5::Core Qt5::Widgets geoutil georaster treetops ui_util)
else()
	add_executable (treetops-app src/treetops.cpp src/treetops_ui.hpp)
	target_include_directories(treetops-app PUBLIC libgeo/include include ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS})
	target_link_libraries (treetops-app PUBLIC treetops geoutil georaster treetops)
endif()

set_target_properties (treetops-app PROPERTIES OUTPUT_NAME treetops)

# Install ###########################################################################################################

install (TARGETS treetops-app DESTINATION bin)
