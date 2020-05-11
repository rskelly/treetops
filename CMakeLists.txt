cmake_minimum_required(VERSION 3.10)

project (treetops)
add_subdirectory (libtreetops)
add_subdirectory (libgeo)

set (BUILD_APPS false CACHE BOOL "If true, build apps in libgeo.")
set (WITH_GUI true CACHE BOOL "If true, will build the Qt GUI.")

if(${WITH_GUI})
	add_definitions(-DWITH_GUI)
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS "-Wall -D_GLIBCXX_PARALLEL")
set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(${WITH_GUI})
	set (CMAKE_AUTOUIC ON)
	set (CMAKE_AUTOMOC ON)
	set (CMAKE_INCLUDE_CURRENT_DIR ON)
	set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")
endif()

set (SQLITE_LIBRARY sqlite3)
set (SPATIALITE_LIBRARY spatialite)

if (${APPLE})
# Applie configs probably not good. Used for specific installation.
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
	set (DEV_DIR "C:\\dev")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\OSGeo4W64")
	set (GEOS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin") #set (GEOS_LIBRARYDIR "${DEV_DIR}\\geos\\lib;${DEV_DIR}\\geos\\bin")
	set (GEOS_INCLUDEDIR "${GIS_DIR}\\include") #set (GEOS_INCLUDEDIR "${DEV_DIR}\\geos\\include")
	set (EIGEN3_INCLUDE_DIR "${GIS_DIR}\\include\\eigen3")
	set (GDAL_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (GIS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GIS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GIS_DIR}\\bin;${GIS_DIR}\\lib;${DEV_DIR}\\Qt\\5.11.0\\msvc2017_64\\lib\\cmake")
	set (SQLITE_LIBRARY sqlite3_i)
	set (SPATIALITE_LIBRARY spatialite_i)
	set (GDAL_LIBRARY gdal_i)
else()
	if(${WITH_GUI})
		set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}:/usr/lib/x86_64-linux-gnu/cmake/Qt5")
	endif()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unknown-pragmas")
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (GDAL 2 REQUIRED)
find_package (OpenMP REQUIRED)

if(${WITH_GUI})
	find_package(Qt5 REQUIRED Widgets)
endif()

find_library (GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

# Get the git rev and create an environment variable to be used in compilation.

execute_process (COMMAND git rev-parse --short HEAD RESULT_VARIABLE GIT_RESULT 
	OUTPUT_VARIABLE GIT_OUTPUT ERROR_VARIABLE GIT_ERROR)
string (STRIP "${GIT_OUTPUT}" GIT_REV)
message (STATUS "Finding Git rev: ${GIT_REV} (Result: ${GIT_RESULT}, Error: ${GIT_ERROR})")

# Configure directories

include_directories(BEFORE SYSTEM
	libgeo/include 
	libtreetops/include 
	${EIGEN_INCLUDEDIR}	
	${GDAL_INCLUDE_DIR}
	${Qt5Core_INCLUDE_DIRS} 
	${Qt5Widgets_INCLUDE_DIRS} 
	${Qt5Gui_INCLUDE_DIRS}
	${GEOS_INCLUDEDIR} 
)
	 
link_directories(
	${GEOS_LIBRARYDIR} 
	${GIS_LIBRARYDIR} 
)
	
# Build support libraries ############################################################################################

set (ui_utilsrc src/ui_util.cpp src/tops_thresholds_ui.cpp src/crowns_thresholds_ui.cpp)
add_library (ui_utilobj OBJECT ${ui_utilsrc})
target_compile_features(ui_utilobj PUBLIC cxx_std_17)
set_property (TARGET ui_utilobj PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library (ui_util SHARED $<TARGET_OBJECTS:ui_utilobj>)
target_link_libraries(ui_util treetops 
	Qt5::Gui Qt5::Core Qt5::Widgets 
	${SQLITE_LIBRARY} 
	${SPATIALITE_LIBRARY} 
)

# Build programs #####################################################################################################

add_executable (treetops-app WIN32 src/treetops.cpp src/treetops_ui.cpp src/settings.cpp)
target_compile_features(treetops-app PUBLIC cxx_std_17)
target_compile_definitions (treetops-app PUBLIC GIT_REV=${GIT_REV})
target_link_libraries (treetops-app PUBLIC ui_util treetops 
	Qt5::Gui Qt5::Core Qt5::Widgets 
	${GEOS_LIBRARY}
)

# Install ###########################################################################################################

install (TARGETS treetops-app RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
