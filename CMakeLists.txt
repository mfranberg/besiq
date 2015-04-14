project( bayesic )
cmake_minimum_required( VERSION 2.8.9 )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/modules/" )

set( BAYESIC_VERSION_MAJOR 0 )
set( BAYESIC_VERSION_MINOR 1 )
set( BAYESIC_VERSION_PATCH 0 )

if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
    # Silently ignore Armadillo template warning
    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -wd488" )
    set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -wd488" )
    set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} -wd488" )
endif( )

find_package( Armadillo REQUIRED )
find_package( BLAS REQUIRED )

include( CheckIncludeFiles )
check_include_files( "tr1/random", HAVE_TR1_RANDOM )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/libs/bayesic/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/libs/bayesic/config.h )
include_directories( ${CMAKE_BINARY_DIR}/libs/ )
if( NOT HAVE_TR1_RANDOM )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif( )

set( LIBS_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/" )

set( PLINKIO_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/libplinkio/src/" )
set( GTEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/gtest/include/" )
set( PLINKIO_LIBRARIES libplinkio-static )

add_definitions( -DGZSTREAM_NAMESPACE=gz )

enable_testing( )

set( DISABLE_TESTS ON CACHE BOOL "Build tests in libplinkio" )

add_subdirectory( libs )
add_subdirectory( src )
add_subdirectory( tests )
