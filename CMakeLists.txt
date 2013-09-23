project( bayesic )
cmake_minimum_required( VERSION 2.8.9 )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/modules/" )

set( BAYESIC_VERSION_MAJOR 0 )
set( BAYESIC_VERSION_MINOR 1 )
set( BAYESIC_VERSION_PATCH 0 )

find_package( Armadillo REQUIRED )
find_package( BLAS REQUIRED )

set( LIBBAYESIC_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/bayesic/" )
set( LIBPLINK_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/plink/" )
set( LIBCPPARGPARSE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/cpp-argparse/" )
set( LIBGLM_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/glm/" )
set( GTEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/gtest/include/" )
set( DCDFLIB_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/dcdflib/" )
set( GZSTREAM_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/gzstream/" )
set( SHAREDPTR_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/shared_ptr/" )

set( PLINKIO_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/libplinkio/src/" )
set( PLINKIO_LIBRARIES libplinkio-static )

add_definitions( -DGZSTREAM_NAMESPACE=gz )

enable_testing( )

set( DISABLE_TESTS ON CACHE BOOL "Build tests in libplinkio" )

add_subdirectory( libs )
add_subdirectory( src )
add_subdirectory( tests )
