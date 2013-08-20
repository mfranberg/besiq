project( bayesic )
cmake_minimum_required( VERSION 2.7 )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/modules/" )

set( BAYESIC_VERSION_MAJOR 0 )
set( BAYESIC_VERSION_MINOR 1 )
set( BAYESIC_VERSION_PATCH 0 )

find_package( Armadillo REQUIRED )
find_package( BLAS REQUIRED )

set( LIBBAYESIC_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/bayesic/" )
set( LIBCPPARGPARSE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/cpp-argparse/" )
set( LIBGLM_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/glm/" )
set( GTEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/gtest/include/" )
#set( LIBCPPTHREAD_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/libcppthread/" )

enable_testing( )

add_subdirectory( libs )
add_subdirectory( src )
add_subdirectory( tests )
