include_directories( "." )

include_directories( ${ARMADILLO_INCLUDE_DIR} )
include_directories( ${LIBPLINK_INCLUDE_DIR} )
include_directories( ${LIBGLM_INCLUDE_DIR} )
include_directories( ${DCDFLIB_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libbayesic ${SRC_LIST} )

target_link_libraries( libbayesic )
SET_TARGET_PROPERTIES( libbayesic PROPERTIES OUTPUT_NAME bayesic )
