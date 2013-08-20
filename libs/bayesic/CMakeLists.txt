include_directories( "." )

include_directories( ${ARMADILLO_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libbayesic ${SRC_LIST} )

target_link_libraries( libbayesic )
SET_TARGET_PROPERTIES( libbayesic PROPERTIES OUTPUT_NAME bayesic )
