include_directories( "." )

include_directories( ${ARMADILLO_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libcpp-argparse ${SRC_LIST} )

target_link_libraries( libcpp-argparse )
SET_TARGET_PROPERTIES( libcpp-argparse PROPERTIES OUTPUT_NAME cpp-argparse )
