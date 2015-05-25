include_directories( ${ARMADILLO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${PLINKIO_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libbesiq ${SRC_LIST} )

target_link_libraries( libbesiq libglm )
SET_TARGET_PROPERTIES( libbesiq PROPERTIES OUTPUT_NAME besiq )
