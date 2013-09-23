include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
message( ${PLINKIO_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libplink ${SRC_LIST} )

target_link_libraries( libplink )
SET_TARGET_PROPERTIES( libplink PROPERTIES OUTPUT_NAME plink )
