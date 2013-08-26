include_directories( "." )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libplink ${SRC_LIST} )

target_link_libraries( libplink )
SET_TARGET_PROPERTIES( libplink PROPERTIES OUTPUT_NAME plink )
