include_directories( ${LIBS_INCLUDE_DIR} )

file( GLOB_RECURSE C_SRC_LIST "*.c" "." )
file( GLOB_RECURSE CPP_SRC_LIST "*.cpp" "." )

add_library( libdcdf ${C_SRC_LIST} ${CPP_SRC_LIST} )

target_link_libraries( libdcdf )
SET_TARGET_PROPERTIES( libdcdf PROPERTIES OUTPUT_NAME dcdf )
