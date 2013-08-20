include_directories( "." )

include_directories( ${ARMADILLO_INCLUDE_DIR} )

file( GLOB_RECURSE SRC_LIST "*.cpp" "." )

add_library( libglm ${SRC_LIST} )

target_link_libraries( libglm )
SET_TARGET_PROPERTIES( libglm PROPERTIES OUTPUT_NAME glm )
