include_directories( ${LIBS_INCLUDE_DIR} )

add_library( libgzstream gzstream.cpp gzutil.cpp )
target_link_libraries( libgzstream -lz )

SET_TARGET_PROPERTIES( libgzstream PROPERTIES OUTPUT_NAME gzstream )
