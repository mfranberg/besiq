include_directories( "." )

add_library( libgzstream gzstream.cpp )
target_link_libraries( libgzstream -lz )

SET_TARGET_PROPERTIES( libgzstream PROPERTIES OUTPUT_NAME gzstream )
