include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBPLINK_INCLUDE_DIR} )
include_directories( ${LIBBAYESIC_INCLUDE_DIR} )
include_directories( ${LIBCPPARGPARSE_INCLUDE_DIR} )
include_directories( ${LIBGLM_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )
include_directories( ${GZSTREAM_INCLUDE_DIR} )
include_directories( ${SHAREDPTR_INCLUDE_DIR} )

add_executable( bayesic bayesic.cpp )
target_link_libraries( bayesic libdcdf libglm libbayesic libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-pairs bayesic_pairs.cpp )
target_link_libraries( bayesic-pairs libplink libcpp-argparse libbayesic libgzstream ${PLINKIO_LIBRARIES} )

