include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBPLINK_INCLUDE_DIR} )
include_directories( ${LIBBAYESIC_INCLUDE_DIR} )
include_directories( ${LIBCPPARGPARSE_INCLUDE_DIR} )
include_directories( ${LIBGLM_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_executable( bayesic bayesic.cpp )
target_link_libraries( bayesic libdcdf libglm libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( generate_pairs generate_pairs.cpp )
target_link_libraries( generate_pairs libplink libcpp-argparse ${PLINKIO_LIBRARIES} )

