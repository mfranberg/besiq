include_directories( "." )

include_directories( ${LIBBAYESIC_INCLUDE_DIR} )
include_directories( ${LIBCPPARGPARSE_INCLUDE_DIR} )
include_directories( ${LIBGLM_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_executable( bayesic bayesic.cpp )
target_link_libraries( bayesic libglm libbayesic libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} )
