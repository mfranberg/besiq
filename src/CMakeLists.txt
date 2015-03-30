include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_executable( bayesic bayesic.cpp )
target_link_libraries( bayesic libdcdf libglm libbayesic libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-env bayesic_env.cpp )
target_link_libraries( bayesic-env libdcdf libglm libbayesic libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-pairs bayesic_pairs.cpp )
target_link_libraries( bayesic-pairs libplink libcpp-argparse libbayesic libgzstream ${PLINKIO_LIBRARIES} )

add_executable( bayesic-view bayesic_view.cpp )
target_link_libraries( bayesic-view libcpp-argparse libbayesic )

add_executable( bayesic-correct bayesic_correct.cpp )
target_link_libraries( bayesic-correct libdcdf libglm libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )
