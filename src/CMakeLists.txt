include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_executable( bayesic bayesic.cpp )
target_link_libraries( bayesic )

add_executable( bayesic-stagewise bayesic_stagewise.cpp )
target_link_libraries( bayesic-stagewise libdcdf libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-bayes bayesic_bayes.cpp )
target_link_libraries( bayesic-bayes libdcdf libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-caseonly bayesic_caseonly.cpp )
target_link_libraries( bayesic-caseonly libdcdf libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-glm bayesic_glm.cpp )
target_link_libraries( bayesic-glm libdcdf libbayesic libglm libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-loglinear bayesic_loglinear.cpp )
target_link_libraries( bayesic-loglinear libdcdf libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-wald bayesic_wald.cpp )
target_link_libraries( bayesic-wald libdcdf libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-env bayesic_env.cpp )
target_link_libraries( bayesic-env libdcdf libglm libbayesic libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( bayesic-pairs bayesic_pairs.cpp )
target_link_libraries( bayesic-pairs libplink libcpp-argparse libbayesic libgzstream ${PLINKIO_LIBRARIES} )

add_executable( bayesic-view bayesic_view.cpp )
target_link_libraries( bayesic-view libcpp-argparse libbayesic )

add_executable( bayesic-correct bayesic_correct.cpp )
target_link_libraries( bayesic-correct libdcdf libglm libbayesic libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

