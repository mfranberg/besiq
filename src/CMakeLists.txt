include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_library( common_options common_options.cpp )

add_executable( besiq besiq.cpp )
target_link_libraries( besiq )

add_executable( besiq-stagewise besiq_stagewise.cpp )
target_link_libraries( besiq-stagewise common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-bayes besiq_bayes.cpp )
target_link_libraries( besiq-bayes common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-caseonly besiq_caseonly.cpp )
target_link_libraries( besiq-caseonly common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-glm besiq_glm.cpp )
target_link_libraries( besiq-glm common_options libdcdf libbesiq libglm libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-scaleinv besiq_scaleinv.cpp )
target_link_libraries( besiq-scaleinv common_options libdcdf libbesiq libglm libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-separate besiq_separate.cpp )
target_link_libraries( besiq-separate common_options libdcdf libbesiq libglm libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-loglinear besiq_loglinear.cpp )
target_link_libraries( besiq-loglinear common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-wald besiq_wald.cpp )
target_link_libraries( besiq-wald common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-env besiq_env.cpp )
target_link_libraries( besiq-env libdcdf libglm libbesiq libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-pairs besiq_pairs.cpp )
target_link_libraries( besiq-pairs libplink libcpp-argparse libbesiq libgzstream ${PLINKIO_LIBRARIES} )

add_executable( besiq-view besiq_view.cpp )
target_link_libraries( besiq-view libcpp-argparse libbesiq )

add_executable( besiq-correct besiq_correct.cpp )
target_link_libraries( besiq-correct libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-imputed besiq_imputed.cpp )
target_link_libraries( besiq-imputed libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

INSTALL( TARGETS besiq besiq-stagewise besiq-bayes besiq-caseonly
    besiq-glm besiq-scaleinv besiq-loglinear besiq-wald besiq-env
    besiq-pairs besiq-view besiq-correct besiq-imputed DESTINATION bin )

