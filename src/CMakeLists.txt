include_directories( "." )

include_directories( ${PLINKIO_INCLUDE_DIR} )
include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )

add_library( common_options common_options.cpp )

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

add_executable( besiq-var besiq_var.cpp )
target_link_libraries( besiq-var common_options libdcdf libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-pairs besiq_pairs.cpp )
target_link_libraries( besiq-pairs libplink libcpp-argparse libbesiq libgzstream ${PLINKIO_LIBRARIES} )

add_executable( besiq-view besiq_view.cpp )
target_link_libraries( besiq-view libcpp-argparse libbesiq )

add_executable( besiq-correct besiq_correct.cpp )
target_link_libraries( besiq-correct libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-imputed besiq_imputed.cpp )
target_link_libraries( besiq-imputed libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-meta besiq_meta.cpp )
target_link_libraries( besiq-meta libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

find_package( OpenMP )

add_executable( besiq besiq.cpp )
target_link_libraries( besiq )

add_library( gene_environment gene_environment.cpp )
set_target_properties( gene_environment PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( gene_environment PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")

add_executable( besiq-env besiq_env.cpp )
set_target_properties( besiq-env PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-env PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-env gene_environment libdcdf libglm libbesiq libplink libcpp-argparse libgzstream ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-lars besiq_lars.cpp )
set_target_properties( besiq-lars PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-lars PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-lars gene_environment libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-predict besiq_predict.cpp )
set_target_properties( besiq-predict PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-predict PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-predict gene_environment libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-predict-lm besiq_predict_lm.cpp )
set_target_properties( besiq-predict-lm PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-predict-lm PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-predict-lm gene_environment libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-extract besiq_extract.cpp )
set_target_properties( besiq-extract PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-extract PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-extract gene_environment libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-gxe besiq_gxe.cpp )
set_target_properties( besiq-gxe PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" )
set_target_properties( besiq-gxe PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
target_link_libraries( besiq-gxe gene_environment libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

add_executable( besiq-mglm besiq_mglm.cpp )
target_link_libraries( besiq-mglm libdcdf libglm libbesiq libplink libcpp-argparse ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )

INSTALL( TARGETS besiq besiq-stagewise besiq-bayes besiq-caseonly
    besiq-glm besiq-scaleinv besiq-loglinear besiq-wald besiq-env
    besiq-pairs besiq-view besiq-correct besiq-imputed besiq-var
    besiq-separate besiq-lars besiq-meta besiq-mglm besiq-predict
    besiq-gxe besiq-predict-lm besiq-extract DESTINATION bin )

