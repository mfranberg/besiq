include_directories( ${LIBS_INCLUDE_DIR} )
include_directories( ${GTEST_INCLUDE_DIR} )
include_directories( ${ARMADILLO_INCLUDE_DIR} )
include_directories( ${DCDFLIB_INCLUDE_DIR} )

file( GLOB TEST_LIST "*_test.cpp" )

foreach( TEST_PATH ${TEST_LIST} )
    get_filename_component( TEST_NAME ${TEST_PATH} NAME_WE )
    add_executable( ${TEST_NAME} ${TEST_PATH} )
    target_link_libraries( ${TEST_NAME} gtest gtest_main libglm libbesiq
        libplink libdcdf ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${PLINKIO_LIBRARIES} )
    add_test( ${TEST_NAME} ${TEST_NAME} )
endforeach( TEST_PATH )
