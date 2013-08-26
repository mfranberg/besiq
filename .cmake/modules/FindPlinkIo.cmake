# - Try to find PlinkIo
# Once done this will define
#  PLINKIO_FOUND - System has libplinkio
#  PLINKIO_INCLUDE_DIRS - The libplinkio include directories
#  PLINKIO_LIBRARIES - The libraries needed to use libplinkio
#  PLINKIO_DEFINITIONS - Compiler switches required for using libplinkio

find_package( PkgConfig )
pkg_check_modules( PC_PLINKIO QUIET plinkio )
set( PLINKIO_DEFINITIONS ${PC_PLINKIO_CFLAGS_OTHER} )

find_path( PLINKIO_INCLUDE_DIR plinkio/plinkio.h
           HINTS ${PC_PLINKIO_INCLUDEDIR} ${PC_PLINKIO_INCLUDE_DIRS}
           PATH_SUFFIXES plinkio )

find_library( PLINKIO_LIBRARY NAMES plinkio
              HINTS ${PC_PLINKIO_LIBDIR} ${PC_PLINKIO_LIBRARY_DIRS} )

set( PLINKIO_LIBRARIES ${PLINKIO_LIBRARY} )
set( PLINKIO_INCLUDE_DIRS ${PLINKIO_INCLUDE_DIR} )

include( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set LIBPLINKIO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( PlinkIo DEFAULT_MSG
                                   PLINKIO_LIBRARY PLINKIO_INCLUDE_DIR)

mark_as_advanced( PLINKIO_INCLUDE_DIR PLINKIO_LIBRARY )
