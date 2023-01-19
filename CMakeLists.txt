# CMakeFiles.txt for the MRISIM library
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
SET(MRISIM_PACKAGE_VERSION_MAJOR 1)
SET(MRISIM_PACKAGE_VERSION_MINOR 0)
SET(MRISIM_PACKAGE_VERSION_PATCH 20)


# Packaging defines
IF(NOT MINC_TOOLKIT_BUILD)
  PROJECT(mrisim)
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${MRISIM_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${MRISIM_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${MRISIM_PACKAGE_VERSION_PATCH})
  OPTION(MRISIM_BUILD_SHARED_LIBS "Build shared libraries." OFF)
  INCLUDE(CPack)
  INCLUDE(CTest)
ENDIF(NOT MINC_TOOLKIT_BUILD)

SET(PACKAGE "MRISIM")
SET(PACKAGE_BUGREPORT "a.janke@gmail.com")

SET(PACKAGE_NAME "MRISIM")
SET(PACKAGE_VERSION "${MRISIM_PACKAGE_VERSION_MAJOR}.${MRISIM_PACKAGE_VERSION_MINOR}.${MRISIM_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

IF(NOT MINC_TOOLKIT_BUILD)
  FIND_PACKAGE( LIBMINC REQUIRED )
ENDIF(NOT MINC_TOOLKIT_BUILD)

INCLUDE( ${LIBMINC_USE_FILE} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
IF(NOT APPLE)
ADD_DEFINITIONS(-D_XOPEN_SOURCE=500)
ENDIF(NOT APPLE)
# aliases
SET(VERSION "${PACKAGE_VERSION}")

# check for prereqs
INCLUDE(CheckIncludeFiles)

CHECK_INCLUDE_FILES(dlfcn.h     HAVE_DLFCN_H)
CHECK_INCLUDE_FILES(sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(stdint.h    HAVE_STDINT_H)
CHECK_INCLUDE_FILES(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(dirent.h    HAVE_DIRENT_H)
CHECK_INCLUDE_FILES(memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES(strings.h   HAVE_STRINGS_H)

# config files for build
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)
# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in ${CMAKE_CURRENT_BINARY_DIR}/epm-header)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} src )


INCLUDE( ${LIBMINC_USE_FILE} )
LINK_LIBRARIES(${LIBMINC_LIBRARIES})


SET(coil_files 
  coils/coil_pd_pn2_1mm.rf 
  coils/coil.rf 
  coils/coil_t1_intrinsic.rf 
  coils/coil_t1_pn2_1mm.rf 
  coils/coil_t2_pn2_1mm.rf )

SET(doc_files 
  doc/mrisim.man 
  doc/mrisim.ps 
  doc/switch.txt)

SET(field_files 
  fields/Philips/makerf.zip 
  fields/Philips/qbc_id.mnc.Z 
  fields/Philips/qbc_q450.mnc.Z 
  fields/Philips/qbc_q86.mnc.Z 
  fields/Philips/qhc.mnc.Z )

SET(sequence_files 
  sequences/pd_ai.seq 
  sequences/pd_icbm.seq 
  sequences/t1_ai.seq 
  sequences/t1_icbm_avg.seq 
  sequences/t1_icbm.seq 
  sequences/t2_ai.seq 
  sequences/t2_icbm_avg.seq
  sequences/t2_icbm.seq)


SET(MRISIM_HEADERS
  src/minc/chirp.h 
  src/minc/fourn.h 
  src/minc/imincfile.h 
  src/minc/iomincfile.h 
  src/minc/mincfile.h 
  src/minc/mincicv.h 
  src/minc/mriimage.h 
  src/minc/mrilabel.h 
  src/minc/mrimatrix.h 
  src/minc/mriminc.h 
  src/minc/mristring.h 
  src/minc/mrivolume.h 
  src/minc/omincfile.h 
  src/minc/time_stamp.h 
  src/mrisim/discrete_label_phantom.h 
  src/mrisim/discrete_phantom.h 
  src/mrisim/discrete_rf_phantom.h 
  src/mrisim/fuzzy_label_phantom.h 
  src/mrisim/fuzzy_phantom.h 
  src/mrisim/fuzzy_rf_phantom.h 
  src/mrisim/image_snr_coil.h 
  src/mrisim/intrinsic_coil.h 
  src/mrisim/mriscanner.h 
  src/mrisim/mrisimargs.h 
  src/mrisim/mrisim.h 
  src/mrisim/mrisim_main.h 
  src/mrisim/paramfile.h 
#  src/mrisim/ParseArgv.h 
  src/mrisim/percent_coil.h 
  src/mrisim/phantom.h 
  src/mrisim/rf_coil.h 
  src/mrisim/rf_tissue_phantom.h 
  src/mrisim/scanner_output.h 
  src/mrisim/tissue_phantom.h 
  src/signal/ce_fast.h 
  src/signal/customseq.h 
  src/signal/event.h 
  src/signal/fast_iso_model.h 
  src/signal/ffe.h 
  src/signal/fisp.h 
  src/signal/flash.h 
  src/signal/ir.h 
  src/signal/isochromat_model.h 
  src/signal/pulseseq.h 
  src/signal/quick_model.h 
  src/signal/quickseq.h 
  src/signal/repeat.h 
  src/signal/rf_pulse.h 
  src/signal/sample.h 
  src/signal/se.h 
  src/signal/signal.h 
  src/signal/spin_model.h 
  src/signal/spoiled_flash.h 
  src/signal/spoiler.h 
  src/signal/tissue.h 
  src/signal/vector.h 
  src/signal/vector_model.h
)

# some variables
SET(MRISIM_SRCS
  ${MRISIM_HEADERS}
  src/minc/chirp.cxx 
  src/minc/fourn.c 
  src/minc/imincfile.cxx 
  src/minc/iomincfile.cxx 
  src/minc/mincfile.cxx 
  src/minc/mincicv.cxx 
  src/minc/mriimage.cxx 
  src/minc/mrilabel.cxx 
  src/minc/mrimatrix.cxx 
  src/minc/mristring.cxx 
  src/minc/mrivolume.cxx 
  src/minc/omincfile.cxx 
  src/minc/time_stamp.c 
  src/mrisim/discrete_label_phantom.cxx 
  src/mrisim/discrete_phantom.cxx 
  src/mrisim/discrete_rf_phantom.cxx 
  src/mrisim/fuzzy_label_phantom.cxx 
  src/mrisim/fuzzy_phantom.cxx 
  src/mrisim/fuzzy_rf_phantom.cxx 
  src/mrisim/image_snr_coil.cxx 
  src/mrisim/intrinsic_coil.cxx 
  src/mrisim/mriscanner.cxx 
  src/mrisim/mrisimargs.cxx 
  src/mrisim/mrisim_main.cxx 
  src/mrisim/paramfile.cxx 
#  src/mrisim/ParseArgv.c 
  src/mrisim/percent_coil.cxx 
  src/mrisim/phantom.cxx 
  src/mrisim/rf_coil.cxx 
  src/mrisim/rf_tissue_phantom.cxx 
  src/mrisim/scanner_output.cxx 
  src/mrisim/tissue_phantom.cxx 
  src/signal/ce_fast.cxx 
  src/signal/customseq.cxx 
  src/signal/event.cxx 
  src/signal/fast_iso_model.cxx 
  src/signal/ffe.cxx 
  src/signal/fisp.cxx 
  src/signal/flash.cxx 
  src/signal/ir.cxx 
  src/signal/isochromat_model.cxx 
  src/signal/pulseseq.cxx 
  src/signal/quick_model.cxx 
  src/signal/quickseq.cxx 
  src/signal/repeat.cxx 
  src/signal/rf_pulse.cxx 
  src/signal/sample.cxx 
  src/signal/se.cxx 
  src/signal/spin_model.cxx 
  src/signal/spoiled_flash.cxx 
  src/signal/spoiler.cxx 
  src/signal/tissue.cxx 
  src/signal/vector.cxx 
  src/signal/vector_model.cxx
)

ADD_EXECUTABLE(mrisim ${MRISIM_SRCS})

INSTALL(TARGETS 
  mrisim 
  DESTINATION bin)

# testing
IF(BUILD_TESTING)
   ADD_SUBDIRECTORY( testdir )
ENDIF(BUILD_TESTING)
