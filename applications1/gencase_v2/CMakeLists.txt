if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
   message(FATAL_ERROR "femus gencase cannot (yet) be compiled standalone, compile it from the app directory.")
endif()

# dependencies

if(DEFINED ENV{PLAT_DIR})
  set(MPI_HOME "$ENV{PLAT_THIRD_PARTY_DIR}/openmpi")
endif()
find_package(MPI
  REQUIRED
  COMPONENTS
    CXX
)

find_package(hdf5
  REQUIRED
  NO_MODULE
  PATHS
    $ENV{PLAT_THIRD_PARTY_DIR}/hdf5/share/cmake/hdf5
)

find_package(PkgConfig REQUIRED)
set(ENV{PKG_CONFIG_PATH} $ENV{PLAT_CODES_DIR}/libmesh/lib/pkgconfig)
pkg_check_modules(libmesh REQUIRED IMPORTED_TARGET libmesh)
# workaround for bogus libmesh include path
get_target_property(
  LIBMESH_INCLUDE_DIRECTORIES PkgConfig::libmesh INTERFACE_INCLUDE_DIRECTORIES)
list(GET LIBMESH_INCLUDE_DIRECTORIES 0 LIBMESH_INCLUDE_ROOT)
list(REMOVE_ITEM LIBMESH_INCLUDE_DIRECTORIES /usr/include/eigen3)
list(REMOVE_ITEM LIBMESH_INCLUDE_DIRECTORIES /usr/include/tirpc)
list(REMOVE_ITEM LIBMESH_INCLUDE_DIRECTORIES $ENV{PLAT_THIRD_PARTY_DIR}/petsc/linux-opt/include)
list(REMOVE_ITEM LIBMESH_INCLUDE_DIRECTORIES ${LIBMESH_INCLUDE_ROOT}/contrib)
set_property(TARGET PkgConfig::libmesh
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${LIBMESH_INCLUDE_DIRECTORIES})
  
# # get_target_property(tmp PkgConfig::libmesh INTERFACE_COMPILE_DEFINITIONS)
get_target_property(COMPILE_DEFINITIONS PkgConfig::libmesh INTERFACE_COMPILE_DEFINITIONS_DEBUG)
message(STATUS "libmesh definitions: ${COMPILE_DEFINITIONS}")

find_package(Eigen3 REQUIRED)

# executable

add_executable(femus_gencase)

set(GENCASE_HEADERS
  ${FEMUS_GENCASE_DIR}/src/gambit_io.h
  ${FEMUS_GENCASE_DIR}/src/gencase_main.h
  ${FEMUS_GENCASE_DIR}/src/Libmeshlib_conf.h
  ${FEMUS_GENCASE_DIR}/src/MGGenCase_conf.h
  ${FEMUS_GENCASE_DIR}/src/MGGenCase.h
  # ${FEMUS_GENCASE_DIR}/src/MGFE_conf.h
  ${FEMUS_GENCASE_DIR}/src/MGGeomEl.h
  ${FEMUS_GENCASE_DIR}/src/MGUtils.h
  ${FEMUS_GENCASE_DIR}/src/namebased_io.h
  ${FEMUS_GENCASE_DIR}/src/salome_io.h
)

set(GENCASE_SOURCES
  ${FEMUS_GENCASE_DIR}/src/gambit_io.C
  ${FEMUS_GENCASE_DIR}/src/gencase_main.C
  ${FEMUS_GENCASE_DIR}/src/MGGenCase.C
  ${FEMUS_GENCASE_DIR}/src/MGGenCase_OP_Matrix.C
  ${FEMUS_GENCASE_DIR}/src/MGGenCase_OP_PROL.C
  ${FEMUS_GENCASE_DIR}/src/MGGenCase_OP_REST.C
  ${FEMUS_GENCASE_DIR}/src/MGGenCase_PRINT.C
  ${FEMUS_GENCASE_DIR}/src/MGGeomEl.C
  ${FEMUS_GENCASE_DIR}/src/MGUtils.C
  ${FEMUS_GENCASE_DIR}/src/namebased_io.C
  ${FEMUS_GENCASE_DIR}/src/salome_io.C
)

target_sources(femus_gencase
  PUBLIC
    ${GENCASE_HEADERS}
  PRIVATE
    ${GENCASE_SOURCES}
)

target_include_directories(femus_gencase
  PUBLIC
    ${FEMUS_GENCASE_DIR}/src
    ${APP_DIR}/DATA
    ${CMAKE_BINARY_DIR}
)

target_link_libraries(femus_gencase
  PUBLIC
    MPI::MPI_CXX
    # hdf5::hdf5-shared
    medcoupling::medcoupling
    # PkgConfig::libmesh
    $ENV{PLAT_CODES_DIR}/libmesh/lib/libmesh_opt.so
    $ENV{PLAT_CODES_DIR}/libmesh/lib/libtimpi_opt.so
    $ENV{PLAT_THIRD_PARTY_DIR}/petsc/lib/libpetsc.so
    Eigen3::Eigen
)
target_include_directories(femus_gencase
  PUBLIC
    $ENV{PLAT_CODES_DIR}/libmesh/include
    $ENV{APP_PATH}/DATA/
    $ENV{PLAT_THIRD_PARTY_DIR}/petsc/include
)

set_target_properties(femus_gencase PROPERTIES LINKER_LANGUAGE CXX)

# add_custom_command(
#   OUTPUT
#     ../RESU/msh1.Matrix0.h5
#   COMMAND
#     ${CMAKE_COMMAND} -E env
#       FEMUS_DIR=${FEMUS_DIR}
#       APP_DIR=${CMAKE_SOURCE_DIR}
#       FM_MYAPP=${CMAKE_SOURCE_DIR}
#     mpiexec -np $ENV{NP} ${CMAKE_BINARY_DIR}/femus_gencase
#   DEPENDS
#     ${CMAKE_BINARY_DIR}/femus_gencase
#   COMMENT
#     "run gencase"
# )
# add_custom_target(run_gencase
#   DEPENDS
#     ../RESU/msh1.Matrix0.h5
# )
# add_dependencies(run_gencase femus_gencase)

configure_file(run_gencase.sh.in ${CMAKE_BINARY_DIR}/run_gencase.sh @ONLY)
