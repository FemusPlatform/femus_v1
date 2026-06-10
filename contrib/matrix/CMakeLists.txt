add_library(femus_matrix SHARED)
set_target_properties(femus_matrix PROPERTIES LINKER_LANGUAGE CXX)
# target_compile_definitions(femus_matrix
#   PUBLIC
#     DIMENSION=2
#     HAVE_PETSCM=1
#     LSOLVER=PETSC_SOLVERSM
#     LIBMESH_DETECTED_PETSC_VERSION_MAJOR=3
#     LIBMESH_DETECTED_PETSC_VERSION_MINOR=13
#     LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR=6
# )
target_include_directories(femus_matrix
  PUBLIC
    ${FEMUS_DIR}/contrib/parallel
    ${FEMUS_DIR}/include
    ${APP_DIR}/DATA
)
target_link_libraries(femus_matrix
  PUBLIC
    femus_parallel
    hdf5::hdf5-shared
    PkgConfig::petsc
    ${PETSC_EXTRA_LIBS}
)

set(FEMUS_MATRIX_HEADERS
  ${FEMUS_DIR}/contrib/matrix/utilityM.h
  # ${FEMUS_DIR}/contrib/matrix/type_vector.h
  ${FEMUS_DIR}/contrib/matrix/sparse_MmatrixM.h
  ${FEMUS_DIR}/contrib/matrix/sparse_matrixM.h
  ${FEMUS_DIR}/contrib/matrix/Solvertype_enum.h
  ${FEMUS_DIR}/contrib/matrix/SolverPackage_enum.h
  ${FEMUS_DIR}/contrib/matrix/Precondtype_enum.h
  ${FEMUS_DIR}/contrib/matrix/preconditionerM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_vectorM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_preconditionerM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_MmatrixM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_matrixM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_macroM.h
  ${FEMUS_DIR}/contrib/matrix/petsc_linear_solverM.h
  ${FEMUS_DIR}/contrib/matrix/Paralleltype_enum.h
  ${FEMUS_DIR}/contrib/matrix/numeric_vectorM.h
  ${FEMUS_DIR}/contrib/matrix/linear_solverM.h
  ${FEMUS_DIR}/contrib/matrix/laspack_vectorM.h
  ${FEMUS_DIR}/contrib/matrix/laspack_MmatrixM.h
  ${FEMUS_DIR}/contrib/matrix/laspack_matrixM.h
  ${FEMUS_DIR}/contrib/matrix/laspack_linear_solverM.h
  ${FEMUS_DIR}/contrib/matrix/enum_subset_solve_modeM.h
  ${FEMUS_DIR}/contrib/matrix/distributed_vectorM.h
  ${FEMUS_DIR}/contrib/matrix/dense_vectorM.h
  ${FEMUS_DIR}/contrib/matrix/dense_vector_baseM.h
  ${FEMUS_DIR}/contrib/matrix/dense_subvectorM.h
  ${FEMUS_DIR}/contrib/matrix/dense_submatrixM.h
  ${FEMUS_DIR}/contrib/matrix/dense_set.h
  ${FEMUS_DIR}/contrib/matrix/dense_matrixM.h
  ${FEMUS_DIR}/contrib/matrix/dense_matrix_baseM.h
)

set(FEMUS_MATRIX_SOURCES
  # ${FEMUS_DIR}/contrib/matrix/type_vector.C
  ${FEMUS_DIR}/contrib/matrix/sparse_MmatrixM.C
  ${FEMUS_DIR}/contrib/matrix/sparse_matrixM.C
  ${FEMUS_DIR}/contrib/matrix/preconditionerM.C
  ${FEMUS_DIR}/contrib/matrix/petsc_vectorM.C
  ${FEMUS_DIR}/contrib/matrix/petsc_preconditionerM.C
  ${FEMUS_DIR}/contrib/matrix/petsc_MmatrixM.C
  ${FEMUS_DIR}/contrib/matrix/petsc_matrixM.C
  ${FEMUS_DIR}/contrib/matrix/petsc_linear_solverM.C
  ${FEMUS_DIR}/contrib/matrix/numeric_vectorM.C
  ${FEMUS_DIR}/contrib/matrix/linear_solverM.C
  ${FEMUS_DIR}/contrib/matrix/laspack_vectorM.C
  ${FEMUS_DIR}/contrib/matrix/laspack_MmatrixM.C
  ${FEMUS_DIR}/contrib/matrix/laspack_matrixM.C
  ${FEMUS_DIR}/contrib/matrix/laspack_linear_solverM.C
  ${FEMUS_DIR}/contrib/matrix/distributed_vectorM.C
  ${FEMUS_DIR}/contrib/matrix/dense_vectorM.C
  ${FEMUS_DIR}/contrib/matrix/dense_vector_baseM.C
  ${FEMUS_DIR}/contrib/matrix/dense_subvectorM.C
  ${FEMUS_DIR}/contrib/matrix/dense_submatrixM.C
  ${FEMUS_DIR}/contrib/matrix/dense_set.C
  ${FEMUS_DIR}/contrib/matrix/dense_matrixM.C
  ${FEMUS_DIR}/contrib/matrix/dense_matrix_baseM.C
)

target_sources(femus_matrix
  PUBLIC ${FEMUS_MATRIX_HEADERS}
  PRIVATE ${FEMUS_MATRIX_SOURCES}
)

set_source_files_properties(${FEMUS_MATRIX_HEADERS} ${FEMUS_MATRIX_SOURCES} PROPERTIES
  COMPILE_DEFINITIONS
    DIMENSION=2
  COMPILE_DEFINITIONS
    HAVE_PETSCM=1
  COMPILE_DEFINITIONS
    LSOLVER=PETSC_SOLVERSM
  COMPILE_DEFINITIONS
    LIBMESH_DETECTED_PETSC_VERSION_MAJOR=3
  COMPILE_DEFINITIONS
    LIBMESH_DETECTED_PETSC_VERSION_MINOR=13
  COMPILE_DEFINITIONS
    LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR=6
)

