add_library(femus_parallel SHARED)
set_target_properties(femus_parallel PROPERTIES LINKER_LANGUAGE CXX)
# target_compile_definitions(femus_parallel PUBLIC DIMENSION=2 HAVE_MPI=1)
target_include_directories(femus_parallel
  PUBLIC
    ${APP_DIR}/DATA
)

set(FEMUS_PARALLEL_HEADERS
  ${FEMUS_DIR}/contrib/parallel/parallel_algebraM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_bin_sorterM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_conversion_utilsM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_elemM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_ghost_syncM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_hilbertM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_histogramM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_implementationM.h
  ${FEMUS_DIR}/contrib/parallel/parallelM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_nodeM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_objectM.h
  ${FEMUS_DIR}/contrib/parallel/parallel_sortM.h
  ${FEMUS_DIR}/contrib/parallel/threads_allocatorsM.h
  ${FEMUS_DIR}/contrib/parallel/threadsM.h
)

set(FEMUS_PARALLEL_SOURCES "")

target_sources(femus_parallel
  PUBLIC ${FEMUS_PARALLEL_HEADERS}
  PRIVATE ${FEMUS_PARALLEL_SOURCES}
)

target_link_libraries(femus_parallel
  PUBLIC
    MPI::MPI_CXX
)

