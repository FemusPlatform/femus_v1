cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(femus_lib
  VERSION 0.0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG main
)
FetchContent_MakeAvailable(fmt)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(FeatureSummary)

find_package(Eigen3 REQUIRED NO_MODULE)
# femus libraries

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

if (EXISTS $ENV{PLAT_CODES_DIR})
  message(STATUS "setting FEMUS_DIR from PLAT_CODES_DIR")
  set(FEMUS_DIR
    $ENV{PLAT_CODES_DIR}/femus
    CACHE PATH "The femus library location."
  )
elseif(EXISTS $ENV{FEMUS_DIR})
message(STATUS "setting FEMUS_DIR from ENV FEMUS_DIR")
  set(FEMUS_DIR
    $ENV{FEMUS_DIR}
    CACHE PATH "The femus library location."
  )
endif()

if (NOT EXISTS ${FEMUS_DIR})
  message(FATAL_ERROR "could not find femus directory at ${FEMUS_DIR}, please set FEMUS_DIR.")
else()
  message(STATUS "femus location: ${FEMUS_DIR}")
endif()

string(REPLACE femus_app_ "" APP_NAME ${PROJECT_NAME})
set(APP_DIR ${CMAKE_SOURCE_DIR}/..)

add_subdirectory(${FEMUS_DIR} femus)


target_link_libraries(femus PUBLIC Eigen3::Eigen)

# add_executable(${PROJECT_NAME})

# # ------------------- gprof -------------------
# option(ENABLE_GPROF "Enable gprof profiling (-pg)" ON)

# if (ENABLE_GPROF)
#   message(STATUS "Building with gprof (-pg)")
#   target_compile_options(${PROJECT_NAME} PRIVATE -pg)
#   target_link_options(${PROJECT_NAME} PRIVATE -pg)

#   # profila anche la libreria femus (importante se il tempo sta lì)
#   target_compile_options(femus PRIVATE -pg)
#   target_link_options(femus PRIVATE -pg)
# endif()
# # ---------------------------------------------


# target_include_directories(${PROJECT_NAME}
#   PUBLIC
#     ${FEMUS_DIR}/include
#     ${FEMUS_DIR}/include/MGfemus  
#     ${FEMUS_DIR}/include/MGmesh    
#     ${FEMUS_DIR}/include/MGsystem  
#     ${FEMUS_DIR}/include/MGutility
#     ${FEMUS_DIR}/include/MGfem         
#     ${FEMUS_DIR}/include/MGgeom   
#     ${FEMUS_DIR}/include/MGmed         
#     ${FEMUS_DIR}/include/Mgsolver  
#     ${FEMUS_DIR}/include/MGturb
#     ${FEMUS_DIR}/contrib/matrix
#     ${FEMUS_DIR}/contrib/parallel
#     ${MEDCOUPLING_INCLUDE_DIRS}
#     $ENV{PLAT_THIRD_PARTY_DIR}/medcoupling/include
#     $ENV{PLAT_THIRD_PARTY_DIR}/openmpi/include
# )


# target_link_libraries(${PROJECT_NAME}
#   PUBLIC
#     femus
#     fmt::fmt
# )

# set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# configure_file(cmake/app_setup.h.in ${CMAKE_BINARY_DIR}/app_setup.h @ONLY)
# configure_file(cmake/run_app.sh.in ${CMAKE_BINARY_DIR}/run_app.sh @ONLY)
# configure_file(cmake/femus_env.sh.in ${CMAKE_BINARY_DIR}/femus_env.sh @ONLY)

# if (EXISTS $ENV{PLAT_CODES_DIR})
#   set(FEMUS_GENCASE_DIR
#     $ENV{PLAT_CODES_DIR}/femus/applications/gencase
#     CACHE PATH "The location of the gencase application."
#   )
# elseif(EXISTS $ENV{FEMUS_GENCASE_DIR})
#   set(FEMUS_GENCASE_DIR
#     $ENV{FEMUS_GENCASE_DIR}
#     CACHE PATH "The location of the gencase application."
#   )
# endif()

# if (NOT EXISTS ${FEMUS_GENCASE_DIR})
#   message(FATAL_ERROR "could not find gencase directory at ${FEMUS_GENCASE_DIR}, please set FEMUS_GENCASE_DIR.")
# else()
#   message(STATUS "gencase location: ${FEMUS_GENCASE_DIR}")
# endif()
# 
# add_subdirectory(${FEMUS_GENCASE_DIR} gencase)

feature_summary(WHAT ALL)
