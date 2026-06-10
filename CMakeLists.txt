if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
   message(FATAL_ERROR "femus cannot (yet) be compiled standalone, compile it from the app directory.")
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

find_package(MEDCoupling
  REQUIRED
  PATHS
    $ENV{SPACK_ENV}/.spack-env/view/cmake_files
    $ENV{PLAT_THIRD_PARTY_DIR}/medcoupling/cmake_files
)
add_library(mymedcoupling INTERFACE)
# fix missing target in medcoupling cmake
target_link_libraries(mymedcoupling
  INTERFACE
    interpkernel
    medcouplingcpp
    medcouplingremapper
    # medicoco
    medloader
)
target_include_directories(mymedcoupling
  INTERFACE
    ${MEDCOUPLING_INCLUDE_DIRS}
)
add_library(medcoupling::medcoupling ALIAS mymedcoupling)

find_package(PkgConfig REQUIRED)
# set(CMAKE_PREFIX_PATH $ENV{PLAT_THIRD_PARTY_DIR}/petsc/lib/pkgconfig)
set(ENV{PKG_CONFIG_PATH} $ENV{PLAT_THIRD_PARTY_DIR}/petsc/lib/pkgconfig)
pkg_check_modules(petsc REQUIRED IMPORTED_TARGET PETSc)

# petsc pkgconf does not include static dependencies
set(PETSC_EXTRA_LIBS)
if(DEFINED ENV{PLAT_DIR})
  list(APPEND PETSC_EXTRA_LIBS
    ${petsc_LIBRARY_DIRS}/libfblas.a
    ${petsc_LIBRARY_DIRS}/libflapack.a
    ${petsc_LIBRARY_DIRS}/libHYPRE.so
  )
elseif(DEFINED ENV{SPACK_ENV})
  list(APPEND PETSC_EXTRA_LIBS
  $ENV{SPACK_ENV}/.spack-env/view/lib/libopenblas.so
)
endif()

# contrib

add_subdirectory(contrib/parallel)
add_subdirectory(contrib/matrix)

add_library(femus)

target_include_directories(femus
  PUBLIC
    ${FEMUS_DIR}/include
    ${FEMUS_DIR}/include/MGfemus  
    ${FEMUS_DIR}/include/MGmesh    
    ${FEMUS_DIR}/include/MGsystem  
    ${FEMUS_DIR}/include/MGutility
    ${FEMUS_DIR}/include/MGfem         
    ${FEMUS_DIR}/include/MGgeom   
    ${FEMUS_DIR}/include/MGmed         
    ${FEMUS_DIR}/include/Mgsolver  
    ${FEMUS_DIR}/include/MGturb
    ${FEMUS_DIR}/include/MGquadrature
    ${FEMUS_DIR}/contrib/matrix
    ${FEMUS_DIR}/contrib/parallel
    
    
    ${APP_DIR}/DATA
    ${APP_DIR}/SRC
    ${APP_DIR}/SRC/NS
    #${APP_DIR}/SRC/NSA
    ${APP_DIR}/SRC/T
    ${APP_DIR}/SRC/TA
    ${APP_DIR}/SRC/CTRL
    ${APP_DIR}/SRC/COL
    ${APP_DIR}/SRC/KK
    ${APP_DIR}/SRC/TFM
    ${APP_DIR}/SRC/CUTFEM
)

target_link_libraries(femus
  PUBLIC
    femus_parallel
    femus_matrix
    medcoupling::medcoupling
)

set_target_properties(femus PROPERTIES LINKER_LANGUAGE CXX)

set(FEMUS_HEADERS
       ${FEMUS_DIR}/include/MGCasts.h  
       ${FEMUS_DIR}/include/MGGraph.h
    # MGgeom   
       ${FEMUS_DIR}/include/MGgeom/MGGeom_L0_El.h  
       ${FEMUS_DIR}/include/MGgeom/MGGeom_L1_El.h
       ${FEMUS_DIR}/include/MGgeom/MGGeom_L1_Map.h
       # MGfem
      ${FEMUS_DIR}/include/MGfem/MGFE_L0_II_enum_order.h  
       ${FEMUS_DIR}/include/MGfem/MGFE_L0_abstract.h  
       ${FEMUS_DIR}/include/MGfem/MGFE_L0_Map.h
       ${FEMUS_DIR}/include/MGfem/MGFE_L0_C.h  
       ${FEMUS_DIR}/include/MGfem/MGFE_L1_base.h           
       ${FEMUS_DIR}/include/MGfem/MGFE_L2_FE_template.h  
       ${FEMUS_DIR}/include/MGfem/MGFE_L2_lagrange_shape_1D.h 
     # MGfemus  
 ${FEMUS_DIR}/include/MGfemus/MGFemus_L0.h  
 ${FEMUS_DIR}/include/MGfemus/MGFemus_L0_Init.h  
 ${FEMUS_DIR}/include/MGfemus/MGFemus_L1.h
#  MGmesh
 ${FEMUS_DIR}/include/MGmesh/MGMesh_L0.h  
 ${FEMUS_DIR}/include/MGmesh/MGMesh_L1_Extended.h
 #         MGquadrature
 ${FEMUS_DIR}/include/MGquadrature/MGenum_quadrature_type.h
 ${FEMUS_DIR}/include/MGquadrature/MGquadrature_L1_conical.h 
 ${FEMUS_DIR}/include/MGquadrature/MGquadrature_L1_jacobi.h
 ${FEMUS_DIR}/include/MGquadrature/MGquadrature_L1_gauss.h  
 ${FEMUS_DIR}/include/MGquadrature/MGquadrature_L0.h        
 ${FEMUS_DIR}/include/MGquadrature/MGquadrature_L1_trap.h 
#         MGsystem
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L0.h  
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L1_II.h  
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L2_EquationsMap.h
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L1.h  
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L2.h     
 ${FEMUS_DIR}/include/MGsystem/MGSystem_L3_TimeLoop.h
 # MGutility
       ${FEMUS_DIR}/include/MGutility/MGUtils_L0.h  
       ${FEMUS_DIR}/include/MGutility/MGUtils_L1_InterfaceFunctionM.h  
       ${FEMUS_DIR}/include/MGutility/MGUtils_L1_InterfaceProjection.h
       # MGmed
       ${FEMUS_DIR}/include/MGmed/MGMed_L0.h
       # Mgsolver
       ${FEMUS_DIR}/include/Mgsolver/MGSolver_L0_Base.h  
       ${FEMUS_DIR}/include/Mgsolver/MGSolver_L1_Base.h  
       ${FEMUS_DIR}/include/Mgsolver/MGSolver_L3_R.h
       ${FEMUS_DIR}/include/MGturb/MGTurbModels_L0.h  
       ${FEMUS_DIR}/include/MGturb/MGTurbModels_L0_II.h  
       ${FEMUS_DIR}/include/MGturb/MGTurbModels_L0_III.h
       
)

set(FEMUS_SOURCES

# MGfem
 ${FEMUS_DIR}/src/MGfem/MGFE_L0_abstract.C       
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange.C           
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_2D.C
 ${FEMUS_DIR}/src/MGfem/MGFE_L0_C.C 
 ${FEMUS_DIR}/src/MGfem/MGFE_L0_Map.C
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_0D.C  
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_3D.C
 ${FEMUS_DIR}/src/MGfem/MGFE_L1_base.C  
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_1D.C  
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_FE_template.C
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_xfem.C 
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_xfem_shape_2D.C
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_xfem_shape_0D.C  
 ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_xfem_shape_3D.C
  ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_xfem_shape_1D.C  
# MGfemus
 ${FEMUS_DIR}/src//MGfemus/MGFemus_L0_Init.C
 ${FEMUS_DIR}/src//MGfemus/MGFemus_L1.C  
 ${FEMUS_DIR}/src//MGfemus/MGFemus_L1_med.C
# MGgeom
 ${FEMUS_DIR}/src/MGgeom/MGGeom_L0_El.C 
 ${FEMUS_DIR}/src/MGgeom/MGGeom_L1_ElLag.C
 ${FEMUS_DIR}/src/MGgeom/MGGeom_L1_ElLag10.C  
 ${FEMUS_DIR}/src/MGgeom/MGGeom_L1_ElLag27.C
# MGmesh
 ${FEMUS_DIR}/src/MGmesh/MGMesh_L0.C  
 ${FEMUS_DIR}/src/MGmesh/MGMesh_L0_PR.C  
 ${FEMUS_DIR}/src/MGmesh/MGMesh_L1_Extended.C  
 # ${FEMUS_DIR}/src/MGmesh/MGMesh_L1_Extended_II.C
 
 #         MGquadrature
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L0.C  
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L0_build.C 
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_jacobi.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_jacobi_1D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_conical.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_conical_2D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_conical_3D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_gauss.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_gauss_1D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_gauss_2D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_gauss_3D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_trap.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_trap_1D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_trap_2D.C
 ${FEMUS_DIR}/src/MGquadrature/MGquadrature_L1_trap_3D.C
# MGsolver
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L0_Base.C         
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L0_Base_VANKA.C   
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L1_Base.C   
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA.C        
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_SOL.C 
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_check.C  
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_ASSEMBLY.C  
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_init.C   
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_BCIC.C   
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L2_DA_PR.C  
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L3_R_0.C
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L3_R_BIC.C
 ${FEMUS_DIR}/src//MGsolver/MGSolver_L3_R_CHECK.C
#  MGsystem
${FEMUS_DIR}/solvers/EqMapFiller/MGSystem_L2_EquationsMap.C
 ${FEMUS_DIR}/src/MGsystem/MGSystem_L0.C  
 ${FEMUS_DIR}/src/MGsystem/MGSystem_L1_eq.C  
 ${FEMUS_DIR}/src/MGsystem/MGSystem_L2_ExtendedM.C  
 ${FEMUS_DIR}/src/MGsystem/MGSystem_L3_TimeLoop.C
# MGutility
 ${FEMUS_DIR}/src/MGutility/MGUtils_L0.C      
 ${FEMUS_DIR}/src/MGutility/MGUtils_L1_Femus_InterfaceFunctionM.C
 ${FEMUS_DIR}/src/MGutility/MGUtils_L0_Med.C  
 ${FEMUS_DIR}/src/MGutility/MGUtils_L1_Femus_InterfaceProjection.C



#   # src/FEMuS
#   ${FEMUS_DIR}/src/MGfemus/MGFemus.C
#   ${FEMUS_DIR}/src/MGfemus/MGFemus_med.C
#    ${FEMUS_DIR}/src/MGfemus/MGFemusInit.C
#   
#   ${FEMUS_DIR}/src/MGgeom/MGGeom_L0_El.C
#   ${FEMUS_DIR}/src/MGGeom/MGGeom_L1_ElLag27.C
#   ${FEMUS_DIR}/src/MGfem/MGGeom_L1_ElLag10.C
#    
#   ${FEMUS_DIR}/src/MGfem/MGFE_L0.C    
#   ${FEMUS_DIR}/src/MGfem/MGFE_L0_C.C  
#   ${FEMUS_DIR}/src/MGfem/MGFE_L1_base.C 
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange.C           
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_2D.C  
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_0D.C  
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_3D.C  
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_lagrange_shape_1D.C  
#   ${FEMUS_DIR}/src/MGfem/MGFE_L2_template.C           
#   
#  ${FEMUS_DIR}/src/MGmesh/MGMesh_L0.C  
#  ${FEMUS_DIR}/src/MGmesh/MGMesh_L0_PR.C  
#  ${FEMUS_DIR}/src/MGmesh/MGMesh_L1_Extended.C
#  ${FEMUS_DIR}/src/MGmesh/MGMesh_L1_Extended_II.C  
 
 
 
 
 # ${FEMUS_DIR}/src/MGfemus/MGFemus_InterfaceFunctionM.C
 #  ${FEMUS_DIR}/src/MGfemus/MGFemus_InterfaceProjection.C
  
 
  # ${FEMUS_DIR}/src/MGfem/MGFE1.C
  # ${FEMUS_DIR}/src/MGfem/MGFE2.C
  # ${FEMUS_DIR}/src/MGfem/MGFE1_GET.C
  #  ${FEMUS_DIR}/src/MGfem/MGFE1_IO.C
  # ${FEMUS_DIR}/src/MGfem/MGFE1_JAC.C
#   ${FEMUS_DIR}/src/MGfem/MGFE0.C   
#   # ${FEMUS_DIR}/src/MGfem/MGFE1_JAC.C  
#   ${FEMUS_DIR}/src/MGfem/MGfe_lagrange.C           
#   ${FEMUS_DIR}/src/MGfem/MGfe_lagrange_shape_3D.C  
#   ${FEMUS_DIR}/src/MGfem/MGGeomElLag27.C
# # ${FEMUS_DIR}/src/MGfem/MGFE1.C      
# ${FEMUS_DIR}/src/MGfem/MGFEbase.C   
# ${FEMUS_DIR}/src/MGfem/MGfe_lagrange_shape_0D.C  
# ${FEMUS_DIR}/src/MGfem/MGGeomEl0.C
# # ${FEMUS_DIR}/src/MGfem/MGFE1_GET.C  
# ${FEMUS_DIR}/src/MGfem/MGfe.C       
# ${FEMUS_DIR}/src/MGfem/MGfe_lagrange_shape_1D.C  
# # ${FEMUS_DIR}/src/MGfem/MGFE1_IO.C   
# ${FEMUS_DIR}/src/MGfem/MGFE_C.C     
# ${FEMUS_DIR}/src/MGfem/MGfe_lagrange_shape_2D.C  
# ${FEMUS_DIR}/src/MGfem/MGGeomElLag10.C
#   # ${FEMUS_DIR}/src/MGfem/MGGeomEl.C
  # ${FEMUS_DIR}/src/MGMesh/MeshExtended.C
  # ${FEMUS_DIR}/src/MGMesh/MGMesh.C
  # ${FEMUS_DIR}/src/MGMesh/MGMesh_PR.C
  # ${FEMUS_DIR}/src/MGSolverBase/MGSolverBase0.C
  # ${FEMUS_DIR}/src/MGSolverBase/MGSolverBase0_VANKA.C
  # ${FEMUS_DIR}/src/MGSolverBase/MGSolverBase.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_ASSEMBLY.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_BCIC.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_check.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_init.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_PR.C
  # ${FEMUS_DIR}/src/MGSolverDA/MGSolverDA_SOL.C
  # ${FEMUS_DIR}/src/MGSolverR/MGSolverR_0.C
  # ${FEMUS_DIR}/src/MGSolverR/MGSolverR_BIC.C
  # ${FEMUS_DIR}/src/MGSolverR/MGSolverR_CHECK.C
  # ${FEMUS_DIR}/src/MGSystem/EquationSystemsExtendedM.C
  # ${FEMUS_DIR}/src/MGSystem/MGEquationsSystem.C
  # ${FEMUS_DIR}/src/MGSystem/MGSystem.C
  # ${FEMUS_DIR}/src/MGSystem/MGTimeLoop.C
  # ${FEMUS_DIR}/solvers/EqMapFiller/EquationsMap.C
  # ${FEMUS_DIR}/src/Utility_lib/MedUtils.C
  # ${FEMUS_DIR}/src/Utility_lib/MGUtils.C
)

target_sources(femus
  PUBLIC
    ${FEMUS_HEADERS}
  PRIVATE
    ${FEMUS_SOURCES}
)

# file(COPY ${CMAKE_CURRENT_BINARY_DIR}/libfemus.so DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib/libfemus.so)
add_custom_command(
 TARGET femus POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_BINARY_DIR}/libfemus.so
   ${FEMUS_DIR}/lib/libfemus.so
)

add_custom_command(
 TARGET femus POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_BINARY_DIR}/contrib/matrix/libfemus_matrix.so
   ${FEMUS_DIR}/lib/libfemus_matrix.so
)

add_custom_command(
 TARGET femus POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_BINARY_DIR}/contrib/parallel/libfemus_parallel.so
   ${FEMUS_DIR}/lib/libfemus_parallel.so
)
