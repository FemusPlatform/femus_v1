// #ifndef __solverlibconf_h__
// #define __solverlibconf_h__

#pragma  once
// ************************************************************************************************
//  SOLVER LIBRARY 
// ************************************************************************************************
// #define NUM_MESH (1)
// #define NDOF_K (1)
// #define MATBC_INTERFACE
// #define XFEM

// PETSCM libs ------------------------------------------------------------------------------------
#define HAVE_PETSCM 
#define HAVE_MPI   //what if I want to use Petsc without MPI?

// LASPACKM libs ----------------------------------------------------------------------------------
//   #define HAVE_LASPACKM
#define LSOLVER LASPACK_SOLVERSM  
#ifdef HAVE_PETSCM
#undef LSOLVER
 #define LSOLVER  PETSC_SOLVERSM
#endif

// Med library ------------------------------------------------------------------------------------
#define    HAVE_MED

// Vanka solver -----------------------------------------------------------------------------------
//  #define VANKA (0)



// ************************************************************************************************
//  OPTIONS (do not change below)
// ************************************************************************************************

//******PETSC VERSION ***88
#ifdef HAVE_PETSCM


/* PETSc's major version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MAJOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MAJOR  3 
#endif

/* PETSc's minor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MINOR  10
#endif

/* PETSc's subminor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR  2 
#endif

#endif


// HDF5 -------------------------------------------------------------------------------------------
// #define HDF5_VERSIONM 1812
// #define HDF5_VERSIONM 1810

// #define HDF5_VERSIONM 188



// #endif
