// ************************************************************************************************
//                 MGSolverR_0.C
// ************************************************************************************************
//     Constructor -> MGSolR::MGSolR
//     Destructor  -> MGSolR::~MGSolR
//     Init
// ************************************************************************************************
//    TIME STEP->  T1) MGSolR::MGTimeStep
//                 T2) MGSolR::MGTimeStep_no_up
//                 T3) MGSolR::MGTimeStep_solve

// ===============================================================
#include "MGSolver_L3_R.h"       // Navier-Stokes class header file

// // #include "MGFE_conf.h"        // FEM approximation
// #include "MGGeomEl.h"        // FEM approximation
// #include "Printinfo_conf.h"  // Print options
// 
// #include "MeshExtended.h"
// // #include "MGFE.h"          // Mesh class
// #include "EquationSystemsExtendedM.h"  // Equation map class
// 
// // local alg lib -----------------------------------------------
// #include "sparse_matrixM.h"  // algebra sparse matrices
// #include "numeric_vectorM.h" // algebra numerical vectors
// #include "linear_solverM.h"  // algebra solvers




// ************************************************************************************************
// ************************************************************************************************
//                Constructor/Destructor/Init
//                c1) MGSolR::MGSolR 
//                c2) MGSolR::~MGSolR
// ************************************************************************************************
// ==================================================================
/// This routine constructs the FSI class:
MGSolR::MGSolR (
 MGEquationsSystem & mg_equations_map_in,
 int             nvars_in[],
 std::string     eqname_in,
 std::string     varname_in
) :  MGSolDA ( mg_equations_map_in, nvars_in, eqname_in, varname_in )
{// ==============================================================================================
   return;
} 


// ================================================================================================
MGSolR::~MGSolR (
){// ==============================================================================================
     delete[] _FeMBC_R;     ///< local rhs          
    delete[] _KeMBC_R;     ///< local matrix
  return;
}

// ************************************************************************************************
//               TIME STEP
// ************************************************************************************************
//    TIME STEP->  T1) MGSolR::MGTimeStep
//                 T2) MGSolR::MGTimeStep_no_up
//                 T3) MGSolR::MGTimeStep_solve
// ************************************************************************************************

// ================================================================================================
int MGSolR::MGTimeStep(
 double const time,     ///< time (double)
 int const    i_step,   ///< time (int)
 int const    max_iter, ///< max num iterations
 int const    mode      ///< mode (0=only solve;1=only up; 2 solve+up
) { // ============================================================================================
  int err=0;                           
//   if(mode%2==0)  err=  MGTimeStep_solve(time,max_iter);
if(mode%2==0)  err=  MGTimeStep_no_up(time,i_step,max_iter);// mode=0 or 2
 if(mode>0) MGUpdateStep();                                // mode=1 or 2)
   return err;
 }

// ================================================================================================
int   MGSolR::MGTimeStep_no_up(// mode=0 or 2
 double const time,    ///< time (double)
 int const  i_step,    ///< time (int)
 int const  max_iter   ///< max num iterations
) {// ============================================================================================
 int err= MGTimeStep_solve(time,max_iter);
  return err;
}


// ================================================================================================
/// 3) This function controls the assembly and the solution of the NS_equation system:
int MGSolR::MGTimeStep_solve(
 const double time,  ///< time
 const int max_iter  // Number of max inter
) { // ============================================================================================

// ========================================================================================= //
//              A) Set up the time step                                                      //
// ========================================================================================= //
 std::cout  << std::endl << "\033[038;5;" << NS_F + 50 << ";1m "
            << "--------------------------------------------------- \n\t"
            <<  _eqname_D.c_str()
            << " solution of problem " << _mgutils_D.get_name() << " with dir " << _dir_D
            << "\n ---------------------------------------------------\n\033[0m";
//                 x[_NoLevels - 1]      -> localize(*x_old[_NoLevels - 1]);       // time step
//      x_old[_NoLevels - 1]  -> localize(*x_oold[_NoLevels - 1]);     // time step -1
 // ========================================================================================= //
 //              B) Assemblying of the Matrix-Rhs                                             //
 // ========================================================================================= //
TIME(std::clock_t start_time = std::clock();)
 GenMatRhs(time, _NoLevels - 1, 1);                                          // matrix and rhs
//     printf("disp norm genmat is %20f %d \n",disp[0]->l2_norm(),_iproc);
 for(int Level = 0 ; Level < _NoLevels - 1; Level++) {
  GenMatRhs(time, Level, 0);      // matrix
 }
TIME( std::clock_t end_time = std::clock();std::cout << "  Assembly time -----> =" 
                    << double(end_time - start_time) / CLOCKS_PER_SEC << " s \n";)
 // ========================================================================================= //
 //              C) Solution of the linear MGsystem (MGSolFSI::MGSolve)                       //
 // ========================================================================================= //
int err= MGSolve_B(1.e-6, 15);
TIME(end_time = std::clock(); std::cout << " Assembly+solution time -----> =" 
                        << double(end_time - start_time)/CLOCKS_PER_SEC<< "s \n";)

 return err;
}


