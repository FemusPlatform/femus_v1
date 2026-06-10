
// c++ class include files ----------------------------------------------------
#include <string>

// alg include-------------------
#include "linear_solverM.h"
// #include "numeric_vectorM.h"
// #include "sparse_MmatrixM.h"
// #include "sparse_matrixM.h"

// config files
#include "Printinfo_conf.h" // for TIME() and TRACKING_FUN

// local inlude -----------------
// #include "MGMesh_L0.h"
#include "MGSolver_L0_Base.h"
#include "MGMesh_L0.h"


// ===============================================================================================
//                                     MULTILEVEL OPERATOR
//================================================================================================
// ============================================================================
/// This function  is the MGSolBase0 constructor :
MGSolBase0::MGSolBase0(
MGMesh&  mgmesh_in  ///<  MG equation map
//   std::string eq_name_in      ///< equation name (default Base0)
) :
  _mgmesh(mgmesh_in),    // mgmesh pointer from equation map pointer
  _NoLevels((int)(_mgmesh._NoLevels)) {

  // allocation of dynamic system ---------------------------------------------
  _Nnodes=_mgmesh._NoNodes[_NoLevels-1];
  _Ndof_lev = new int[_NoLevels];   // matrix and vect  dim
  _node_dof = new int*[_NoLevels];  // dof (+1)
  A.resize(_NoLevels);         // matrix
  x.resize(_NoLevels);         // matrix vect sol
  b.resize(_NoLevels);         // rhs
  res.resize(_NoLevels);       // residual
  // restr and prol operators -------------------------------------------------
  Rst.resize(_NoLevels);  Prl.resize(_NoLevels);  // Projector and restrictor
  // solver -------------------------------------------------------------------
  _iproc = _mgmesh._iproc;  //  processor number (parallel >1CPU)
  _solver = new LinearSolverM*[_NoLevels];  // one solver for each level
  for(int l = 0; l < _NoLevels; l++) _solver[l] = LinearSolverM::build(_mgmesh._comm.comm(), LSOLVER).release();
  _Solve_B=true;
  return;
}


// ===================================================
/// This function  is the MGSolBase0 destructor.
//-------------------------------------------------------------------------
MGSolBase0::~MGSolBase0() {
  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::MGSolBase (Destructor ******) \n  ");)
// ======================================================================================
  // clear substructrures
  clear();
  A.clear();  x.clear();  b.clear();    //  A and x and b
  res.clear();  //  rhs and residual vector
  Rst.clear();  Prl.clear();  // Restrictor and projector
  delete[] _Ndof_lev;  // dimension system Ax=b
  delete[] _solver;     // delete solver;
  delete[] _node_dof;
// ======================================================================================
  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::MGSolBase (Destructor ******) \n  ");)
}  ///< Destructor (level structure)




// =================================================================
/// This function  is the substructure MGSolBase0 destructor (clear).
void MGSolBase0::clear() {
  TRACKING_FUN(printf(" --> MGSolverBase0.C: MGSolBase0::clear \n  ");)
// ======================================================================================
  for(int Level = 0; Level < _NoLevels; Level++) {
    delete A[Level];  delete x[Level]; delete b[Level]; //  A and x b at Level
    delete res[Level];  //  residual  at Level
    if(Level < _NoLevels - 1)  delete Rst[Level];    // Restrictor
      if(Level > 0)       delete Prl[Level];    // projector
      delete _solver[Level];                         // delete solver  at Level
      delete[] _node_dof[Level]; // dof distribution at Level
    }
// ======================================================================================
    TRACKING_FUN(printf(" <-- MGSolverBase0.C: MGSolBase0::clear \n  ");)
  }  ///< Substructure destructor






// ===============================================================================================
//                                     MULTILEVEL OPERATOR
//================================================================================================


// ====================================================================
/// This function solves the discrete problem with multigrid solver
  int MGSolBase0::MGSolve_B(
  double Eps1,  // tolerance
  int MaxIter,  // n iterations
// -----------------------------------------------------------
  const int clearing,   ///< 1= clean the init matrix and precond
  const int Gamma,      ///< Control V W cycle
  const int Nc_pre,     ///< n pre-smoothing cycles
  const int Nc_coarse,  ///< n coarse cycles
  const int Nc_post     ///< n post-smoothing cycles
  ) {
    TRACKING_FUN(printf(" --> MGSolverBase0.C: MGSolBase0::MGSolve \n  ");)
    // ===================================================================
    double rest = 0.;
    b[_NoLevels - 1]->close();
    double bNorm = 0.;
    bNorm = b[_NoLevels - 1]->l2_norm();

//    P_CONV(std::cout <<" bNorm2 "<<bNorm<</*" bNorm1 "<<bNorm1<<" bNorm inf "<<bNorm0<<*/"\n";)
    if(bNorm != bNorm) return -2; // check if bNorm is NaN
    if(bNorm>1.e+15) return -1;

    x[_NoLevels - 1]->close();
    x[_NoLevels - 1]->zero();

    // FAS Multigrid (Nested)
    int NestedMG = 1;
    if(NestedMG == 0) { // NestedMG = 1 -> no FAS Multigrid
    x[_NoLevels - 1]->close();
    x[_NoLevels - 1]->zero();
    rest = MGStep_B(0, 1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);
    for(int Level = 1; Level < _NoLevels; Level++) {
       x[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);                                 // projection
       rest = MGStep_B(Level,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);  // MGStepsolution
     }
   }

    // V or W cycle
    int cycle = 0;  int exit_mg = 0;
    while(exit_mg == 0 && cycle < MaxIter) {   // multigrid step start
#ifndef VANKA
      rest = MGStep_B(_NoLevels - 1, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post, clearing);  // MGStep
      DEBUG_MGSOLB2(rest=0.; MGCheck(_NoLevels-1);)
#else
      rest = MGStep_Vanka_B(_NoLevels - 1, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse,Nc_post);  // MGStep_Vanka
#endif
      if(rest>1.e+20) {std::cout << "  residual is too high " << std::endl; return-1;}
      if(rest < Eps1*(1. + bNorm)) exit_mg = 1;   // exit test
      cycle++;
      std::cout << " cycle= " << cycle << " residual= " << rest << " \n";
    }  // --------------------- end multigrid step
    TRACKING_FUN(printf(" <-- MGSolverBase0.C: MGSolBase0::MGSolve \n  ");)
    return 0;
  }

// ================================================================================================
/// This function does one multigrid step
  double MGSolBase0::MGStep_B(
  int Level,            // Level
  double Eps1,          // Tolerance
  int MaxIter,          // n iterations
  const int Gamma,      // Control V W cycle
  const int Nc_pre,     // n pre-smoothing cycles
  const int Nc_coarse,  // n coarse cycles
  const int Nc_post,    // n post-smoothing cycles
  const int clearing    // clean the init matrix and precond
  ) {
    TRACKING_FUN(printf(" --> MGSolverBase0.C: MGSolBase0::MGStep \n  ");)
// ================================================================================================
    std::pair<int, double> rest(0, 0.);
    if(Level == 0) {   // coarse level
      rest = _solver[0]->solve(*A[0], *x[0], *b[0], Eps1, 200, clearing);// solver
      res[0]->resid(*b[0], *x[0], *A[0]); // coarse residual

      P_CONV(std::cout << "Coarse res " << rest.second << " " << rest.first << std::endl;)
      DEBUG_MGSOLB1(std::cout<< "\n coarse sol[" << Level <<"] " << *x[Level] << "\n  coarse res[" << Level <<"] " << *res[Level] << std::endl;)
    }
    else {  // fine levels
      TIME(std::clock_t start_time = std::clock();)
      // presmoothing (Nu1)
      int Nc_pre1 = Nc_pre;
      if(Level < _NoLevels - 1) Nc_pre1 *= 2;
      rest = _solver[Level]->solve(*A[Level], *x[Level], *b[Level], Eps1, Nc_pre1, clearing);
      res[Level]->resid(*b[Level], *x[Level], *A[Level]);  // presmoothing residual
      b[Level - 1]->matrix_mult(*res[Level], *Rst[Level - 1]); // restriction

      P_CONV(std::cout << " Pre Lev " << Level << " res " << rest.second << " " << rest.first;)
      TIME(std::clock_t end_time = std::clock(); std::cout << " time ="
           << double(end_time - start_time)/CLOCKS_PER_SEC << std::endl;)
      DEBUG_MGSOLB1(std::cout<< "\n  Pre-smooth x[" << Level <<"] " << *x[Level] <<
                    "\n  Pre-smooth (Level) b-Ax=e[" << Level <<"] " << *res[Level] <<
                    "\n  Restriction (on Level-1) Re=b[" << Level-1 <<"] " << *b[Level-1] << std::endl;)
      // --------- end presmoothing (Nc_pre)

      // -----------------------------------------------
      //  solving of system of equations for the residual on the coarser grid
      x[Level - 1]->close();
      x[Level - 1]->zero();
      double coarser_rest=0.;
      for(int g=1; g <= Gamma; g++) coarser_rest=MGStep_B(Level-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);
      // -------------------------------------------------------------------
      // interpolation of the solution from the coarser grid (projection)
      res[Level]->matrix_mult(*x[Level-1],*Prl[Level]);
      if(coarser_rest < rest.second) x[Level]->add(1., *res[Level]);   // add the coarser solution only if it helps

      //  res[L]= e[L-1] Prl[L] ->  e[L]= e[L]+res[L]  debug only
      DEBUG_MGSOLB1(std::cout<< "\n  (Level-1)  x[" << Level-1 <<"] " << *x[Level-1] << std::endl;
                    std::cout<< "\n   Proj (on Level) Pe[" << Level <<"] " << *res[Level] << std::endl;
                    std::cout<< "\n   Sol (on Level)  x+Pe[" << Level <<"] " << *x[Level] << std::endl;)

      // end projection
      // postsmoothing (Nc_post)
      TIME(start_time = std::clock();)

      int Nc_post1 = Nc_post;
      if(Level < _NoLevels - 1) Nc_post1 *= 2;
      rest = _solver[Level]->solve(*A[Level], *x[Level], *b[Level], Eps1, Nc_post1, clearing);
      res[Level]->resid(*b[Level], *x[Level], *A[Level]);//  postsmoothing residual


      // b[L]=res[L] from pre-smooth -> A[L]e[L]=b[L] ->  res[L]=-A[L]e[L]+b[L]  debug only
      P_CONV(std::cout << " Post Lev " << Level << " res " << rest.second << " " << rest.first;)
      TIME(end_time = std::clock(); std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;)
      DEBUG_MGSOLB1(std::cout<< "\n  post Level e[" << Level <<"] "<<   *x[Level] << std::endl;
                    std::cout<< "\n  post Level Ae-b=res[" << Level <<"] "<<   *res[Level] << std::endl;)
      // ----------------  end postsmoothing
    }
    // end cycle -------------------------------------
    res[Level]->close();
// ================================================================================================
    TRACKING_FUN(printf(" <-- MGSolverBase0.C: MGSolBase0::MGStep \n  ");)
    return rest.second;
  }

// ================================================================================================
/// Check for Prolong and Restr Operators
  void MGSolBase0::MGCheck_B(
  int Level
  ) const {
    TRACKING_FUN(printf(" --> MGSolverBase0.C: MGSolBase0::MGCheck \n  ");)
// ================================================================================================
    std::cout << "\nxlevel-1 before rest\n";  x[Level - 1]->print();
    std::cout << "\n x level before rest\n";  x[Level]->print();
    x[Level - 1]->matrix_mult(*x[Level], *Rst[Level - 1]);
    std::cout << "\nxlevel-1 after rest\n";  x[Level - 1]->print();
    x[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);
    std::cout << "\n x level after prol\n";  x[Level]->print();
    //   x[Level-1]->matrix_mult(*x[Level],*Rst[Level-1]);
    //   x[Level]  ->matrix_mult(*x[Level-1],*Prl[Level]);
    TRACKING_FUN(printf(" <-- MGSolverBase0.C: MGSolBase0::MGCheck \n  ");)
// ================================================================================================
    return;
  }





