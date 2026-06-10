#ifndef __mgsolbase0__
#define __mgsolbase0__


// c++ class include files ----------------------------------------------------
#include <string>
#include <vector>

// algebra  matrix-vector class files (contrib) -------------------------------
class SparseMatrixM;   // algebra sparse matrices
class SparseMMatrixM;  // algebra sparse quadrangolar matrices
class NumericVectorM;  // algebra numerical vector
class LinearSolverM;   // algebra linear solver

// FEMus  include files 
class MGMesh;

// ===============================================================================================
// define DEBUG_SOLVER only to DEBUG
//   DEBUG_MGSOLB1 to check the cycle MGSolBase0::MGStep_B  one multigrid step
//   DEBUG_MGSOLB2  to check the cycle  MGSolBase0::MGSolve V cycle
//   DEBUG_MGSOLB3  to check the cycle  MGSolBase0::MGSolve W cycle
  #define DEBUG_MGSOLB1(x)
//   #define DEBUG_MGSOLB1(x) x
#define DEBUG_MGSOLB2(x)
// #define DEBUG_MGSOLB2(x) x
// ===============================================================================================


// =======================================================
/// MGSolBase0 zero solver. It consists of:
/// - Matrix/vector dimension (Mesh data pointer and dof)
/// - Multigrid matrices and rhs's and  Multigrid operators
/// - labeling
/// - solver
// =======================================================
class MGSolBase0{
    public:
// const int MAX_DIMENSION=3;  // for element local only
//   const int MAX_NDOF_FEM=27;  // max -> LAGRANGE quad 27
//   const int MAX_NDOF_P=8;  // max -> LAGRANGE linear 8 
//    const int MAX_NDOF_FEMB=9; // // max -> LAGRANGE 9NDOF_P
//    const int MAX_NDOF_PB=4;
/// ***********************************************************************************************
///                         DATA class MGSolBase0
/// ***********************************************************************************************
/// // Matrix/vector dimension (Mesh data pointer and dof) -----------------------------------------
///    MGMesh&  _mgmesh;     ==  mesh pointer
///     const int _NoLevels;  ==  level number
///     int* _Dim;            ==  dimension number of dofs per level _Dim[level]=
///      int** _node_dof;      ///< _node_dof: Mesh[node] -> Matrix dof [unknown](built by init_dof())
///     // Parallel-multilevel solver     ---------------------------------------------------------
///     int _iproc;               == processor  (working processor)
///     LinearSolverM** _solver;  == linear system solver type (each level) *LinearSolverM[level]
    
/// //  Multigrid matrices and rhs's and  Multigrid operators -------------------------------------
///     std::vector<SparseMatrixM*>  A;   == A[level] (Matrix) _Dim[level]x_Dim[level] A=vector[_NoLevels]
///     std::vector<NumericVectorM*> b;   == b[level] (rhs) = 1x_Dim[level]
///     std::vector<NumericVectorM*> x;   == x[level] (solution)(parallel structures *LinearSolverM[level])
///     std::vector<NumericVectorM*> res; == r[level]=b-Ax
/// //   Multigrid operators    x(level+1)= P x(level); x[level-1]=R x[level]
///     std::vector<SparseMMatrixM*> Rst;  == PT(level->level-1)=R=P^T Restrictor (prologation transp)
///     std::vector<SparseMMatrixM*> Prl;  == P (level->level+1)= prologation
/// ***********************************************************************************************    
 public:
//  -------------------------------------------------------------------------------------------
// Matrix/vector dimension (Mesh data pointer and dof)
     // _mgmesh -> dof through _node_dof[_NoLevels][_Ndof_lev]
    MGMesh&   _mgmesh;    ///<  mesh pointer
    int       _Nnodes;    ///<  mesh number of nodes
    const int _NoLevels;  ///< level number (from _mgmesh._Nolevels)
    int*      _Ndof_lev;  ///< dimension number of dofs per level (from _mgmesh._Nolevels)
    int**     _node_dof;  ///< dof map global MGMesh[point]->Matrix[dof]  (built by init_dof())
    int  _type_solver=0; ///< type solver (default 0)           
     
protected:
//  -------------------------------------------------------------------------------------------
    // Parallel-multilevel solver
    int _iproc;               ///< processor
    LinearSolverM** _solver;  ///< linear system solver type (each level)
    bool _Solve_B;
//  -------------------------------------------------------------------------------------------
//  Multigrid matrices and rhs's and  Multigrid operators
    /// Matrix A  A[level] (Matrix) _Dim[level]x_Dim[level] A=vector[_NoLevels]
    std::vector<SparseMatrixM*>  A;  
    /// solution x  x[level] (Parallel numerical vector) 1x_Dim[level] x=vector[_NoLevels]
    std::vector<NumericVectorM*> x;    ///< solution x
    std::vector<NumericVectorM*> b;    ///< rhs b
    std::vector<NumericVectorM*> res;  ///< residual
//   Multigrid operators
    std::vector<SparseMMatrixM*> Prl;  ///< Prolongation
    std::vector<SparseMMatrixM*> Rst;  ///< Restrictor

/// *********************************************************************************************** 
/// ***********************************************************************************************
///                                 FUNCTIONS class MGSolBase0
/// ***********************************************************************************************
  // functions MGSolbase0  -> _B
   //  virtual void MGSol_init_B() =0;
   //  virtual int MGSolve_B
   //  virtual double MGStep_B 
   //  virtual void MGCheck_B 
   //  void Vanka_solve_B
   // virtual double MGStep_Vanka_B (

///     MGSolBase0 (MGMesh&)     MG mesh Constructor
///     -------------------------------    
///     ~MGSolBase0();  ///< Destructor (level structure)
///     -------------------------------    
///     void clear();  ///< Substructure destructor
///     -------------------------------    
///     virtual void MGSol_init() =0; // This function read the Operators (defined in the MGSolverDA)
///     ===================================
///     ///  PARALLEL MULTILEVEL SOLUTION  ( if no virtual =0  -> MGSolverBase0.C) --------------------
///     -------------------------------    
///     virtual void GenMatRhs ( const double,const intconst int) = 0; // Assemblying A matrix function
///     -------------------------------    
///     virtual void MGSolve (double Eps,int MaxIter, 0 or 1, 1,8,40,8);  /// This function performes the V,W and Fast MultiGrid scheme
///     -------------------------------    
///     virtual double MGStep (int,double,int,int,int,int,int,1); /// This function solves a MultiGrid step
///     -------------------------------    
///     virtual void MGCheck (int Level) const;   ///< Check Operators
///     ===============================    
///     Vanka - (MGSolverBase_VANKA.C) ---------------------------------------------------------------
///     -------------------------------    
///     void Vanka_solve (int,SparseMatrixM&, NumericVectorM&, NumericVectorM& rhs_in,double,int,false );
///     -------------------------------    
///     double Vanka_test (int Level);
///     -------------------------------    
///     virtual double MGStep_Vanka (int,double,int,int,int,int,int); ///< MultiGrid Step
/// ***********************************************************************************************    
public:
// ==================================================================================================
// CONSTRUCTOR-DESTRUCTOR
// ==================================================================================================
  
    
    MGSolBase0(MGMesh& mgmesh_in);  ///<  MG mesh
    ~MGSolBase0();                  ///< Destructor (level structure)
    void clear();                   ///< Substructure destructor
    
    void set_type_solver(int type_solver){_type_solver= type_solver;  } 
  
    virtual void MGSol_init_B() =0;   ///< This function read the Operators (defined in the MGSolverDA)
    
    // --------------------------------------------------------------------
    ///  PARALLEL MULTILEVEL SOLUTION  (MGSolverBase0.C)
    

    // --------------------------------------------------------------------
    /// This function performes the V,W and Fast MultiGrid scheme
     int MGSolve_B (
        double Eps,                ///< tolerance <-
        int MaxIter,               ///< n iterations <-
        const int clearing = 0,    ///< 1= clean the init matrix and precond <-
        const int Gamma = 1,       ///< Control V W cycle <-
        const int Nc_pre = 8,      ///< n pre-smoothing cycles <-
        const int Nc_coarse = 40,  ///< n coarse cycles <-
        const int Nc_post = 8      ///< n post-smoothing cycles <-
    );
    // ------------------------------------------------------------------
    /// This function solves a MultiGrid step
    double MGStep_B (
        int Level,            ///<  MG level <-
        double Eps1,          ///< tolerance <-
        int MaxIter,          ///< n iteratio <-
        const int Gamma,      ///< Control V W cycle  <-
        const int Nc_pre,     ///< n pre-smoothing cycles <-
        const int Nc_coarse,  ///< n coarse cycles   <-
        const int Nc_post,    ///< n post-smoothing cycles  <-
        const int clearing    ///< 1= clean the init matrix and precond <-
    );
    // --------------------------------------------------------------------
     void MGCheck_B (
        int Level  /// \param[in] <>   MG level
    ) const;   ///< Check Operators

// ************************************************************************ 
#ifdef VANKA
  // --------------------------------------------------------------------
    void Vanka_solve_B (           ///<  (MGSolverBase_VANKA.C)
        int Level,                                      ///< Level <-
        SparseMatrixM& matrix_in,                        ///< A matrix <-
        NumericVectorM& solution_in,                     ///< x solution <-
        NumericVectorM& rhs_in,                          ///< b rhs <-
        const double tol,                                ///< tolerance <-
        const int m_its,                                 ///< n cycles <-
        bool isdirect = false 
    );                        
    // --------------------------------------------------------------------
    double Vanka_test ( 
    int Level           ///< Level
    );
    // --------------------------------------------------------------------
     double MGStep_Vanka_B ( ///< MultiGrid Step
        int Level,              ///<  MG level <-
        double Eps1,            ///< tolerance <-
        int MaxIter,            ///< n iteratio  <-
        const int Gamma,        ///< Control V W cycle <-
        const int Nc_pre,       ///< n pre-smoothing cycles<-
        const int Nc_coarse,    ///< n coarse cycles <-
        const int Nc_post       ///< n post-smoothing cycles <-
                 );                       
#endif  // ****************************************************************   


};


#endif

