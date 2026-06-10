#ifndef __mgsolverrt_h__
#define __mgsolverrt_h__


// ================================================================================================
// --------------   Vector system with rotation R -------------------------------------------------
// ================================================================================================
// ************************************************************************************************
// ================================================================================================
//    Constructor (MGSolverR_0.C)-> MGSolR::MGSolR
//    Destructor                 -> MGSolR::~MGSolR
//    Init
// ================================================================================================
//   GET/SET(Virtual)  ->  G1) get_el_vel
//                         G2)  get_el_nl_vel
//                         G3)  get_el_press
//    SET
// ================================================================================================
// ================================================================================================
// TIME STEP (MGSolverR_0.C)->  T1) MGSolR::MGTimeStep       MGSolverR_0.C
//                              T2) MGSolR::MGTimeStep_no_up
//                              T3) MGSolR::MGTimeStep_solve
//                              T4) matrix_RART_R (Virtual)
// ================================================================================================
// BCuniT (MGSolverR_BIC.C) ->  B1) MGSolR::set_unit_vectors_along_stream_
//                              B2) MGSolR::set_test_fun_bound_normals_R
//                              B3) MGSolR::set_unit_vectors_on_boundary_R
// ================================================================================================
//  Check (MGSolverR_CHECK.C)->  C1) MGSolR::check_file_R
//                               C2) MGSolR::check_KeMBC_R
//                               C3) MGSolR::check_FeMBC_R
//                               C4) MGSolR::check_res_R
// ================================================================================================
//  Virtual not here
// ================================================================================================
//  Print u  (virtual) -> P1) print_u_hdf5_D
//                        P2) print_u_xdmf_D
//  Read bc (virtual) ->  R1) bc_read_D
//  Read ic (virtual) ->  R2) ic_read_D
// ================================================================================================
 
 
 

// classe include ---------
#include "MGSolver_L2_DA.h"
// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;

// debug normal-tan-bitan structure
#define DEBUG_1(x) 


// ===========================================================================
//                                Navier-Stokes equation class
//=============================================================================
/// Navier-Stokes equation  (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolR: public MGSolDA {
/// Class for mg Navier-Stokes equation  with name NS_EQUATIONS.
/// Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)



protected:
  // Element Volume data ---------------------------------------------------------

  // a:) basic element topology for computation
  double _xx_qnds_R [27 * 3];  /// element point coord vector
  double _xxb_qnds_R[9 *3]; ///<  el-point coords on boundary

  // ======================  local matrices ===================================
  double _InvJac2_R [3* 3]; ///< Quadratic Jac at a gaussian point  (nabla x(psi))
  double _InvJac1_R [3 * 3]; ///< Linear Jac at a gaussian point (nabla x(psi))

  double * _FeMBC_R;     ///< local rhs
  double * _KeMBC_R;     ///< local matrix

  int _el_mat_nrows_R ;  ///< local dofs in the matrix (matrix dimension)
  int _mtr_aux_nrows_R;  ///< global matrix dimension
  double _matrix_R[27][3][3]; // local rotation matrix

public:


// ============================================================================================
// ============================================================================================
// CON/DESTRUCTORS: these functions are defined in MGSolverR_0.C or online here
// ============================================================================================
//                 MGSolverR_0.C
// ================================================================================================
//     Constructor -> MGSolR::MGSolR
//     Destructor  -> MGSolR::~MGSolR
//     Init   
// ============================================================================================
    
// ------------------------------------------------------------------------------------------------   
MGSolR(               ///< Constructor
  MGEquationsSystem & mg_equations_map, ///< equation map class (Mesh and parameters)
  int nvars_in[],                       ///< KLQ number of variables
  std::string eqname_in = "R0",         ///< base name system
  std::string varname_in = "u"          ///< base name variable
  );
// ------------------------------------------------------------------------------------------------
 ~MGSolR();            ///< Destructor
  
  
  
  
  
  
  
// ================================================================================================
//              Time step functions (solution)
// ================================================================================================
// ================================================================================================
// TIME STEP (MGSolverR_0.C)->  T1) MGSolR::MGTimeStep     
//                              T2) MGSolR::MGTimeStep_no_up
//                              T3) MGSolR::MGTimeStep_solve
//                              T4) MGUpdateStep
// Time step function:
//  mode  = 0 only  MGTimeStep_no_up
//  mode  = 1 only  MGUpdateStep
//  mode  = 2 both time step + update
// ================================================================================================    
  virtual    int MGTimeStep(
  double const time,  ///< time
  int const  i_step,  ///< time step
  int const  n_iter,  ///< n iterations
  int const  mode     ///< mode = 0 only  MGTimeStep_no_up 1 only  MGUpdateStep 2 both
  ) ;

// ------------------------------------------------------------------------------------------------
// This function does the time step but no update
  virtual   int  MGTimeStep_no_up(
  double const time,             ///< time
  int const  delta_t_step_in,    ///< number of time steps
  int const  n_iter              ///< number max of iterations
  );
// ------------------------------------------------------------------------------------------------
/// This function is the solver routine (called by MGTimeStep...)
  virtual int MGTimeStep_solve(
  const double time,     ///< time
  const int    max_iter  ///< number max of iterations
  );
// ------------------------------------------------------------------------------------------------
/// This function does the  update
  virtual   void  MGUpdateStep()=0;
  
void init_ntgbg();
protected:  
// ================================================================================================
//                     Oxy (stream, normal , tang)
// ================================================================================================
// ================================================================================================
// BCuniT (MGSolverR_BIC.C) ->  B1) MGSolverR::set_unit_vectors_along_stream_R
//                              B2) MGSolverR::set_test_fun_bound_normals_R
//                              B3) MGSolverR::set_unit_vectors_on_boundary_R
//                              B4) MGSolverR::set_unit_vectors_along_stream_R
// ------------------------------------------------------------------------------------------------
// This function sets the unit coords  system  with geometry construction
  void set_unit_vectors_on_boundary_R(
  const int option,         ///< tengent option
  const int itime_flag,     ///< computes only on time=itime_flag
  int coupled               ///< segrgate (0) or coupled normal storage
  );
// ------------------------------------------------------------------------------------------------
// This function sets the unit coords  system  with geometry construction
  void set_test_fun_bound_normals_R(
  double normi_tmp[],  ///< basic boundary normal
  int coupled          ///< segrgate (0) or coupled normal storage
  );
  // ------------------------------------------------------------------------------------------
  // This function sets the unit coords  system with n defineed by test function
  void set_unit_vectors_on_boundary_R(
  double  normi_tmp[],    ///< basic boundary normal
  const  int option,      ///< options
  const  int itime_flag,  ///< time  integer
  int coupled             ///< coupled (1) or segregated (0) storage
  );
  // -----------------------------------------------------------------------------------------
  // This function sets the unit coords  system along stream unit
  // -----------------------------------------------------------------------------------------
  void  set_unit_vectors_along_stream_R(
  const int i_step,   ///< time step
  const double vel[], ///< unit velocity
  double stream[],    ///< unit stream vector
  double normal[],    ///< unit normal
  double tg2[]        ///< second unit vector on tangent plane (3D)
  );
  
  
// ============================================================================================
// SET/GET FUNCTIONS: these functions are defined in inherited classes
// ============================================================================================
//   GET/SET(Virtual)  ->  G1) get_el_vel
//                         G2)  get_el_nl_vel
//                         G3)  get_el_press
//    SET
// ================================================================================================
  // ---------------------------------------------------------------------------------------------
  /// This function gets the old velocity field.
  /// It is implemented in the son classes (coupled and projection solvers)
  // --------------------------------------------------------------------------------------------
  virtual   void get_el_vel(
//   const int Level,    
//         const int el_ndof[], ///< element dof for quad,lin,konst dofs
  const int el_conn[], ///< element node connectivity (local-> global)
  double u_xyz[]       ///< element velocity
  ) =0;/*{
        std::cout << "  get_el_vel Virtual: not implemented in MGSolver class  \n"; abort();
    }*/
  // ---------------------------------------------------------------------------------------------
  /// This function gets the old non linear velocity field.
  /// It is implemented in the son classes (coupled and projection solvers)
  // --------------------------------------------------------------------------------------------   
  virtual   void    get_el_nl_vel(
      const int Level,
        const int el_conn[], ///< element node connectivity (local-> global)
  double u_nlxyz[]      ///< element velocity
      ) =0;
  // --------------------------------------------------------------------------------------------
  // This function get the old pressure field
  // It is a virtual function implemented in the son classes (coupled and projection solvers)
  // --------------------------------------------------------------------------------------------
  virtual  void get_el_press(
//         const int el_ndof[],   ///< element dof for quad,lin,konst dofs
  const int Level,  
  const int el_conn[],   ///< element node connectivity (local-> global)
  double p_xyz[]     ///< element pressure
  )=0;
  // --------------------------------------------------------------------------------------------
  // This function rotates the local matrix
  // It is a virtual function implemented in the son classes (coupled and projection solvers)
  // --------------------------------------------------------------------------------------------
  virtual  void matrix_RART_R(
  double FeMBC[],                         ///< rhs with no rotation (<-)
  double KeMB[],                          ///< matrix with no rotation  (<-)
  double matrix_R[][3][3] ///< rotation matrix  (->)
  )=0;                    

public:
// ==============================================================================================
// ==============================================================================================
//  Input/Output: these functions are defined in MGolverNS0_IOH.C   or in inherited classes
// ==============================================================================================

    
// ================================================================================================
//  USER:(Virtual in MGSolver_R) these functions are defined in USER_NS.C or in inherited classes
// ================================================================================================
//  Read bc (virtual) ->  R1) bc_read_D
//  Read ic (virtual) ->  R2) ic_read_D   
// ==============================================================================================
/// This function  reads boundary conditions
  virtual    void ic_read_D(
      int bc_gam, int bc_mat, ///< boundary and material type 
      double xp[],            ///< point coordinates
      int iel,                ///< element
      double u_value[]        ///< point field values
                         )=0;
// --------------------------------------------------------------------------------------------
/// This function reads initial conditionsc
 virtual    void bc_read_D(
     int bc_gam, int bc_mat,  ///< boundary and material type 
     double x[],              ///< point coordinates
     int bc_Neu[],            ///<  bc volume integral flag
     int bc_bd[]              ///<  bc surface integral flag
 )=0;


// ============================================================
//  PRINT SOL (Virtual in MGSolver_R)
// ============================================================
//  Print u (virtual) -> P1) print_u_hdf5_D
//                       P2) print_u_xdmf_D                         
// ============================================================
/// Print solution in hdf5 format
  virtual  void print_u_hdf5_D(
  std::string namefile, ///<  hdf5 file name
  const int Level       ///<  level (multi-level)
  )=0;
  /// Print xdmf attrib
  virtual void print_u_xdmf_D(
  std::ofstream& out,   ///<  xdmf file
  int nodes,            ///<  number of nodes
  int nelems,           ///<  number of elements
  std::string file_name ///<  xdmf file name
  )=0;  // ================================      */

  
 
  
  
  
  
  
  
  
  
  
protected:
  
  // This function prints tan-normal-bitan system
  void print_normal_hdf5_R(
    const int nodes,                  ///< nodes
    const hid_t file_id,             ///< file name to print
    hsize_t  dimsf[],
    double sol[]
); 

// This function  prints only tan-normal-bitan system
void print_normal_xdmf_R(
    std::ofstream& out,   ///<  xdmf file
    int nodes,            ///<  number of nodes
    std::string file_name ///<  xdmf file name
) ;
  
  private:
  
// ================================================================================================
//  Check (MGSolverR_CHECK.C)->  C1) MGSolR::check_file_R
//                               C2) MGSolR::check_KeMBC_R
//                               C3) MGSolR::check_FeMBC_R
//                               C4) MGSolR::check_res_R
// ================================================================================================
  int check_file_R(
  const std::string name_file
  );
  int check_KeMBC_R(
  const std::string name_file,
  const std::string name_dataset
  );

  int check_FeMBC_R(
  const std::string name_file,
  const std::string name_dataset
  );

  int check_res_R(double res[]
               );
  int check_KeMBC1_R(
  const std::string name_file,
  const std::string name_dataset
  );

  int check_FeMBC1_R(
  const std::string name_file,
  const std::string name_dataset
  );
  
void CHECK_R1();
};


#endif // endif _mgsolverr_h















