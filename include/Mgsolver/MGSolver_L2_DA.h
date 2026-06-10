#ifndef __mgsolverda_h__
#define __mgsolverda_h__

// conf includes --------------------------------------------------------------
// #include "MGFE_conf.h"
#include "Printinfo_conf.h"  //
#include "Solverlib_conf.h"

// algebra class --------------------------------------------------------------
#include "dense_matrixM.h"
#include "dense_vectorM.h"
class NumericVectorM;

// local class ----------------------------------------------------------------
#include "MGFE_L0_Map.h"
#include "MGSolver_L1_Base.h"  //for the inherited class
class MGUtils;
class MGSystem;
class MGMesh;
class MGEquationsMap;
class MGEquationsSystem;
class MGFEMap;
class MGMesh;

#include "hdf5.h"

// ================================================================================================
// ================================================================================================
/// This Class contains all the pointer table of the possible external fields for the system
class external_field
/// This class is based on enum FIELDS ( to define in Equation_conf.h)
/// FIELDS: equation symbol (*_F) -> preassigned order (0)
/// enum  FIELDS{
///    NS_F  =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSX_F =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSY_F =1,     // [1] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSZ_F =2,     // [2] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
///    .......................
///    ........................ see  Equation_conf.h
///  };
/// This class defines:
/// 1)  tab_eqs : FIELDS -> Eqs order (system appearance order)  (.)_F=0) ->  0
/// 2)  indx_ub : Eqs order -> data pos in  ub
/// 3)  ub      : I -> data
/// The  tab_eqs + indx_ub are defined in the
///                    MGSolverDA function   set_ext_fields(const std::vector<FIELDS> &pbName)
/// The vector ub is filled during assemblying by
/// MGSolverDA function get_el_sol(#block data,KLQ,el_ndof,el_conn,offset,starting,_data_eq[1].ub)
/// Each  MGSolverXX class has a _data  external_field in its data
// ================================================================================================
{// ===============================================================================================
public:
    ///@{ \name EXTERNAL FIELDS INFORMATION
    int n_eqs;                 ///< number of equations in the system
    const int max_neqs = 60;   ///< max number of equations (it can be changed)
    MGSolDA* mg_eqs[60];     ///< equation system pointer (mg_eqs)
    int tab_eqs[60];           ///< map external system and index
    int indx_ub[60];           ///< index of external equations(const/linear/quad):
    double ub[60*27];    ///< element external field old solution= max_neqs*27
    

    /// CONSTRUCTOR-DESTRUCTOR
    external_field() {}   //< Empty Constructor
    ~external_field() {}  //< Empty destructor
    ///@}
};
// =====================================================


// ================================================================================================
/// Class for MG Diffusion-Advection equation generic solvers
// ================================================================================================
class MGSolDA : public MGSolBase{
protected:
// ========================================================================
//                             DATA
// ========================================================================
//  -------------------------------------------------------------------------
// DATA POINTER
    MGEquationsSystem& _mgeqnmap_D;  ///<  equation map  pointerz
    MGUtils&           _mgutils_D;   ///<  utility class pointer
    MGFEMap&           _mgfemap_D;   ///<  FEM class pointer

public:
    // ************************************************************************************************
    //                   DATA
    // ************************************************************************************************
    
    // 4. geometry ---------------------------------------------------------------------------------------
    
    // Mesh
    int  _top_offset_D;  // top mesh nodes
    
    // element (volume)
    double _xx_qnds_D[27 *3];             ///< elem coords (volume) 
    double _xyz_D[3];                         ///< 1 point (middle)
    double _InvJac_D[3 *3];
    int  _el_conn_D[27];                        // map local->global (connectivity)
 
    // element (boundary)
    double _xxb_qnds_D[9 *3];          // boundary coordinates  (surface order)
    int  _elb_conn_D[9];                     // map local->global (boundary connectivity)
    int  _sur_toply_D[9];                    // boundary topology map surface-> local
    int  _el_sides_D;                                // element sides  (Boundary)
    int  _el_neigh_D[27];                       // bd element connectivity
    

    // 1. Common data MGSolDA -------------------------------------------------------------
    // Labeling
    const std::string _eqname_D;  ///< equation name
    //   Variables internal fields
  
    const int    _n_vars_D;         ///< total number of variables
    int          _varst_D[3];        ///< number of variables for each type [0],[1],[3]
    std::vector<std::string> _var_names_D;      ///< variable names
    std::vector<double>      _var_ref_D;        /// reference values (distance[0] mass[1] time[2])
    
    int          _dir_D;            ///< dir= eq index for vector system of equations
    int          _AxiSym_D;         ///< _AxiSym_D =0 no  _AxiSym_D =1 axis (x-axis); _AxiSym_D =2   fixed plane;
    int          _offset_level_flag_D[3];
    
  
    // 2. FEM/dof---------------------------------------------------------------------------------------
    // MGFE* _fe[3]  = fem  (piecewise linear,piecewise quadratic)
    /// int _el_dof[3]= number of dof in each element  (piecewise, linear, quadratic)
    //          _el_dof[0]=volume piecewice   _el_dof[1]=volume linear  _el_dof[2]=volume quad
    ///         _el_dof[3]=face piecewice    _el_dof[4]=face linear  _el_dof[5]=face quad
    int _el_dof_D[3]; // for each type variable vol dof 
    int _face_dof_D[3]; // for each type variable bd (face) dofs
    MGFEbase* _fe_D[5];    ///< fem  (piecewise (_fe_D[2]) linear(_fe_D[1]),piecewise quadratic(_fe_D[0]),...)
    MGFEbase* _fe_face_D[5];  ///< fem  (piecewise (_fe_D[2]) linear(_fe_D[1]),piecewise quadratic(_fe_D[0]),..)
    
    int _var_index_D[3]={0,1,2};
    int _var_index_Nonodes_D[3]={0,1,3};
    
    // 3. Local element storage MATRIX, VECTOR, shape functions ----------------------------------------
    // element local shape and derivative function at a  gaussian point
    double _xyz_g_D[3];                          ///< 1 gaussian point
    double _phi_g_D[3]  [27];                      ///< shape field (0-1-2 degree)
    double _dphi_g_D[3] [27*3];            ///< shape derivative  (0-1-2 degree)
    double _ddphi_g_D[3][27*3*3];  ///< shape second derivative  (0-1-2 degree)

 
    // 5. DenseMatrixM _KeM  = local  matrix;  DenseVectorM _FeM = local  rhs ----------------------------
    DenseMatrixM _KeM_D;  ///< local  matrix
    DenseVectorM _FeM_D;  ///< local  rhs
    int _mtr_nrows_D;     ///< local matrix dimension
   
    //6. External field to rhe MGSolverBase  (see section  below)

    external_field _DAdata_eq_D[3];  ///< external data structure
    int _FF_idx_D[60];      ///< field equation flag
    FIELDS _field_D;

    // ========================================================================
    // 1. Constructor destructor memory allocation
    // ========================================================================
    // const int nvars_in[],    = number of piecewise[0], linear[1], quadratic[2] variables
    //  std::string eq_name_in  = Equation name  ("DA")
    //  std::string varname_in  = Variable suffix name (u)
    //  const int Level,        = MG Level
    //  const int vb_0 = 0,     = mesh element fem type (vb_0=0 volume  (vb_0=1 surface)
    //  const int n_vb = 1      = mesh element groups (n_vb=1 only one fem type)
    // ========================================================================
    /// Level constructor  I
    
/// ----> *************************Constructor destructor *****************************************    
    MGSolDA (
        MGEquationsSystem& mg_equations_map,
        const int nvars_in[],           ///< number of piecewise[0], linear[1], quadratic[2] variables
        std::string eq_name_in = "DA",  ///< Equation name  ("DA")
        std::string varname_in = "u"    ///< Variable suffix name (u)
    );
    // Destructor------------------------------------------------------------------
    virtual ~MGSolDA();    ///< Destructor (level structure)
private:    
    void clean();  ///< Clean all substructures
/// <---- *************************Constructor destructor *****************************************    

public:
    virtual void SetValue(double value) {}

/// ----> *******************************  Init  MGSolDA solver ********************************************     
    ///< Contruction Setting  init_dof->init_fAx ->init_bc-> init_ic
    // This function read the Operators (defined in the MGSolverDA)
    void MGSol_init_B(); //
private:
///    ------> Setting operator dimensions r=f-Ax and Prl, Rst -----------------------------------------------
    void init_fAx (
        const int Level  ///<  init_fAx-> MULTILEVEL OPERATORS -> Reading/Writing MG operators
    );
    void init_PR(const int Level) ;///<  init_fAx-> MULTILEVEL OPERATORS  P R 
    
    virtual void init_ntgbg(){}
    
    ///  init_fAx-> MULTILEVEL OPERATORS -> Reading/Writing MG operators
    /// Reading/Writing M (Matrix) operators
    virtual void ReadMatrix (     ///< Reading Matrix
        const int Level,          // MG Level <-
        const std::string& name,  // file name for M <-
        SparseMatrixM& Mat, const int* nvars_in );
    /// Reading/Writing P(Prolongation) operators
    virtual void ReadProl (       ///< Read Prolongation Op
        const int Level,          // MG Level <-
        const std::string& name,  // file name for P <-
        SparseMMatrixM& Mat,      //  Matrix for P <-
        const int nvars_in[]); // , int node_dof_c[], int node_dof_f[] );
    /// Reading/Writing R (Restrictor) operators
    virtual void ReadRest (       ///< Restriction Op.
        const int Level,          // MG Level
        const std::string& name,  // file name (reading from)
        SparseMMatrixM& Rest,     // Restriction Matrix
        const int nvars_in[]     // # cost,linear quad variables
//         int node_dof_f[],         // dof map fine mesh
//         int node_dof_c[],         // dof map coarse
//         int _node_dof_top[]       // dof map top level
    );                            // ----
///   <------ Setting operator dimensions r=f-Ax and Prl, Rst -----------------------------------------------

///   ------>    system degrees of freedom 
///  This function initializes the system degrees of freedom (dof).  virtual void init_dof()
///  It Builds  _node_dof[Level]: map local-global (with multiple variable)
    void init_dof (
        const int Level,  ///< MG Level
        const int vb_0 = 0,
        const int n_vb = 1
    );
///   <------    system degrees of freedom 
    
///   ------> Setting boundary conditions -> init_bc
    ///   _bc[2];  ///< boundary conditions map (top level
    void init_bcic();    ///< Initial conditions   -> init_ic
    ///  init_bcic-> BOUNDARY/INITIAL CONDITIONS -> Reading/Writing BC/IC
    /// Reading BC
    virtual void GenBc_loop (
        const int vb, const int ndof_femv, const int n_ub_dofs, const int n_pb_dofs, const int n_kb_dofs,
        int bc_id[], int mat_id[] );
    /// Reading IC
    virtual void GenIc_loop (
        const int vb, const int ndof_femv, const int n_ub_dofs, const int n_pb_dofs, const int n_kb_dofs,
        int bc_id[], int mat_id[] );


///    <------ Setting boundary conditions -> init_bc
/// <----- *******************************  Init  MGSolDA solver ********************************************     
 
public:
    // ========================================================================
    //    MultiGrid Solver (defined in the MGSolver__)
    // ========================================================================
       ///  MULTILEVEL SOLUTION/ASSEMBLYING  (MGSolverBase.C)
    virtual int MGTimeStep (
        const double time,  /// \param[in] <>  time
        const int i_step, 
        const int i_iter,  
        const int mode     /// \param[in] <>   rhs assembler flag
    );                      ///< MG time step solver (backward Euler)
    void read_u_D(std::string namefile, std::string var_name, std::vector<double> &target_sol);
protected:    
    // ---------------------------------------------------------------------------------------------------------
    // This function assembles the matrix   (defined in the MGSolver__)
    virtual void GenMatRhs ( ///< Volume Assemblying  matrix-rhs
        const double time,   // time
        const int Lev,       // Level
        const int mode          // rhs assembly control
    );
    virtual double CalcFUpwind (double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes );

public:
    // ========================================================================
    /// RETURN FUNCTIONS
    // ========================================================================
    // **********************************************************************************************
    //  2) External field to rhe MGSolverBase  (see section  below)
    // **********************************************************************************************
    // Functions ===================================================================================
    //   virtual void set_ext_fields(const std::vector<FIELDS> & pbName)=0;   ///< set external fields
    virtual void setUpExtFieldData();
    // printf("\n \n Wrong use of function print_ext_data in SolverBase for coupled mesh \n \n") =0*/
    virtual void print_ext_data ( double /*vect_data*/[] ) {};
    // ------------------------------------------------------------------------
    void ActivateVectField ( int Order, int Field, std::string SystemFieldName, int& n_index, int coupled );
    void ActivateScalar ( int Order, int Field, std::string SystemFieldName, int& n_index );
    void ActivateControl ( int Order, int Field, std::string SystemFieldName, int& n_index, int vector, int neq );
    void ActivateCoupled (
        int Order, int Field, std::string SystemFieldName, int& n_index, std::string SystemFieldName2 );
    void ActivateEquation ( int Order, int Field, std::string SystemFieldName, int& n_index );
    void ActivateDA (
        int Order,    ///< max order to storage (quad>lin>konst)
        int Field,        ///< Field to activate
        std::string sfn,  ///< SystemFieldName
        int& n_index     ///< n_index (collecting index)
    );


// =================================================------
// set functions
// // =================================================------
// //  Get functions
    // -------------------------------------------------------------------------------
    /// Return element dof indices using local-to-global map
    void get_el_dof_indices (
        const int Level,      ///< \param[in] <Level>   level
        const int iel,        ///< \param[in] <iel>     eLement number
        const int el_conn[],  ///< \param[in] <el_conn> connectivity
        const int el_dof[],   ///< \param[in] <el_dof>  quadratic[2] linear[1] const[0] dofs
        const int offset,     ///< \param[in] <offset>  offset for connectivity
        std::map<int, std::vector<int>>& el_dof_indices,  ///< \param[out]<el_dof-indices>  dof indices
        int  nvars[]
    ) const;
    //  ------------------------------------------------------
    /// This function interpolates a vector field over the fem element
    void interp_el (
        const double uold[],    // node values <-
        const int ivar0,        // init variable  <-
        const int nvars,        // # of variables  <-
        const double phi[],     // shape functions  <-
        const int n_shape,      // # of shape functions  <-
        double u_int[],         // interpolated function ->
        int dim,                // deriv or
        const int sur_tpgly[],  // surface nodes topology <-
        const int el_ndof       // surface nodes topology <-
    ) const;
    // -------------------------------------------------------------------------
    // This function interpolates the solution at gaussian point (phi[])
    void interp_el_sol (
        const double uold_b[],  // node values <-
        const int ivar0,        // init variable  <-
        const int nvars,        // # of variables  <-
        const double phi[],     // shape functions  <-
        const int n_shape,      // # of shape functions  <-
        double uold[]           // interpolated function ->
    ) const;
    // -------------------------------------------------------------------------
    // This function interpolates the 1derivative  at gaussian point (dphi[])
    void interp_el_gdx (
        double uold_b[],      // node values <-
        const int ivar0,      // init variable  <-
        const int nvars,      // # of variables  <-
        const double dphi[],  // derivatives of the shape functions  <-
        const int n_shape,    // # of shape functions  <-
        double uold_dx[]      // interpolated derivatives ->
    ) const;
    // -------------------------------------------------------------------------
    // This function interpolates the  2derivative at gaussian point (dphi[])
    void interp_el_gddx (
        double uold_b[],      // node values <-
        const int ivar0,      // init variable  <-
        const int nvars,      // # of variables  <-
        const double dphi[],  // derivatives of the shape functions  <-
        const int n_shape,    // # of shape functions  <-
        double uold_dx[]      // interpolated derivatives ->
    ) const;
    // -------------------------------------------------------------------------
    // This function interpolates the  1derivative
    // on the boundary (-> sur_tpgly[]) at gaussian point (dphi[])
    void interp_el_bd_gdx (
        const double uold_b[],  // node values <-
        const int sur_tpgly[],  // surface nodes topology <-
        const int el_ndof,      // surface nodes topology <-
        const int ivar0,        // init variable  <-
        const int nvars,        // # of variables  <-
        const double dphi[],    // derivatives of the shape functions  <-
        const int n_shape,      // # of shape functions  <-
        double uold_dx[]        // interpolated derivatives ->
    ) const;                // =======================================
    // -------------------------------------------------------------------------
    // This function interpolates the  solution
    // on the boundary (-> sur_tpgly[]) at gaussian point (dphi[])
    void interp_el_bd_sol (
        const double uold_b[],  // node values <-
        const int sur_tpgly[],  // surface nodes topology <-
        const int el_ndof,      // surface nodes topology <-
        const int ivar0,        // init variable  <-
        const int nvars,        // # of variables  <-
        const double phi[],     // shape functions  <-
        const int n_shape,      // # of shape functions  <-
        double uold[]           // interpolated function ->
    ) const;                // =======================================
    // -------------------------------------------------------------------------

    // -----------------------------------------------------------------------
    /// Dof , the bc and the solution  vector at the nodes of  an element
    void get_el (
        const int Level,                   // level <-
        const int nvar0,                   // intial varables number based on 27 <-
        const int nvar,                    // final varables number based on 27 <-
        const int el_nds,                  // number of nodes for element
        const int el_conn[],               // connectivity <-
        const int offset,                  // offset <-
        std::vector<int>& el_dof_indices,  // DOF indices ->
        int bc_dofs[][27],           // boudary conditions ->
        double uold[]                      // solution ->
    ) const;
    // -----------------------------------------------------------------------
    void get_el_dof_bc (
        const int Level,  // level
        const int iel,
        //   const int nvars[],       // # of variables to get  <-
//         const int el_nds[],                // # of element nodes for this variable  <-
        const int el_conn[],               // connectivity <-
        const int offset,                  // offset for connectivity <-
        std::vector<int>& el_dof_indices  // element connectivity ->
//         int bc_vol[],                      // element boundary cond flags ->
        // int bc_bd[]                        // element boundary cond flags ->
    ) const;
  // ----------------------------------------------------------------------  
      void get_el_dof_bc (
        const int Level,  // level
        const int iel,
        //   const int nvars[],       // # of variables to get  <-
//         const int el_nds[],                // # of element nodes for this variable  <-
        const int el_conn[],               // connectivity <-
        const int offset,                  // offset for connectivity <-
        std::vector<int>& el_dof_indices,  // element connectivity ->
//         int bc_vol[],                      // element boundary cond flags ->
         int bc_bd[]                        // element boundary cond flags ->
    ) const;
    
    

    // ========================================================================
    /// SET FUNCTIONS
    // ========================================================================
    //    virtual void set_cp_vector(const int&, const int&);
    // -----------------------------------------------------------------------
    virtual void set_xooold2x();
    // -----------------------------------------------------------------------
    void set_el_dof_bc (
        const int Level,  // level
        const int iel,
        //   const int nvars[],       // # of variables to get  <-
//         const int el_nds[],                // # of element nodes for this variable  <-
        const int el_conn[],               // connectivity <-
        const int offset,                  // offset for connectivity <-
        std::vector<int>& el_dof_indices,  // element connectivity ->
        int bc_vol[],                      // element boundary cond flags ->
        int bc_bd[]                        // element boundary cond flags ->
    );
    // -----------------------------------------------------------------------

    void set_dt ( double dt )  {   _dt_B = dt;  };  // computation val3=musker

    // ========================================================================
    /// COMPUTE FUNCTIONS
    // ========================================================================
    // -------------------------------------------------------------------------
    void compute_jac (
        const int j, const int idim,
        double uold_b[],      // node values <-
        const int nvars,      // # of variables  <-
        const double phi[],   // shape functions  <-
        const double dphi[],  // derivatives of the shape functions  <-
        const int n_shape,    // # of shape functions  <-
        double u_forw[],      // interpolated function ->
        double u_back[],      // interpolated function ->
        double u_forw_dx[],   // interpolated derivatives ->
        double u_back_dx[]    // interpolated derivatives ->
    ) const;

    // // -------------------------------------------------------------------------
    virtual double MGFunctional_B ( double, double& );
    virtual double MGFunctional_B (){ return 0.0;}

    // ========================================================================
    ///  EXTERNAL FIELDS
    // ========================================================================
    //    int Order                       = pol interpolation order
    //    int Field,                      = Field
    //    std::string SystemFieldName,    = Field name
    //    std::string SystemFieldName2,   = Coupled Field name
    //    int& n_index,                   = equation index
    //    int coupled                     = coupled (1) or segregated (0)  solver
    //    int vector                      = vector
    //    int neq                         = number of equations
    //  =========================================================================================
    // ------------------------------------------------------------------------
    //  virtual
    // auxiliary function fro user purpose
    virtual double eval_var1 ( double[] )    {        return 0.;    }
    virtual double eval_var2 ( double[] )    {        return 0.;    }
    virtual double eval_var3 ( double[] )    {        return 0.;    }

    // ========================================================================
    /// SOLUTION WRITE AND READ
    // ========================================================================
 

    // -------------------------------------------------------------------------
    // Read
    virtual void read_u_D (  /// Read solution from a xdmf file
        std::string name,  // filename <-
        int Level          // restart level <-
    );
    // -------------------------------------------------------------------------
    /// This function  defines the boundary conditions for the system:
    virtual void bc_intern_read_D (
        int /*face_id_node*/,  ///<  face identity           (in)
        int /*mat_flag*/,      ///<  volume identity         (in)
        double /*xp*/[],       ///< xp[] node coordinates    (in)
        int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
        int bc_flag[]          ///< boundary condition flag  (out)
    );
    // -------------------------------------------------------------------------
    /// This function reads Boundary conditions  from function
    virtual void bc_read_D (
        int face_id_node,  ///<  face identity          (in)
        int mat_flag,      ///<  volume identity         (in)
        double *xp,       ///< xp[] node coordinates    (in)
        int *bc_Neum,         ///< Neuman (1)/Dirichlet(0)  (out)
        int *bc_flag          ///< boundary condition flag  (out)
    ) {std::cout << "not defined" ; abort();}
    // -------------------------------------------------------------------------
    virtual void ic_read_D ( /// initial conditions  from function
        int k,             // bc from gambit    <-
        int m,             // material from gambit    <-
        double *xp,       // point coordinates <-
        int iel,           // element           <-
        double *value     //  point values     ->
    ) {std::cout << "not defined" ; abort();}
    
   // Print bc  =========================================================
    virtual void print_bc_D (    ///< Print boundary conditions to a xdmf file.
        std::string namefile,  // filename <-
        const int Level        // MGLevel  <-
    );
    // Print   =========================================================
    void print_u_D (             /// Print solution to a xdmf file.
        std::string namefile,  // filename <-
        const int Level        // MGLevel  <-
    );

    // Print solution in hdf5 =========================================================
    virtual void print_u_hdf5_D ( std::string namefile, const int Level );
    /// This function prints the solution for quad in hdf5
    virtual void print_u_hdf5_quad_D (
        const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[] );
     virtual void print_u_hdf5_quad_D (
        const  int ivar, const  std::string var_name,
        const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[] );
    /// This function prints the solution for  linear fem in  hdf5
    virtual void print_u_hdf5_lin_D (
        const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[] );
    virtual void print_u_hdf5_lin_D (
        const  int ivar, const  std::string var_name,
        const int Level, const int istep, const int offset, hid_t file_id,
        hsize_t dimsf[] );
    /// This function prints the solution: for quad and linear fem
    virtual void print_u_hdf5_const_D (
        const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[] );
    // =================================================================================
    
    /// Print xml attrib ================================================================
    virtual void print_u_xdmf_D (
        std::ofstream& out,  // file stream to print ->
        int nodes,
        int nelems,
        std::string file_name
    );
    void print_u_xdmf_pts_D (
        std::ofstream& out,     ///<  xdmf file//  file xdmf
        int dimension,          ///<  number of nodes
        std::string var_name,   ///<  var names
        std::string file_name,  ///<  hdf5 file name
        std::string datatype="Node",   ///<  datatype=Float, Int
        std::string center="Float"      ///<  center= Cell, Nodes
    ) const ;  // =========================================================================   
    
    

    // ==================== check =======================
    void check_multiphysics_D ( int mode );
    void check_mtrx_rhs_D ( int mode );
    void check_ass_nan_inf_D ( int mode );
    void check_bc_D ( int mode );

    // ************************************************************************
    //                                  MED
    // ************************************************************************
#ifdef HAVE_MED
    // Data ========================================================================================
    // Functions ===================================================================================
    MEDCoupling::MEDCouplingFieldDouble* _ExtField;
    std::map<int,int> _map_bc2g;
    // ------------------------------------------------------------------------
    // --------------------------------------------------------------------
    virtual void print_u_med_D ( /// Print solution to a med file.
        std::string namefile,  // filename <-
        const int Level ) ; // MGLevel  <-
    // --------------------------------------------------------------------
    virtual void print_weight_med_D ( /// Print weight for control problems to a med file.
        std::string namefile,       // filename <-
        const int Level );      // MGLevel  <-
#endif


// Two-phase coupling =============================================================================
#ifdef TWO_PHASE_LIB
    MGSolCC* _msolcc;
    void set_mgcc ( MGSolCC& cc );
#endif



};
#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
