#ifndef __mgSolBase__
#define __mgSolBase__


class MGEquationsSystem;
#ifdef HAVE_MED  // ---------------------------------------------------------
#include "MEDCouplingFieldDouble.hxx"
#endif
#ifdef TWO_PHASE_LIB
class MGSolCC;

#endif
#include "MGSolver_L0_Base.h"
#include <memory>
// ===============================================================================================
/// MGSolBase is a basic solver adding  specialized structures for
/// 1) time dependent solvers
/// 2) non linear solvers
/// 3) control solvers
// ===============================================================================================
  // functions MGSolbase0  -> _B ----------------------------------------
   //  virtual void MGSol_init_B() =0;
   //  virtual int MGSolve_B
   //  virtual double MGStep_B 
   //  virtual void MGCheck_B 
   // void Vanka_solve_B
   // virtual double MGStep_Vanka_B 
  // functions MGSolbase  -> _B -----------------------------------------
  // virtual void print_u_xdmf_B 
  // virtual double MGFunctional_B 
  //   virtual void  set_ctrl_dom_B 



class MGSolBase: public MGSolBase0
{


public:

// ***********************************************************************************************
//                                        DATA
// ***********************************************************************************************
//  ----------------------------------------------------------------------------------------
    //  Old VECTOR see SOLUTION SECTION  MGSolBase_SOL.C  (see section  below)
    //   std::vector<NumericVectorM*> x_old[3];  ///< Multi level [l=0..NoL] old solution x old x[0] oold x[1]
    //   ooold x[2] std::vector<NumericVectorM*> x_nonl;    ///<  Multi level [l=0..NoL] non linear solution x
    //   std::vector<NumericVectorM*> x_aux;     ///< vector for multiple uses on top
    //   std::vector<NumericVectorM*> d_aux;     ///< vector for multiple uses
    //   std::vector<double> _weight_ctrl;      ///< controlled region for optimal control problems
// _NumRestartSol= Number of solutions needed for restart
    // int  *_bc[2];  ///< boundary conditions map (top level) see section below
    // _bc[0]=volume
    //                 This number is used for printing and reading solutions for system restart.
    //                 Default value is equal to 1. With value _NumRestartSol=2 the solution
    //                 stored in _x_oold vector is printed with _old suffix, while with _NumRestartSol=3
    //             also the _x_ooold stored solution is printed, with _oold suffix (visible inside .h5 file)

// 1) Data  time dependent solver
/// ***********************************************************************************************
     // static constexpr int MDIM_FEM=3; // max geom dim FEM (volume(3))
     //  static constexpr int NDOF_FEM=27; // max dof FEM (volume(3) LAGRANGE quad=27)
     //    static constexpr int NDOF_FEMB=9; // max dof FEM (volume(2) LAGRANGE quad=9)
     //   static constexpr int MAX_NDOF_P=8;
    
    
    int _ndim_B;       ///<  = _mgmesh.dimension
    double _dt_B;  ///< time step
    int _NumRestartSol;  ///< Number of solutions needed for restart

    std::vector<std::unique_ptr<NumericVectorM>> x_old[3];  ///< Multi level [l=0..NoL] old solution x old x[0] oold x[1]  ooold x[2]
    std::vector<std::unique_ptr<NumericVectorM>> x_aux;     ///< vector for multiple uses on top mesh
    std::array<std::vector<int>, 2> _bc;  ///< boundary conditions map (top level)see MGSolDA_BCIC.C

// 2) Data Non linear solvers
    std::vector<std::unique_ptr<NumericVectorM>> x_nonl;  ///<  Multi level [l=0..NoL] non linear solution x

// 3) Data control  solvers
    std::vector<std::unique_ptr<NumericVectorM>> d_aux;   ///< vector for multiple uses on top mesh
    double _control;                   ///< control flag
    std::vector<std::vector<double>> _weight_ctrl;  ///< controlled region for optimal control
    std::vector<double> _weight_controlled;  ///< controlled region for optimal control
    std::vector<double> _weight_lift;  ///< lift region for optimal control

// ***********************************************************************************************
// CONSTRUCTOR-DESTRUCTOR
// ***********************************************************************************************
    MGSolBase ( MGMesh& mgmesh_in );                 ///<  MG mesh
    virtual ~MGSolBase();                                    ///< Destructor (level structure)
    void clear();                                     ///< Substructure destructor
  static std::unique_ptr< MGSolBase> build ( MGMesh& mgmesh_in, FIELDS field_type,
            MGEquationsSystem& mg_equations_map_in,  //
            const int nvars_in[],                    // # of quad variables
            std::string eqname_in,       // equation name
            std::string /*varname_in*/);

// ***********************************************************************************************
//  1) TIME DEPENDENT SOLVERS
// ***********************************************************************************************

    // set/get parameters from solvers
    virtual void   set_dt ( double dt ); ///< MG time step solver (backward Euler)
    virtual void   SetValue ( double value ) = 0; ///< This function sets a value in the solver solution
    virtual double GetValue ( int flag ); ///< This function returns a value from the solvers
    virtual void   SetValueVector ( std::vector<double> value ); ///< This function sets a value in the solvers
  

   // -----------------------------------------------
   // set functions
//     virtual void set_xooold2x() = 0;
    void localize_xooold();
    void set_sol ( int i, int kdofs, double value );   ///< set x_old[level](kdofs)=value
    void set_x_aux ( int i, int kdofs, double value ); ///< set x_aux[i](kdofs)=value at toplev

    // -----------------------------------------------
 
// -----------------------------------------------
// get functions
    double get_sol ( int i,int kdofs ); ///< This function gets the old_sol[0,1,2](kdofs) over all mesh levels
    double get_x_aux ( int i,int kdofs ); ///< This function gets the x_aux(kdofs) on the top mesh
    
    
    void get_el_sol_F(
    const int level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_faces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    
    ) const; 
    
void get_el_sol_C(
    const int level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_faces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    
    ) const;     
    
    
    void get_el_nonl_sol_F(
     const int Level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_faces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    
    ) const; 
    
void get_el_nonl_sol_C(
      const int Level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_faces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    
    ) const;     
    
    
 
    // ----------------------------------------------------------------
    void get_el_sol (
        const int i_step,     ///<  i-step
        const int ivar0,      ///< \param[in]   <ivar0>   initial variable
        const int nvars,      ///< \param[in]   <nvars>   number of variables to get  <-
        const int el_nds,     ///< \param[in]   <el_nds>  number  of element nodes for this variable
        const int el_conn[],  ///< \param[in]  <el_conn> connectivity
        const int offset,     ///< \param[in]  <offset>  offset for connectivity
        const int kvar0,      ///< \param[in]  <kvar0>   offset  variable for  uol
        double uold[]         ///< \param[out]  <uold>   solution
        ///<
    ) const;  ///< Return element solution
    // ----------------------------------------------------------------
    void get_el_sol_piece (
        const int ivar0,   ///< \param[in]   <ivar0>   initial variable
        const int nvars,   ///< \param[in]   <nvars>   number of variables to get  <-
        const int el_nds,  ///< \param[in]   <el_nds>  number  of element nodes for this variable
        const int iel,     ///< \param[in]  <iel> element id (ofset)
        const int offset,  ///< \param[in]  <offset>  offset for connectivity
        const int kvar0,   ///< \param[in]  <kvar0>   offset  variable for  uold
        double uold[]      ///< \param[out]  <uold>   solution
        ///<
    ) const;  ///< Return element solution

// ---------------------------------------------------------------- ---> aggiunta per p0
    void get_el_sol_p0 (
        const int Level,
        const int ivar0,   ///< \param[in]   <ivar0>   initial variable
        const int nvars,   ///< \param[in]   <nvars>   number of variables to get  <-
        const int el_nds,  ///< \param[in]   <el_nds>  number  of element nodes for this variable
        const int iel,     ///< \param[in]  <iel> element id (ofset)
        const int offset,  ///< \param[in]  <offset>  offset for connectivity
        const int kvar0,   ///< \param[in]  <kvar0>   offset  variable for  uold
        double uold[]      ///< \param[out]  <uold>   solution
        ///<
    ) const;  ///< Return element solution


    // --------------------------------------
 void get_el_nonl_sol(
    const int Level,       ///< Level (usually _NoLevels-1) <-
    const int in_step,     ///< step for x_old  0=sol 1=final 2=old <-
    const int in_var0,     ///< initial variable to get  <-
    const int in_nvar,     ///< final variable to get  <-
    const int in_dof0,     ///< initial dof  <- 
    const int in_dofn,     ///< final dof  <- 
    const int el_conn[],   ///< connectivity  compatible with Level<- element global node
    const int offset,      /// offset=mesh nodes  (_node_dof) for connectivity <-
    const int out_voffset, ///< offset  variable for  uold <-
     const int out_idoffset, ///< offset id variable for  uold <-
    double *uold           /// element node values ->
    ) const;   
    
    
    
    
    
    
    /// Return no linear  solution on element
    void get_el_nonl_sol (
        const int Level,
        const int ivar0,      ///< \param[in]   <ivar0>   initial variable
        const int nvars,      ///< \param[in]   <nvars>   number of variables to get  <-
        const int el_nds,     ///< \param[in]   <el_nds>  number  of element nodes for this variable
        const int el_conn[],  ///< \param[in]  <el_conn> connectivity
        const int offset,     ///< \param[in]  <offset>  offset for connectivity
        const int kvar0,      ///< \param[in]  <kvar0>   offset  variable for  uold
        double uold[]         ///< \param[out]  <uold>   solution
    ) const;
    // --------------------------------------
    /// Return d_aux  solution on element
    void get_el_d_aux (
        const int istep,
        const int ivar0,      // initial variable  <-
        const int nvars,      // # of variables to get  <-
        const int el_nds,     // # of element nodes for this variable  <-
        const int el_conn[],  // connectivity <-
        const int offset,     // offset for connectivity <-
        const int kvar0,      // offset  variable for  uold <-
        double uold[]         // element node values ->
    ) const;
    
    void set_cp_vector ( int nvars[],const int& vec_from, const int& vec_to ) ;

 

    // Computational Functions ========================================================================

     // Print/Read soution  =========================================================================
    // void print_ext_data() is not implemented (to do)
    // void print_xml_attrib() is defined
    // the other are defined in MGSolverDA
    // Print
    // fname=  file name    // Level=  MG level
    // of_out=  out stream file
    // n_nodes=  number of nodes n_elems=  number of elements


    // --------------------------------------------------------------------
//     virtual void print_u_xdmf_B (
//         std::ofstream& of_out,  ///< \param[in] <>  out stream file
//         int n_nodes,            ///< \param[in] <>  number of nodes
//         int n_elems,            ///< \param[in] <>  number of elements
//         std::string fname,       ///< \param[in] <>  file name
//         int n_vars,
//         std::string *var_names
//     ) const;                ///< print xml file
//     // --------------------------------------------------------------------
//     virtual void print_u (
//         std::string fname,  ///< \param[in] <>  file name
//         const int Level     ///< \param[in] <> MG level
//     ) = 0;              ///< print solution
//     // --------------------------------------------------------------------
//     virtual void read_u (
//         std::string fname,  ///< \param[in] <>  file name
//         int Level           ///< \param[in] <>  restart MG level
//     ) = 0;              ///< Read solution
// 
//     // all virtual defined in MGSOLDA_BCIC.C
//     // reading print ic
//     // ------------------------------------------------------------------------------------
//     virtual void ic_read ( /// initial conditions  from function
//         int k,             // bc from gambit    <-
//         int m,             // material from gambit    <-
//         double xp[],       // point coordinates <-
//         int iel,           // element  <-
//         double valueic[]   //  point values     ->
//     ) = 0;             ///< Reading IC  function (element)
// 
// 
//     // reading print bc ==================================================================
//     // ----------------------------------------------------------------------------------
//     virtual void bc_read (
//         int k,
//         int m,         ///< global node id
//         double x[],    ///< point vector
//         int bc_vol[],  ///< Volume flag
//         int bc_sur[]   ///< Surface flag
//     ) = 0;         ///< Reading sur BC function (element)
//     // ------------------------------------------------------------------------------------
//     /// Reading vol BC function (element)
//     virtual void bc_intern_read (
//         int k, int m,
//         double x[],    // point vector
//         int bc_vol[],  // Volume flag
//         int u[]        // value vector
//     ) = 0;
//     // -------------------------------------------------------------------------------------
//     virtual void print_bc (
//         std::string fname,  ///< \param[in] <>  file name
//         const int Level     ///< \param[in] <> MG level
//     ) = 0;              ///< print boundary conditions

// ***********************************************************************************************
//  2) NONLINEAR SOLVERS
// ***********************************************************************************************
// -----------------------------------------------------------------------------------------------
// get/set function
    double get_nonl ( int kdofs ); ///< It gets the x_nonl(kdofs) on the top mesh

// ***********************************************************************************************
// 2) CONTROL  SOLVERS
// ***********************************************************************************************
// -----------------------------------------------------------------------------------------------
// get/set function
      void set_nonl (int kdofs, double value ); 
    void set_d_aux ( int i, int kdofs, double value ); ///< It sets the d_aux(kdofs) on the top mesh
    double get_d_aux ( int i,int kdofs ); ///< It gets the d_aux(kdofs) on the top mesh
// ------------------------------------------------------------------------------------------------
// Computing functions
    virtual double MGFunctional_B (
        double parameter,  /// Use of the function: (0) compute functional OR (1) set _eta
        double& control    /// \param[in] <>  eta multiplier for optimal method
    );
// ------------------------------------------------------------------------------------------------
    virtual void  set_ctrl_dom_B (const std::vector<std::vector<double>> coord
        // const double x_min, const double x_max, const double y_min, const double y_max, const double z_min,const double z_max 
    );
    
    // virtual void  set_controlled_dom_B (
    //     const double x_min, const double x_max, const double y_min, const double y_max, const double z_min,
    //     const double z_max );
    virtual void  set_controlled_dom_B (const std::vector<std::vector<double>> coord);
    virtual void  set_lift_dom_B (const std::vector<std::vector<double>> coord);
};







#endif


