// alg include-------------------
#include "numeric_vectorM.h"

// local inlude -----------------
#include "MGMesh_L0.h"
#include "MGSystem_L1.h"
// #include "MGSolver_L1_Base.h"
#include "MGSolver_L1_Base.h"
#include "MGSolver_L2_DA.h"
#include "MGGeom_L1_El.h"
// ========================================================================
// CONSTRUCTOR-DESTRUCTOR
// ========================================================================
MGSolBase::MGSolBase(
  MGMesh& mgmesh_in                    ///<  MG mesh
//   std::string eq_name_in       ///< equation name (default Base0)
) :
  MGSolBase0(mgmesh_in/*,eq_name_in*/) {
//   TRACKING_FUN(printf(" --> MGSolverBase.C: MGSolBase::MGSolBase (Constructor ******) \n  ");)
// ================================================================================================
  // allocation of dynamic system ---------------------------------------------
  x_old[0].resize(_NoLevels);  // old solution 1 step (old)
  x_old[1].resize(_NoLevels);  // old solution 2 step (oold)
  x_old[2].resize(_NoLevels);  // old solution 3 step (ooold)
  x_nonl.resize(_NoLevels);    // non linear solution
  x_aux.resize(0);
  d_aux.resize(0);

  // Set up boundary conditions ++++++++++++++++++++++++++++
//   if(Level == _NoLevels - 1) {
    int top_nodes=  _mgmesh._NoNodes[_NoLevels-1];
    _bc[0].resize(top_nodes); _bc[1].resize(top_nodes);
    for(int k1 = 0; k1 < top_nodes; k1++) { _bc[0][k1] = 11; _bc[1][k1] = 1;}
//   }
  // class parameters ----------------------------------------------------------
  _control = 0.;
  _dt_B = 0.1;  // time step
  _ndim_B=_mgmesh._dim;
  _NumRestartSol=1; // euler discretization (1 x_old)
// ================================================================================================
//   TRACKING_FUN(printf(" <-- MGSolverBas1e.C: MGSolBase::MGSolBase (Constructor ******) \n  ");)
}

//  ==========================================================
// This build function allows the use of template 
// class MGGeomEl<unsigned const int DIM, GEOFamily geofam > calling 
// thos father MGGeomElBase(const unsigned int dim,     const GEOFamily geofam)
std::unique_ptr<MGSolBase> MGSolBase::build (
  MGMesh& mgmesh_in,
    FIELDS field_type_in,
     MGEquationsSystem& mg_equations_map_in,  //
            const int nvars_in[],                    // # of quad variables
            std::string eqname_in,       // equation name
            std::string varname_in
    
) {//  ==========================================================
 
        switch (field_type_in){
          case MG_DA:
            return std::unique_ptr<MGSolBase> (new MGSolDA( 
            mg_equations_map_in,  //
            nvars_in,                    // # of quad variables
            eqname_in,       // equation name
            varname_in)
            );
          // case MG_Temperature:
          //    return std::unique_ptr<MGSolBase> (new MGSolT( 
          //    MGEquationsSystem& mg_equations_map_in,  //
          //    const int nvars_in[],                    // # of quad variables
          //    std::string eqname_in,       // equation name
          //    std::string /*varname_in*/)
          //    );            //  case LAGRANGE14:
              // return std::unique_ptr<MGGeomElBase>(new MGGeomElLag14_1D(1, LAGRANGE14));
                default:  std::cout <<" MGSolSolver not known "; abort();
          }
 
    return nullptr;
}


// ===============================================================================================
/// Destructor (level structure)
MGSolBase::~MGSolBase() {
  TRACKING_FUN(printf(" --> MGSolverBase.C: MGSolBase::MGSolBase (Destructor ******) \n  ");)
// ===============================================================================================
  clear(); // clear substructrures
  x_nonl.clear();  // nonlinear solution tmp
  x_old[0].clear();  x_old[1].clear();  x_old[2].clear();  // old solutions
  x_aux.clear();     d_aux.clear();

// ===============================================================================================
  TRACKING_FUN(printf(" <-- MGSolverBase.C: MGSolBase::MGSolBase (Destructor ******) \n  ");)
}

// ===============================================================================================
/// Substructure destructor
void  MGSolBase::clear() {
  TRACKING_FUN(printf(" --> MGSolverBase.C: MGSolBase::clear \n  ");)
// ===============================================================================================
// ===============================================================================================
  TRACKING_FUN(printf(" <-- MGSolverBase.C: MGSolBase::clear \n  ");)

}

// ***********************************************************************
//  TIME DEPENDENT SOLVER
// ***********************************************************************
/// ======================================================
/// This function change current dt
/// ======================================================
void MGSolBase::set_dt(double dt) {_dt_B=dt;} ///< MG time step solver (backward Euler)
double MGSolBase::GetValue(int /*flag*/) {return 2;}    //just to fix warning, func not implemented
void MGSolBase::SetValueVector(std::vector<double> /*value*/) {}

// ========================================================================
//  SET FUNCTIONS
// ========================================================================
void MGSolBase::set_nonl(int k, double val) {
  x_nonl[_NoLevels - 1]->set(k, val);
  return;
}  // set the field
// ========================================================================
void MGSolBase::set_x_aux(int i, int k, double val) {
  x_aux[i]->set(k, val);
  return;
}  // set the field
// ========================================================================
void MGSolBase::set_d_aux(int i, int k, double val) {
  d_aux[i]->set(k, val);
  return;
}
// ========================================================================
void MGSolBase::set_sol(int i, int k, double val) {
  x_old[i][_NoLevels - 1]->set(k, val); //< i: 0=sol[0]  1=sol[1] 2=sol[2]
  return;
}
//  GET FUNCTIONS
// ========================================================================
double MGSolBase::get_x_aux(int i, int k) { return ((*x_aux[i])(k)); }
double MGSolBase::get_d_aux(int i, int k) { return ((*d_aux[i])(k)); }
double MGSolBase::get_sol(int i, int k) { return ((*x_old[i][_NoLevels - 1])(k)); }

// ***********************************************************************
//  NONLINEAR SOLVER
// ***********************************************************************
double MGSolBase::get_nonl(int kdofs) { return (*x_nonl[_NoLevels-1])(kdofs); }

// ***********************************************************************
// control  SOLVER
// *********************************************************************

/// ======================================================
/// This function  computes the  functional
/// ======================================================
double MGSolBase::MGFunctional_B(
  double /*parameter*/,  /// Use of the function: (0) compute functional OR (1) set _eta
  double& /*control*/    /// \param[in] <>  eta multiplier for optimal method
) {
  std::cout << "Wrong use of MGFunctional from MGSolBase.C, aborting";
  abort();
  return 1;
}

/// =======================================================================
/// This function sets the controlled domain for optimal control problems
/// =======================================================================
void MGSolBase::set_ctrl_dom_B(const std::vector<std::vector<double>> coord
  // const double xMin, const double xMax,
  // const double yMin, const double yMax,
  // const double /*zMin*/, const double /*zMax*/
) {
  TRACKING_FUN(printf(" --> MGSolBase.C: MGSolBase::set_ctrl_dom \n  ");)
// ===============================================================================================
  int el_conn[27];
  double x_m[3];
  double xx_qnds[27 * 3];
  const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * (_mgmesh._n_subdom - 1)];
  double eps = 1.e-6;                                                            // tolerance for coordinates
  const int nel_e = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc];      // stop element
  int n_target_region = coord[0].size();
  _weight_ctrl.resize(n_target_region);
  for(int iel = 0; iel < (nel_e - nel_b); iel++) {
    _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds);  // gets element coordinates
    for(int idim = 0; idim < _ndim_B; idim++) {
      x_m[idim] = 0;
      for(int d = 0; d < n_geom_nodes; d++) {
        x_m[idim] += xx_qnds[idim * n_geom_nodes + d] / n_geom_nodes;
      }                             // end d loop
    }                               // end idim loop

    for(uint ii = 0; ii < n_target_region; ii++) {
      _weight_ctrl[ii].resize(n_elem);
      _weight_ctrl[ii][iel + nel_b] = 0;  // default is 0
      double xMin = coord[0][ii];
      double xMax = coord[1][ii];
      double yMin = coord[2][ii];
      double yMax = coord[3][ii];
      double zMin = coord[4][ii];
      double zMax = coord[5][ii];

      if(x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps)
if(_ndim_B == 3) if( x_m[2] > zMin - eps && x_m[2] < zMax + eps ) {
        _weight_ctrl[ii][iel + nel_b] = 1;  // 1 inside control region
      }
      else _weight_ctrl[ii][iel + nel_b] = 1;
    }
  }  // end iel loop
// ===============================================================================================
  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::set_ctrl_dom \n  ");)
  return;
}
/// =======================================================================
/// This function sets the controllable domain for optimal control problems
/// =======================================================================
void MGSolBase::set_controlled_dom_B(const std::vector<std::vector<double>> coord
  // const double xMin, const double xMax,
  // const double yMin, const double yMax,
  // const double /*zMin*/, const double /*zMax*/
) {
  TRACKING_FUN(printf(" --> MGSolBase.C: MGSolBase::set_controlled_dom \n  ");)
// ===============================================================================================
  int el_conn[27];
  double x_m[3];
  double xx_qnds[27 *3];
      const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * (_mgmesh._n_subdom - 1)];
  double eps = 1.e-6;                                                            // tolerance for coordinates
  const int nel_e = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc];      // stop element
  int n_ctrl_region = coord[0].size();
  
  _weight_controlled.resize(n_elem);
  for(int iel = 0; iel < (nel_e - nel_b); iel++) {
    _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds);  // gets element coordinates
    for(int idim = 0; idim < _ndim_B; idim++) {
      x_m[idim] = 0;
      for(int d = 0; d < n_geom_nodes; d++) {
        x_m[idim] += xx_qnds[idim *n_geom_nodes + d] / n_geom_nodes;
      }                             // end d loop
    }                               // end idim loop
    _weight_controlled[iel + nel_b] = 0;  // default is 0
 if(_ndim_B == 3)   for(uint ii = 0; ii < n_ctrl_region; ii++) {
      double xMin = coord[0][ii];
      double xMax = coord[1][ii];
      double yMin = coord[2][ii];
      double yMax = coord[3][ii];
      double zMin = coord[4][ii];
      double zMax = coord[5][ii];
      
      if(x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps   && x_m[2] > zMin - eps && x_m[2] < zMax + eps
        ) {
        _weight_controlled[iel + nel_b] = 1;  // 1 inside control region
      }
    }
 if(_ndim_B == 2)   for(uint ii = 0; ii < n_ctrl_region; ii++) {
      double xMin = coord[0][ii];
      double xMax = coord[1][ii];
      double yMin = coord[2][ii];
      double yMax = coord[3][ii];
      
      if(x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps) {
        _weight_controlled[iel + nel_b] = 1;  // 1 inside control region
      }
    }    
    
    
    
    
  }  // end iel loop
// ===============================================================================================
  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::set_controlled_dom \n  ");)
  return;
}

/// =======================================================================
/// This function sets the lift domain for optimal control problems
/// =======================================================================
void MGSolBase::set_lift_dom_B(const std::vector<std::vector<double>> coord
  // const double xMin, const double xMax,
  // const double yMin, const double yMax,
  // const double /*zMin*/, const double /*zMax*/
) {
  TRACKING_FUN(printf(" --> MGSolBase.C: MGSolBase::set_controlled_dom \n  ");)
// ===============================================================================================
  int el_conn[27];
  double x_m[3];
  double xx_qnds[27 * 3];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * (_mgmesh._n_subdom - 1)];
   const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  double eps = 1.e-6;                                                            // tolerance for coordinates
  const int nel_e = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc];      // stop element
  int n_ctrl_region = coord[0].size();

  _weight_lift.resize(n_elem);
  for(int iel = 0; iel < (nel_e - nel_b); iel++) {
    _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds);  // gets element coordinates
    for(int idim = 0; idim < _ndim_B; idim++) {
      x_m[idim] = 0;
      for(int d = 0; d < n_geom_nodes ; d++) {
        x_m[idim] += xx_qnds[idim * n_geom_nodes  + d] / n_geom_nodes ;
      }                             // end d loop
    }                               // end idim loop
    _weight_lift[iel + nel_b] = 0;  // default is 0
 if(_ndim_B == 3)   for(uint ii = 0; ii < n_ctrl_region; ii++) {
      double xMin = coord[0][ii];
      double xMax = coord[1][ii];
      double yMin = coord[2][ii];
      double yMax = coord[3][ii];
      double zMin = coord[4][ii];
      double zMax = coord[5][ii];

      if(x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps  && x_m[2] > zMin - eps && x_m[2] < zMax + eps
        ) {
        _weight_lift[iel + nel_b] = 1;  // 1 inside control region
      }
    }
    
 if(_ndim_B == 2)   
    for(uint ii = 0; ii < n_ctrl_region; ii++) {
      double xMin = coord[0][ii];
      double xMax = coord[1][ii];
      double yMin = coord[2][ii];
      double yMax = coord[3][ii];
     

      if(x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps ) {
        _weight_lift[iel + nel_b] = 1;  // 1 inside control region
      }
    }
    
    
  }  // end iel loop
// ===============================================================================================
  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::set_controlled_dom \n  ");)
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_sol(
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    ) const {
       const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nds; id++) {
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top = _node_dof[_NoLevels - 1][el_conn[id] + (ivar + ivar0) * offset];  // dof from top level
      uold[id + (ivar + kvar0) * n_geom_nodes ] = ((*x_old[i_step][_NoLevels-1])(kdof_top)); // element sol
    }  // end quadratic
  }
  return;
}
 // ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_sol_F(
    const int level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nfaces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->

    ) const {
       const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nfaces; id++) {
      int id1=id+n_geom_nodes-1-el_nfaces;
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top = _node_dof[level][el_conn[id1] + (ivar + ivar0) * offset];  // dof from top level
          assert(kdof_top>-1);
      uold[id1 + (ivar + kvar0) * n_geom_nodes] = ((*x_old[i_step][level])(kdof_top));  // element sol
    }  // end quadratic
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_sol_C(
    const int level,
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int /*el_nctr*/,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    ) const {
 const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
      int id = 0;
      int id1=id+ n_geom_nodes-1;
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top = _node_dof[level][el_conn[id1] + (ivar + ivar0) * offset];  // dof from top level
          assert(kdof_top>-1);
      uold[id1 + (ivar + kvar0) *  n_geom_nodes] = ((*x_old[i_step][_NoLevels-1])(kdof_top));  // element sol

    }  // end quadratic
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
/// For       Q (quad) V (vertex)   F(face)                     C(center)
/// in_dof0   0        0          _el_dof_D[0]-1-NDOF_FACE;        _el_dof_D[0]-1
/// in_dofn  _el_dof_D[0]  _el_dof_D[0]P   NDOF_FACE                       1
///
void MGSolBase::get_el_nonl_sol(
    const int Level,       ///< Level (usually _NoLevels-1) <-
    const int /*in_step*/,     ///< step for x_old  0=sol 1=final 2=old <-
    const int in_var0,     ///< initial variable to get  <-
    const int in_nvar,     ///< number of variables to get  <-
    const int in_dof0,     ///< initial dof  <-
    const int in_dofn,     ///< number of dofs  <-
    const int el_conn[],   ///< connectivity  compatible with Level<- element global node
    const int offset,      ///< offset=_mgmesh._NoNodes[Level] mandatory !!!!!!!
    const int out_voffset, ///< offset  variable for  uold <-
    const int out_idoffset, ///< offset id variable for  uold <-
    double *uold           /// element node values ->
    ) const {
const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id1 =0; id1<in_dofn; id1++) {
    const int id=id1+in_dof0;
    // quadratic
    for (int ivar1 = 0; ivar1 < in_nvar; ivar1++) {  // ivar is like idim
      const int indx_in= el_conn[id]+(in_var0+ivar1)*offset;// offset[Level]
      const int kdof_top=_node_dof[Level][indx_in];  // dof from  Level
      const int indx_out= id+out_idoffset+(ivar1+out_voffset)*n_geom_nodes; // local index []uold
      assert(kdof_top>-1); //  el_conn[] incompatible with Level
      uold[indx_out] =((*x_nonl[Level])(kdof_top));  // element sol
    }  // end quadratic
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_nonl_sol_F(
     const int Level,
    const int /*i_step*/,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nfaces,     // # of faces  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    ) const {
const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nfaces; id++) {
      int id1=id+n_geom_nodes-1-el_nfaces;
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top=_node_dof[Level][el_conn[id1]+(ivar+ivar0)*offset];  // dof from top level
      assert(kdof_top>-1);
      uold[id1+(ivar+kvar0)*n_geom_nodes] =((*x_nonl[_NoLevels-1])(kdof_top));  // element sol
    }  // end quadratic
  }
  return;
}
 // ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_nonl_sol_C(
    const int Level,     ///< Level (usually _NoLevels-1) <-
    const int /*i_step*/,    ///< step for x_old  0=sol 1=final 2=old <-
    const int ivar0,      ///< initial variable  <-
    const int nvars,      ///< # of variables to get  <-
    const int /*el_nctr*/,    ///< # of faces  <-
    const int el_conn[],  // connectivity <- element global node
    const int offset,     // offset=mesh nodes  (_node_dof) for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double *uold         // element node values ->
    ) const {
const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
    int id1=n_geom_nodes-1; // center points
    //  center points
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top = _node_dof[Level][el_conn[id1]+(ivar+ivar0)*offset];  // dof from top level
      assert(kdof_top>-1);
      uold[id1 + (ivar + kvar0) * n_geom_nodes] = ((*x_nonl[Level])(kdof_top));  // element sol
    }  // end quadratic
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_sol_piece(
    const int ivar0,   // initial variable  <-
    const int nvars,   // # of variables to get  <-
    const int el_nds,  // # of element nodes for this variable  <-
    const int iel,     // connectivity <-
    const int offset,  // offset for connectivity <-
    const int kvar0,   // offset  variable for  uold <-
    double uold[]      // element node values ->
    ) const {
      const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nds; id++) {
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top =
          _node_dof[_NoLevels - 1][iel * el_nds + id + (ivar + ivar0) * offset];      // dof from top level
      uold[id + (ivar + kvar0) *n_geom_nodes] = ((*x_old[0][_NoLevels - 1])(kdof_top));  // element sol
    }  // end quadratic
  }
  return;
}

void MGSolBase::get_el_sol_p0(
    const int Level,
    const int ivar0,   // initial variable  <-
    const int nvars,   // # of variables to get  <-
    const int el_nds,  // # of element nodes for this variable  <-
    const int iel,     // connectivity <-
    const int offset,  // offset for connectivity <-
    const int kvar0,   // offset  variable for  uold <-
    double uold[]      // element node values ->
    ) const {
      const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nds; id++) {
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      // const int kdof_top =
      //     _node_dof[_NoLevels - 1][iel * el_nds + id + (ivar + ivar0) * offset];      // dof from top level
      uold[id + (ivar + kvar0) *n_geom_nodes] = ((*x_old[1][Level])(iel));  // element sol
    }  // end quadratic
  }
  return;
}

// ==========================================================================================
/// This function gets  the solution  vector at the nodes of  an element.
void MGSolBase::get_el_nonl_sol(
  const int Level,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double uold[]         // element node values ->
    ) const {
      const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nds; id++) {
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivarq is like idim
      const int kdof_top =
          _node_dof[Level][el_conn[id] + (ivar + ivar0) * offset];  // dof from top level
      double val = ((*x_nonl[_NoLevels - 1])(kdof_top));
      uold[id + (kvar0 + ivar) * n_geom_nodes] = val;  // element sol
    }  // end quadratic
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*_el_dof_D[0] with _el_dof_D[0] max dof (quad)
void MGSolBase::get_el_d_aux(
    const int istep,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double uold[]         // element node values ->
    ) const {
      const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
  for (int id = 0; id < el_nds; id++) {
    // quadratic
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivarq is like idim
      const int kdof_top =
          _node_dof[_NoLevels - 1][el_conn[id] + (ivar + ivar0) * offset];               // dof from top level
      uold[id + (kvar0 + ivar) * n_geom_nodes] = (*d_aux[istep])(kdof_top);  // element sol
    }  // end quadratic
  }
  return;
}
// ==========================================================================
/// This function copies values from "vec_from" to "vec_to" vectors.
/// Available vectors are ordered as: x->0, x_old->1, x_oold->2, x_ooold->3,
/// x_nonl->4, disp->5, disp_old->6,  disp_oold->7
// ==========================================================================
void MGSolBase::set_cp_vector(int nvars[],const int& vec_from, const int& vec_to) {
  for (int Level = 0; Level <= _NoLevels - 1; Level++) {
    const int offset = _mgmesh._NoNodes[Level];  // fine level # of nodes
const int n_geom_nodes = _mgmesh._GeomEl.n_q[0];
const int n_geom_lin_nodes = _mgmesh._GeomEl.n_q[1];
    for (int ivar = 0; ivar < nvars[2] + nvars[1]; ivar++) {
      int el_nds = n_geom_nodes;
      if (ivar >= nvars[2])   el_nds = n_geom_lin_nodes;  // quad and linear
      for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++)
        for (int iel = 0; iel < _mgmesh._off_el[0][iproc * _NoLevels + Level + 1] -
                                    _mgmesh._off_el[0][iproc * _NoLevels + Level];
             iel++) {
          int elem_gidx = (iel + _mgmesh._off_el[0][iproc * _NoLevels + Level]) * n_geom_nodes;

          for (int i = 0; i < el_nds; i++) {            // linear and quad
            int k = _mgmesh._el_map[0][elem_gidx + i];  // the global node
            double value = 0.;

            switch (vec_from) {
              case 1: value = (*x_old[0][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 2: value = (*x_old[1][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 3: value = (*x_old[2][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 4: value = (*x_nonl[Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 5: value = (*d_aux[0])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 6: value = (*d_aux[1])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 7: value = (*d_aux[2])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 8:
                value = (*d_aux[3])(_node_dof[_NoLevels - 1][k + ivar * offset]);
                break;
                //                case 9: value = (*x_oooold[Level])(_node_dof[_NoLevels - 1][k + ivar *
                //                offset]); break;
              default: std::cout << "Incorrect vec_from number in set_uoold function" << endl; break;
            }

            switch (vec_to) {
              case 1: x_old[0][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 2: x_old[1][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 3: x_old[2][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 4: x_nonl[Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 5:
                d_aux[0]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 6:
                d_aux[1]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 7:
                d_aux[2]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 8:
                d_aux[3]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
                //
                //               case 9:
                //                 x_oooold[Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);
                //                 // set the field break;

              default: cout << "Incorrect vec_to number in set_uoold function" << endl; break;
            }
          }
        }
    }

    if (vec_from == 1 && vec_to == 3) x[Level]->localize(*x_old[1][Level]);  // for backward compatibility
  }

  int ndof_lev = 0;
  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    int delta =
        _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
    ndof_lev += delta;
  }
  return;
}
