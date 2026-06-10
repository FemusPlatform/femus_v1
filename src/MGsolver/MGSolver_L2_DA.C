// #include "Equations_tab.h"
// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ----------
// #include "Printinfo_conf.h"
// #include "MGSclass_conf.h"
// #include "Equations_conf.h"



// class files ----------------
#include "MGSolver_L2_DA.h"


// // local includes --------------
#include "MGSystem_L2.h"
// // #include "MGFE_conf.h"
#include "MGGeom_L1_El.h"
// // #include "MGSystem.h"
#include "MGMesh_L0.h"
// // #include "dense_set.h"
// 
// // #include "MGEquationsSystem.h"
// // #include "MGFE.h"
// #include "MGFEMap_L0.h"
// #include "MGGraph_L0.h"
// #include "MGUtils_L0.h"

// #ifdef HAVE_MED
// #include "InterfaceFunctionM.h"
// #include "MEDCouplingFieldDouble.hxx"
// #include "MEDCouplingUMesh.hxx"
// #include "MEDLoader.hxx"
// #endif

// algebric includes -----------
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "numeric_vectorM.h"
#include "sparse_MmatrixM.h"
#include "sparse_matrixM.h"

// ==========================================================
//               MGSolDA functions
// ==========================================================

/// This function is the MGSolDA constructor I
MGSolDA::MGSolDA(
  MGEquationsSystem& mg_equations_map_in,  //
  const int nvars_in[],                    // # of quad variables
  std::string eqname_in,       // equation name
  std::string /*varname_in*/)  // basic variable name
  :   MGSolBase(mg_equations_map_in._mgmesh),                     ///< father constructor
    _mgeqnmap_D(mg_equations_map_in),          // equation map pointer
    _mgutils_D(mg_equations_map_in._mgutils),  // mgutils pointer from equation map pointer
    _mgfemap_D(mg_equations_map_in._mgfemap),  // mgfemap pointer from equation map pointer
    _eqname_D(eqname_in),                               ///< equation name
    _n_vars_D(( nvars_in[0]>0? nvars_in[0]:0  ) + nvars_in[1] + nvars_in[2]){   ///< n ofvars (quad,lin,const)

//  int ns_type= stoi ( _mgutils_D._dict_config["MG_NavierStokes"]);

  // Reading parameters from  GeometrySetting.in
  int err=_mgutils_D._dict_geom.size();
  // label[0,1,2] defined by  VART0,1,2
  int type_nonodes[3];
  type_nonodes[0]=(int)(_mgutils_D._dict_geom["VART0"]);
  type_nonodes[1]=(int)(_mgutils_D._dict_geom["VART1"]);
  type_nonodes[2]=(int)(_mgutils_D._dict_geom["VART2"]);
  std::string label[3]; std::string dict_label[4]={"Q","V","F","C",};
//  for(int ktype=0;ktype<3;ktype++) label[ktype]= dict_label[type_nonodes[ktype]];
 // restart
//  const int restart = stoi(_mgutils._dict_config["restart"]);
if(_mgutils_D._dict_geom.size()-err) {
  std::cout<< "MGMesh::MGMesh: Missing parameter in GeometrySetting.in"; std::exit(EXIT_FAILURE);}

 _dt_B=stod(_mgutils_D._dict_config["dt"]);                        // parameter  dt
  //  Lagrange piecewise  variable types (type 2(const),1(linear),0(quad))
  _varst_D[2] = nvars_in[0];  // Lagrange piecewise constant variables
  _varst_D[1] = nvars_in[1];  // Lagrange piecewise linear variables
  _varst_D[0] = nvars_in[2];  // Lagrange piecewise linear variables
  _dir_D = 0;   // number equation for segregated system
  _NumRestartSol = 1;
   _top_offset_D = _mgmesh._NoNodes[_NoLevels-1];   // mesh nodes
  _el_sides_D= _mgmesh._GeomEl._n_sides[0];    // element sides
  // piecewise constant linear, quadratic and (...cubic): volume[0-2] and face element[3-5]
  for(int klq=0; klq<3; klq ++) { //quad-> 0 lin ->1 const->2
    _fe_D[klq] = _mgfemap_D.get_FE(klq); // Lagrange constant(0),linear(1),quad(2)
    _el_dof_D[2-klq] =  _fe_D[klq]->get_NoShape() ;  // Lagrange piecewise constant variables
     if(_ndim_B> 1) {
       _fe_face_D[klq]=_mgfemap_D.get_FE(klq+3);
       _face_dof_D[2-klq] = _fe_face_D[klq]->get_NoShape();}
     else {_fe_face_D[2-klq]=NULL   ;_face_dof_D[2-klq] =0;}

//     if(_ndim_B> 1) { _el_doft_D[2-klq+3] =   _fe_D[klq]->get_NoShape();} // face
//     else {_el_doft_D[2-klq+3] =0;}
  }

  // Raviart Thomas
//   for(int klq=1; klq<2; klq++) { //linear -> 1
//     _fe_D[4] = _mgfemap_D.get_FE(6); // Lagrange constant(0),linear(1),quad(2)
//     _fe_face_D[1] = _mgfemap_D.get_FE(4);
//   s
//     if(_ndim_B> 1) {
//          _el_dof_D[4] =  _fe_D[4]->get_NoShape() ;  // Lagrange piecewise constant variable
//        _fe_face_D[4]=_mgfemap_D.get_FE(7);
//       _face_dof_D[4] =   _fe_face_D[4]->get_NoShape();} // face
//     else {_face_dof_D[4] =0;}
//    }
  _offset_level_flag_D[0]=1; _offset_level_flag_D[1]=0; _offset_level_flag_D[2]=0;

if (type_nonodes[1]==2){// Raviart Thomas
     _fe_D[4] = _mgfemap_D.get_FE(6); // Lagrange constant(0),linear(1),quad(2)
//        _fe_face_D[4]=_mgfemap_D.get_FE(7);
   _offset_level_flag_D[1] = 1;
     _var_index_D[1]=3;
    _var_index_Nonodes_D[1]=2;
    // -----------------   RT
    _el_dof_D[1] = (nvars_in[1] > 0) ? _el_sides_D  : 0;
     if(_ndim_B > 1) _face_dof_D[1]=1;
}

  // System names and units
  _var_names_D.resize(_n_vars_D);  // names
  _var_ref_D.resize(_n_vars_D);

  // output string stream
  int icount=0;
  for(int iname2 = 0; iname2 < _varst_D[0]; iname2++) {   // type 0 quadratic
    std::ostringstream ostr; ostr << "q" << iname2; _var_names_D[icount] = ostr.str();
    _var_ref_D[icount] = 1; icount++;
  }
  for(int iname1 = 0; iname1 <  _varst_D[1]; iname1++) {   // type 1 linear
    std::ostringstream ostr;  ostr << "l" << iname1; _var_names_D[icount] = ostr.str();
    _var_ref_D[icount] = 1; icount++;
  }
  for(int iname0 = 0; iname0 < abs(_varst_D[2]); iname0++) {   // type 0 piecewise
    std::ostringstream ostr; ostr<<"k"<<iname0; _var_names_D[icount]=ostr.str();
    _var_ref_D[icount] = 1; icount++;
  }
   _mtr_nrows_D=0; for(int deg=0;deg<3;deg++) _mtr_nrows_D +=_varst_D[deg]*_el_dof_D[deg];

 // fields -> initalization
 for(int k = 0; k < 60; k++) {    // coupling  basic system fields
//   const int idx = _DAdata_eq_D[2].tab_eqs[k];
  _FF_idx_D[k] =  -1;
 }
  _AxiSym_D=0;// _AxiSym_D =0 no  _AxiSym_D =1 axis (x-axis); _AxiSym_D =2   fixed plane (x or y or z);
#ifdef TWO_PHASE_LIB
      _msolcc(NULL),
#endif
  return;
}
// =========================================
/// This function is the destructor
// =========================================
MGSolDA::~MGSolDA() {}
// ===============================================================================================
/// This function controls the time step operations:
// ===============================================================================================
int MGSolDA::MGTimeStep(const double time, const int  /*i_step*/, const int /*n_iter*/, const int mode) {
  TRACKING_FUN(printf(" --> MGSolBase.C: MGSolBase::MGTimeStep \n  ");)
  int error = 0;

  // [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
  if(mode%2 == 0){
    std::cout<< "INIT GEN_MAT_RHS ASSEMBLY FOR EQUATION: " <<_eqname_D<< std::endl;
    TIME(std::clock_t start_time = std::clock();)
    GenMatRhs(time, _NoLevels - 1, 1);

    // [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
    for(int Level = 0; Level < _NoLevels - 1; Level++) GenMatRhs(time, Level, 0);

    TIME(std::clock_t end_time=std::clock();std::cout<<"Assembly time ="<< double(end_time-start_time)/CLOCKS_PER_SEC<<" s "<<std::endl;)
    // [c] Solution of the linear system (MGSolverBase::MGSolve).
    error = MGSolve_B(1.e-8, 40);

    TIME(std::clock_t end_timef=std::clock();std::cout<<"Assembly + solution time ="<<double(end_timef-start_time)/CLOCKS_PER_SEC<<"s "<< std::endl;)
  }

  // [d] Update of the old solution at the top Level
  if(mode>0){
    std::cout << std::endl << "UPDATE SOLUTION " <<_eqname_D<< std::endl;
    std::cout<<"-------------------------------------\n";
    x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);
    x[_NoLevels - 1]->localize(*x_old[1][_NoLevels - 1]);
    // std::cout<<std::endl;

    // x_old[1][_NoLevels-1]->localize(*x_old[2][_NoLevels-1]);
    // x[_NoLevels-1]->localize(*x_old[0][_NoLevels-1]);
    // int nnodes = _mgmesh._NoNodes[_NoLevels-1];
    // for(int i = 0; i < nnodes; i++) {
    //   int k = _node_dof[_NoLevels - 1][i];
    //   if(k >- 1) {
    //       double value0 = get_sol(0,k);
    //       double value1 = get_sol(1,k);
    //       double value  = value0 + value1;
    //       set_sol(1,k,value);
    //   }
    // }
  }

  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::MGTimeStep \n  ");)
  return error;
}
// =====================================================================
double MGSolDA::CalcFUpwind(double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes) {

  double vel_modulus = 1.e-10; for(int i = 0; i < Dim; i++) { vel_modulus += VelOnGauss[i] * VelOnGauss[i]; }
  vel_modulus = sqrt(vel_modulus);
  double h_eff = 1.e-20; for(int i = 0; i < NbOfNodes; i++) {
    double hh = 1.e-20;   for(int idim = 0; idim < Dim; idim++)  hh += VelOnGauss[idim] * PhiDer[i + idim * NbOfNodes];
    h_eff += fabs(hh);
  }
  h_eff = 2. * vel_modulus/ h_eff;
  if(h_eff < 1.e-10) { h_eff = 1.;  std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n"; }
  // STANDARD SUPG
  const double Pe_h = 0.5 * vel_modulus * h_eff / (Diffusivity+ 1.e-10);
  const double a_opt = (1. / tanh(Pe_h) - 1. / Pe_h);
  if(a_opt > 1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
  //   f_upwind = 0.5*a_opt*h_eff/ (vel_modulus);
  double f_upwind = 0.5 * a_opt * h_eff / (vel_modulus);
  return f_upwind;
}
// =====================================================================
double MGSolDA::MGFunctional_B(double, double&) {
  std::cout << "Not implemented in MsolverDA";
  return 0.;
}
// =====================================================================
void MGSolDA::compute_jac(
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
) const {
  // All data vectors are _el_dof_D[0] long
  const double alfa = 1.e-08;
  // variable loop
  for(int ivar = 0; ivar < nvars; ivar++) {
    u_forw[ivar] = 0.;
    u_back[ivar] = 0.;  // set zero
    for(int jdim = 0; jdim < _ndim_B; jdim++) {
      u_forw_dx[ivar * _ndim_B + jdim] = 0.;  // set zero
      u_back_dx[ivar * _ndim_B + jdim] = 0.;  // set zero
    }
    // interpolation with shape functions
    for(int eln = 0; eln < n_shape; eln++) {
      const int indx = eln + ivar * _el_dof_D[0];
      if(indx == j + idim * _el_dof_D[0]) {
        u_forw[ivar] += phi[eln] * (uold_b[indx] + alfa);
        u_back[ivar] += phi[eln] * (uold_b[indx] - alfa);
      }
      else {
        u_forw[ivar] += phi[eln] * uold_b[indx];
        u_back[ivar] += phi[eln] * uold_b[indx];
      }

      for(int jdim = 0; jdim < _ndim_B; jdim++) {
        if(indx == j + idim * _el_dof_D[0]) {
          u_forw_dx[ivar * _ndim_B + jdim] += dphi[eln + jdim * n_shape] * (uold_b[indx] + alfa);
          u_back_dx[ivar * _ndim_B + jdim] += dphi[eln + jdim * n_shape] * (uold_b[indx] - alfa);
        }
        else {
          u_forw_dx[ivar * _ndim_B + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
          u_back_dx[ivar * _ndim_B + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
        }
      }
    }
  }

  return;
}
// =========================================
void MGSolDA::set_xooold2x() {
  for(int Level = 1; Level <= _NoLevels; Level++) {
    /// A. Setup
    const int offset = _mgmesh._NoNodes[Level - 1];  // fine level # of nodes
    //     const int pie_offset = (NSUBDOM) * (_mgmesh._NoElements[0][Level - 1]);  // fine level # of nodes

    // reading loop over system varables
    const int iel_min=_mgmesh._off_el[0][_iproc*_NoLevels+Level-1];
    const int n_elem_pl=_mgmesh._off_el[0][_iproc*_NoLevels+Level]-iel_min;
    for(int ivar = 0; ivar < _varst_D[0] + _varst_D[1]; ivar++) {
      int el_nds = (ivar >= _varst_D[0])? _el_dof_D[1]:_el_dof_D[0];

      // storing  ivar variables (in parallell)
      for(int iel=0; iel<n_elem_pl; iel++) {
        int elem_gidx = (iel+iel_min)*_el_dof_D[0];
        for(int i = 0; i < el_nds; i++) {             // linear and quad
          int kdof = _node_dof[_NoLevels - 1][_mgmesh._el_map[0][elem_gidx + i] + ivar * offset];  // the global node
          const double value = (*x_old[2][Level - 1])(kdof);
          x[Level - 1]->set(kdof, value);  // set the field
        }
      }
    }

    int ndof_lev = 0;
    for(int pr = 0; pr < _mgmesh._iproc; pr++) {
      int delta =_mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
      ndof_lev += delta;
    }

    /// D. delocalization and clean
    x[Level - 1]->localize(*x_old[0][Level - 1]);    x[Level - 1]->localize(*x_old[1][Level - 1]);
  }

  return;
}
//  =========================================================================================
// EQUATION ACTIVATION (START THE CLASS POINTER DURING EXECUTION
//  =========================================================================================
// int Order,      = order: const,linear,quad
// int Field,      = Field to activate
// std::string sfn,= system field name
// int& n_index,  = n_index (collecting index)
// int vector,    = vector
// int neqs       = dimension (number of eqs) neqs
// int coupled    = coupled (1) or segregated (0)  solver
// ==========================================================================================
void MGSolDA::ActivateVectField(
  int Order,        ///< order: const,linear,quad
  int Field,        ///< Field to activate
  std::string sfn,  ///< SystemFieldName
  int& n_index,     ///< n_index (collecting index)
  int coupled       ///< coupled (1) or segregated (0)  solver
) {
  std::string FieldX = sfn+"X";  std::string FieldY =sfn+"Y";  std::string FieldZ =sfn+"Z";

  if(coupled == 0) {   // flag 0 in SimulationConfiguration -> UNcoupled
    if(sfn.compare("TAU") == 0) {//   if (sfn.compare("TAU") == 0)
      std::string FieldX = sfn + "XX"; std::string FieldY = sfn + "XY"; std::string FieldZ = sfn + "YY";
      ActivateEquation(Order, Field, FieldX, n_index);
      ActivateEquation(Order, Field + 1, FieldY, n_index);
      ActivateEquation(Order, Field + 2, FieldZ, n_index);
    }
    else {    //  if (sfn.compare("TAU") != 0)
      ActivateEquation(Order, Field, FieldX, n_index);
      ActivateEquation(Order, Field + 1, FieldY, n_index);
if(_ndim_B == 3)  ActivateEquation(Order, Field + 2, FieldZ, n_index);
    }
  }
  else {

    // flag 1 in SimulationConfiguration -> coupled
    _DAdata_eq_D[Order].tab_eqs[Field] = n_index;                  // table
    _DAdata_eq_D[Order].mg_eqs[n_index] = _mgeqnmap_D.get_eqs(sfn);  // FSI equation pointer
    _DAdata_eq_D[Order].indx_ub[n_index + 1] =_DAdata_eq_D[Order].indx_ub[n_index] + _ndim_B;  // _DAdata_eq_D[2].ub index
    _DAdata_eq_D[Order].n_eqs++;                           // number of quadratic system
    n_index++;                                         // update counter

//     if(Order == 1){
//     // flag 1 in SimulationConfiguration -> coupled
//     _DAdata_eq_D[Order].tab_eqs[Field] = n_index;                  // table
//     _DAdata_eq_D[Order].mg_eqs[n_index] = _mgeqnmap_D.get_eqs(sfn);  // FSI equation pointer
//     _DAdata_eq_D[Order].indx_ub[n_index + 1] =_DAdata_eq_D[Order].indx_ub[n_index] + 1;  // _DAdata_eq_D[2].ub index
//     _DAdata_eq_D[Order].n_eqs++;                           // number of quadratic system
//     n_index++;                                         // update counter
//     }
  }

  return;
}
// ==========================================================================================
void MGSolDA::ActivateDA(
   int Order,        ///< order: const,linear,quad
   int Field,        ///< Field to activate
   std::string sfn,  ///< SystemFieldName
   int& n_index     ///< n_index (collecting index)
//   int coupled       ///< coupled (1) or segregated (0)  solver
) {
      if(_varst_D[2-Order]>0){
     _DAdata_eq_D[Order].tab_eqs[Field] = n_index;                  // table
     _DAdata_eq_D[Order].mg_eqs[n_index] = _mgeqnmap_D.get_eqs(sfn);  //  equation pointer
     _DAdata_eq_D[Order].indx_ub[n_index + 1] = _DAdata_eq_D[Order].indx_ub[n_index] + _varst_D[2-Order];  // _DAdata_eq_D[2].ub index
     _DAdata_eq_D[Order].n_eqs++;      // number of quadratic/lin/konst equation
  }
  n_index++;
  return;
}
// ================================================================================================
void MGSolDA::ActivateControl(
  int Order,        ///< order: const,linear,quad
  int Field,        ///< Field to activate             // Field to activate
  std::string sfn,  ///< SystemFieldName // system field name
  int& n_index,     ///< n_index (collecting index)             // n_index (collecting index)
  int vector,       ///< coupled (1) or segregated (0)  solver
  int neqs          ///< dimension (number of eqs)
) {
  // flag 2 in SimulationConfiguration.in ->  split
  // flag 1 in SimulationConfiguration -> coupled
  std::string FieldX = sfn + "X";  std::string FieldY = sfn + "Y";  std::string FieldZ = sfn + "Z";
  ActivateEquation(Order, Field, FieldX, n_index);
  if(vector == 1) {
    ActivateEquation(Order, Field + 1, FieldY, n_index);
    if(neqs == 3) { ActivateEquation(Order, Field + 2, FieldZ, n_index); }
  }

  return;
}
// ==========================================================================================
// This function activates a scalar equation
void MGSolDA::ActivateScalar(
  int Order,        ///< order: const,linear,quad
  int Field,        ///< Field to activate
  std::string sfn,  ///< System Field Name
  int& n_index      ///< n_index (collecting index)
) {
  ActivateEquation(Order, Field, sfn, n_index);
  return;
}
// ==========================================================================================
// This function activates coupled equations
void MGSolDA::ActivateCoupled(
  int Order,        ///< order: const,linear,quad
  int Field,        ///< Field 1 to activate
  std::string sfn,  ///< System Field Name
  int& n_index,     ///< n_index (collecting index)
  std::string sfn2  ///< Field 2 to activate
) {
  ActivateEquation(Order, Field, sfn, n_index);
  ActivateEquation(Order, Field + 1, sfn2, n_index);
  return;
}
// ==========================================================================================
// This function activates a single  equation
void MGSolDA::ActivateEquation(
  int Order,        ///< order: const,linear,quad
  int Field,        ///< Field to activate
  std::string sfn,  ///< System Field Name
  int& n_index      ///< n_index (collecting index)
) {
  _DAdata_eq_D[Order].tab_eqs[Field] = n_index;                  // table
  _DAdata_eq_D[Order].mg_eqs[n_index] = _mgeqnmap_D.get_eqs(sfn);  //  equation pointer
  _DAdata_eq_D[Order].indx_ub[n_index + 1] = _DAdata_eq_D[Order].indx_ub[n_index] + _varst_D[2-Order];  // _DAdata_eq_D[2].ub index
  _DAdata_eq_D[Order].n_eqs++;                                                      // number of quadratic system
  n_index++;  // update counter
  return;
}
// ===============================================================================================
/// This function sets up data structures for each problem class
void MGSolDA::setUpExtFieldData() {

  TRACKING_FUN(printf(" --> MGSolBase.C: MGSolBase::setUpExtFieldData \n  ");)

  /// A) set up _mg_eqs
  // external system and index vectors
  for (int deg = 0; deg < 3; deg++) {
    for (int kl = 0; kl < _DAdata_eq_D[deg].max_neqs; kl++) {
      _DAdata_eq_D[deg].n_eqs = 0;
      _DAdata_eq_D[deg].mg_eqs[kl] = NULL;
      _DAdata_eq_D[deg].indx_ub[kl] = -1;
      _DAdata_eq_D[deg].tab_eqs[kl] = -1;

      for (int kk = 0; kk < _el_dof_D[2-deg]; ++kk) {
        _DAdata_eq_D[deg].ub[kk + kl *_el_dof_D[2-deg]] = 0.;  // data
      }
    }
  }
  // start index K from 0, L from 0, Q from _ndim_B (coordinates+ q variable)
  _DAdata_eq_D[0].indx_ub[0] = 0;  //_DAdata_eq_D[0].n_eqs=0; // piecewice constant  (0)
  _DAdata_eq_D[1].indx_ub[0] = 0;  //_DAdata_eq_D[1].n_eqs=0; // piecewice linear    (1)
  _DAdata_eq_D[2].indx_ub[0] = 0;  //_DAdata_eq_D[2].n_eqs=0; // piecewice quadratic (2)

  TRACKING_FUN(printf(" <-- MGSolBase.C: MGSolBase::setUpExtFieldData \n  ");)
  return;
}
//============================================================================
#ifdef TWO_PHASE_LIB
void MGSolBase::set_mgcc(MGSolCC& cc) { _msolcc = &cc; }
#endif

// Two-phase solver
#ifdef TWO_PHASE
#include "MGSolverCC.h"
#endif

#ifdef TWO_PHASE_LIB
//============================================================================
void MGSolBase::set_mgcc(MGSolCC& cc) { _msolcc = &cc; }
#endif
