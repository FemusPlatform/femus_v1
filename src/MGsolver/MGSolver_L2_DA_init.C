// #include "Equations_tab.h"
// #include "Equations_conf.h"
// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ----------
// #include "Printinfo_conf.h"

// class files ----------------
#include "MGSolver_L2_DA.h"
// #include "MGSclass_conf.h"

// local includes --------------
#include "MGSystem_L2.h"
// #include "MGFE_conf.h"
#include "MGGeom_L1_El.h"
// #include "MGSystem.h"
#include "MGMesh_L0.h"
// #include "dense_set.h"

// #include "MGEquationsSystem.h"
// #include "MGFE.h"
// #include "MGFEMap.h"
// #include "MGGraph.h"
// #include "MGUtils.h"

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
//  #if _ndim_B==2
// using namespace DIM2;
// #endif
// #if _ndim_B==3
// using namespace DIM3;
// #endif
// ==========================================================
//               MGSolDA functions
// ==========================================================



// ===============================================================================================
/// This  function reads all the Operators from files
void MGSolDA::MGSol_init_B(
) {
  TRACKING_FUN(printf(" --> MGSolverDA.C: MGSolDA::MGSol_init: -> init_dof-> init_fAx-> init_bc \n");)
// ===============================================================================================
//  initialization of all levels: dofs and matrices;
  for(int Level = 0; Level < _NoLevels; Level++) {
    init_dof(Level);  // init the dofmap (_node_dof)
  }

  for(int Level = 0; Level < _NoLevels; Level++) {
    init_fAx(Level);      // allocation fAx
  }
  //   DEBUG(check_multiphysics(0);)
  init_bcic();  // on top level only init the bc-map
  for(int Level = 0; Level < _NoLevels; Level++) {
    init_PR(Level);      // allocation fAx
  }
  init_ntgbg();

//   init_ic();  // on top level only initial solution
// ===============================================================================================
  TRACKING_FUN(printf(" <--- MGSolverDA.C: MGSolDA::MGSol_init: -> init_dof->init_fAx-> init_bc \n  ");)
  return;
}



// ================================================================================================
/// This function is the main initialization function.
/// It sets operator dimensions r=f-Ax and Prl, Rst on the level (Level)
// ================================================================================================
void MGSolDA::init_fAx(const int Level) {
  // name directories to read
  std::string f_matrix = _mgutils_D.get_file("F_MATRIX");
  std::string f_rest   = _mgutils_D.get_file("F_REST");
  std::string f_prol   = _mgutils_D.get_file("F_PROL");
  const ParallelM::Communicator& comm1 = _mgmesh._comm.comm();
  
  // number of partial dofs ---------------------------------------------------------
  int off_proc = _NoLevels * _iproc;
  int n_local_q =(_el_dof_D[0]>0)? _mgmesh._max_nd[0][Level+off_proc]-_mgmesh._min_nd[0][off_proc]: 0;  // 0=QUADRATIC
  int n_local_l =(_el_dof_D[1]>0)? _mgmesh._max_nd[1][Level+off_proc]-_mgmesh._min_nd[1][off_proc+_offset_level_flag_D[1]*Level]:0;  // 1=LINEAR
//   int n_local_k =(_el_dof_D[2]>0)? (_mgmesh._off_el[0][Level+1+off_proc]-_mgmesh._off_el[0][Level+off_proc])*_el_dof_D[2]: 0;  // 0=Volume
  int n_local_k =(_el_dof_D[2]>0)? (_mgmesh._max_nd[2][Level+off_proc]-_mgmesh._min_nd[2][_offset_level_flag_D[0]*Level+off_proc])*_el_dof_D[2]: 0;  // 0=Volume
  //   int nel_local=_mgmesh._off_el[0][Level+off_proc+1]-_mgmesh._off_el[0][Level+off_proc];
  // global dofs  ---------------------------------------------------------------------
  int klev=Level;
  int n_glob = _Ndof_lev[klev];
  int k0_vars=(_varst_D[2]>0?_varst_D[2]:0);//int k1_vars=(_varst_D[2]<0?-_varst_D[2]:0);
  int n_local = _varst_D[1] * n_local_l + _varst_D[0] * n_local_q + k0_vars * n_local_k;
  int n0_vars[3]; n0_vars[0]=k0_vars; n0_vars[1]=_varst_D[1]; n0_vars[2]=_varst_D[0];
  //  int  nel_glob=_mgmesh._NoElements[0][Level];
  
  //   matrix ------------------------------------------------------------------------------
  std::ostringstream filename("");  filename << _mgutils_D._inout_dir << f_matrix << Level << ".h5";
  A[Level] = SparseMatrixM::build(_mgmesh._comm.comm()).release();
  A[Level]->init(n_glob, n_glob, n_local, n_local);
  
  ReadMatrix(Level,filename.str(),*A[Level],n0_vars);

  

  // vectors -------------------------------------------------------------------------------------
  b[Level]  =NumericVectorM::build(comm1).release();
  b[Level]->init(n_glob,n_local,false,AUTOMATICM);
  res[Level]=NumericVectorM::build(comm1).release(); res[Level]->init(n_glob,n_local,false,AUTOMATICM);
  x[Level]  =NumericVectorM::build(comm1).release(); x[Level]->init(n_glob,n_local,false,AUTOMATICM);
  // old solution/non linear vectors --------------------------------------------------------------
  for(int ieq = 0; ieq < 3; ieq++) {
    x_old[ieq][Level]=NumericVectorM::build(comm1); x_old[ieq][Level]->init(n_glob,false,SERIALM);
  }
  x_nonl[Level]= NumericVectorM::build(comm1);     x_nonl[Level]->init(n_glob,false,SERIALM);


//   // Restrictor -------------------------------------------------------------------------
//   if(Level < _NoLevels - 1) {
//     Rst[Level] = SparseMMatrixM::build(comm1).release(); Rst[Level]->init(n_glob,_Ndof_lev[Level+1],n_glob,_Ndof_lev[Level+1]);
//     //     int n_nodes_qp1= _mgmesh._NoNodes[Level+1];
//     filename.str(""); filename << _mgutils_D._inout_dir <<f_rest<< Level+1 <<"_"<<Level<<".h5";
//     ReadRest(Level,filename.str(),*Rst[Level],n0_vars);
//     Rst[Level]->close();
//   }
//   // Prolongation ----------------------------------------------------------------------
//   if(Level > 0) {
//     Prl[Level] = SparseMMatrixM::build(comm1).release();
//     Prl[Level]->init(n_glob,_Ndof_lev[Level-1],_Ndof_lev[Level],_Ndof_lev[Level-1]);
//     filename.str("");
//     filename << _mgutils_D._inout_dir<<f_prol<<Level-1<<"_"<<Level<<".h5";
//     ReadProl(Level, filename.str(),*Prl[Level],n0_vars);
//     Prl[Level]->close();
//     //  Prl[Level]->print_personal(filename);//print on screen Prol
//   }

  return;
}
// ================================================================================================
/// This function is the main initialization function.
/// It sets operator dimensions r=f-Ax and Prl, Rst on the level (Level)
// ================================================================================================
void MGSolDA::init_PR(const int Level) {
  // name directories to read
  std::string f_matrix = _mgutils_D.get_file("F_MATRIX");
  std::string f_rest   = _mgutils_D.get_file("F_REST");
  std::string f_prol   = _mgutils_D.get_file("F_PROL");
  const ParallelM::Communicator& comm1 = _mgmesh._comm.comm();
  // number of partial dofs ---------------------------------------------------------
//   int off_proc = _NoLevels * _iproc;
//   int n_local_q =(_el_dof_D[0]>0)? _mgmesh._max_nd[0][Level+off_proc]-_mgmesh._min_nd[0][off_proc]: 0;  // 0=QUADRATIC
//   int n_local_l =(_el_dof_D[1]>0)? _mgmesh._max_nd[1][Level+off_proc]-_mgmesh._min_nd[1][off_proc+_offset_level_flag_D[1]*Level]:0;  // 1=LINEAR
//   int n_local_k =(_el_dof_D[2]>0)? (_mgmesh._off_el[0][Level+1+off_proc]-_mgmesh._off_el[0][Level+off_proc])*_el_dof_D[2]: 0;  // 0=Volume
//   int n_local_k =(_el_dof_D[2]>0)? (_mgmesh._max_nd[2][Level+off_proc]-_mgmesh._min_nd[2][_offset_level_flag_D[0]*Level+off_proc])*_el_dof_D[2]: 0;  // 0=Volume
  //   int nel_local=_mgmesh._off_el[0][Level+off_proc+1]-_mgmesh._off_el[0][Level+off_proc];
  // global dofs  ---------------------------------------------------------------------
  int klev=Level;
  int n_glob = _Ndof_lev[klev];
  int k0_vars=(_varst_D[2]>0?_varst_D[2]:0);//int k1_vars=(_varst_D[2]<0?-_varst_D[2]:0);
//   int n_local = _varst_D[1] * n_local_l + _varst_D[0] * n_local_q + k0_vars * n_local_k;
  int n0_vars[3]; n0_vars[0]=k0_vars; n0_vars[1]=_varst_D[1]; n0_vars[2]=_varst_D[0];
  //  int  nel_glob=_mgmesh._NoElements[0][Level];

  std::ostringstream filename("");  filename << _mgutils_D._inout_dir << f_matrix << Level << ".h5";

  // Restrictor -------------------------------------------------------------------------
  if(Level < _NoLevels - 1) {
    Rst[Level] = SparseMMatrixM::build(comm1).release(); Rst[Level]->init(n_glob,_Ndof_lev[Level+1],n_glob,_Ndof_lev[Level+1]);
    //     int n_nodes_qp1= _mgmesh._NoNodes[Level+1];
    filename.str(""); filename << _mgutils_D._inout_dir <<f_rest<< Level+1 <<"_"<<Level<<".h5";
    ReadRest(Level,filename.str(),*Rst[Level],n0_vars);
    Rst[Level]->close();
  }
  // Prolongation ----------------------------------------------------------------------
  if(Level > 0) {
    Prl[Level] = SparseMMatrixM::build(comm1).release();
    Prl[Level]->init(n_glob,_Ndof_lev[Level-1],_Ndof_lev[Level],_Ndof_lev[Level-1]);
    filename.str("");
    filename << _mgutils_D._inout_dir<<f_prol<<Level-1<<"_"<<Level<<".h5";
    ReadProl(Level, filename.str(),*Prl[Level],n0_vars);
    Prl[Level]->close();
    //  Prl[Level]->print_personal(filename);//print on screen Prol
  }

  return;
}

// ================================================================================================
/// This function initializes the system degrees of freedom (dof) on the level (Level).
///  It Builds  _node_dof[Level]: map local-global (with multiple variable)
void MGSolDA::init_dof(
const int Level,    ///< Level
const int /*vb_0*/,     ///< type of element
const int /*n_vb*/  ///< number of type element
) {
// TRACKING_FUN(std::cout << " -->> MGSolDA::init_dof   Level= " << Level << std::endl;)
// ================================================================================================
  // Set up from mesh -----------------------------
  const int n_subdom = _mgmesh._n_subdom;
//   const int n_nodes  = _mgmesh._NoNodes[Level];
//   const int n_elem   = _mgmesh._NoElements[vb_0][Level];  // 0=volume

//   const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  const int* min_nd_0 = _mgmesh._min_nd[0]; //variable node/dof (0)
  const int* max_nd_0 = _mgmesh._max_nd[0]; //variable node/dof (0)
  const int* min_nd_1 = _mgmesh._min_nd[1]; //variable node/dof (1)
  const int* max_nd_1 = _mgmesh._max_nd[1]; //variable node/dof (1)
  const int* min_nd_2 = _mgmesh._min_nd[2]; //variable node/dof (2)
  const int* max_nd_2 = _mgmesh._max_nd[2]; //variable node/dof (2)
//   const int* off_el = _mgmesh._off_el[0];  // volume elements

  // number of total dofs
//   int n_nodes_1 = _mgmesh._NoNodes[_NoLevels+Level];
//   int n_nodes_2 = _mgmesh._NoNodes[3*_NoLevels+Level];
//   if(Level > 0) { n_nodes_1 = _mgmesh._NoNodes[Level - 1]; }  // pressure coarse=top+1
//   int k0_vars= (_varst_D[2]>0 ? _varst_D[2]:0);//  int k1_vars= (_varst_D[2]<0 ? -_varst_D[2]:_varst_D[0]);
  _Ndof_lev[Level] = _varst_D[2]*_mgmesh._NoNodes[_var_index_Nonodes_D[2]*_NoLevels+Level]+_varst_D[1]*_mgmesh._NoNodes[Level+_var_index_Nonodes_D[1]*_NoLevels]+_varst_D[0]*_mgmesh._NoNodes[Level];

  P_INFO(std::cout<<Level<<" total "<<_Ndof_lev[Level] << " dof q "<<_mgmesh._NoNodes[_var_index_Nonodes_D[0]*_NoLevels+Level] <<" linear dof "<< _mgmesh._NoNodes[Level+_var_index_Nonodes_D[1]*_NoLevels] <<" konst dof "<<_mgmesh._NoNodes[_var_index_Nonodes_D[2]*_NoLevels+Level] << std::endl;)



  // construction dof node vector(_node_dof) +++++++++++++++++++++
  _node_dof[Level] = new int[_n_vars_D *  _top_offset_D];
  for(int k1 = 0; k1 <  _top_offset_D * _n_vars_D; k1++) { _node_dof[Level][k1] = -1; }
  int count = 0;  int count_q=0; int count_p=0; int count_k=0;
  for(int isubdom = 0; isubdom < n_subdom; isubdom++) {
    // quadratic -----------------------------------
    int off_proc = isubdom * _NoLevels;
    for(int ivar = 0; ivar < _varst_D[0]; ivar++) {
      for(int k1 = _mgmesh._min_nd[0][off_proc+_offset_level_flag_D[0]*Level]; k1 < _mgmesh._max_nd[0][off_proc + Level]; k1++) {
        _node_dof[Level][k1 + ivar *  _top_offset_D] = count;
//         std::cout<< max_nd_0[off_proc + Level]  << " "<<  k1 << " "<<  ivar *  _top_offset_D<< " k1  _top_offset_D dof  "<<count  << " "<< min_nd_0[off_proc] << "\n ";;
        count++;
      }
    }
    count_q +=  max_nd_0[off_proc + Level]-min_nd_0[off_proc+_offset_level_flag_D[0]*Level];
    // linear -----------------------------------
    for(int ivar = 0; ivar <_varst_D[1]; ivar++) {
      for(int k1 = min_nd_1[off_proc+_offset_level_flag_D[1]*Level];  k1 < max_nd_1[Level  + off_proc]; k1++) {
        _node_dof[Level][k1 + (_varst_D[0] + ivar) *  _top_offset_D] = count;
        count++;
//         std::cout<< max_nd_0[off_proc + Level]  << " "<<  k1 << " "<<  ivar *  _top_offset_D<< " k1  _top_offset_D dof  "<<count  << " "<< min_nd_0[off_proc] << "\n ";;
      }
    }
    count_p += max_nd_1[Level  + off_proc] - min_nd_1[off_proc+_offset_level_flag_D[1]*Level];


    // constant polynomial of order _el_dof_D[2] ----------------------------------- BACKUP
    for(int ivar = 0; ivar <_varst_D[2]; ivar++) {
      for(int k1 = min_nd_2[off_proc+_offset_level_flag_D[0]*Level];  k1 < max_nd_2[Level  + off_proc]; k1++) {
        _node_dof[Level][k1 + (_varst_D[0] + _varst_D[1]+ivar) *  _top_offset_D] = count;
        count++;
      }
    }

    //     // constant polynomial of order _el_dof_D[2] -----------------------------------
    // for(int ivar = 0; ivar <_varst_D[2]; ivar++) {
    //   for(int k1 = min_nd_2[off_proc+_offset_level_flag_D[0]*Level];  k1 < max_nd_2[Level  + off_proc]; k1++) {
    //     int elem_gidx = (k1 +  _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]) * _el_dof_D[0];
    //     int k2 = _mgmesh._el_map[0][elem_gidx + 8];
    //     _node_dof[Level][k2 + (_varst_D[0] + _varst_D[1]+ivar) *  _top_offset_D] = count;
    //     count++;
    //   }
    // }

    count_k += max_nd_2[Level  + off_proc] - min_nd_2[off_proc+_offset_level_flag_D[1]*Level];  
  }
  std::cout<< " dof quad lin face k "<< count_q << " " << count_p << " " << count_q-count_p-count_k << " " << count_k << std::endl;
// =========================================================================================================
// TRACKING_FUN(std::cout << " <-- MGSolDA::init_dof   Level= " << Level << std::endl;)
  return;
}

// ============================================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::init_bcic() {  // ========================================================================

  /// A) Set up: mesh,dof, bc
  // mesh ----------------------------------------------------------------------
//   const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int ntot_el = 0; for(int ilev=0; ilev<_NoLevels; ilev++) {ntot_el +=_mgmesh._NoElements[0][ilev];}

  // Dof ----------------------------------------------------------------------
//   const int n_kb_dofs=((_varst_D[2]>0)? _face_dof_D[2]: 0);  // surface dofs
//   const int n_pb_dofs=((_varst_D[1]>0)? _face_dof_D[1]: 0); // get_n_shapes(_ndim_B-2);
//   const int n_ub_dofs=((_varst_D[0]>0)? _face_dof_D[0] :0); // get_n_shapes(_ndim_B-2);
  //   const int  n_dofs =  n_pb_dofs*_varst_D[1] + n_ub_dofs*_varst_D[2];
  const int n_k_dofs = ((_varst_D[2] > 0) ? _el_dof_D[2] : 0);                    // volume dofs
  const int n_l_dofs = ((_varst_D[1] > 0) ? _el_dof_D[1] : 0);  // get_n_shapes(_ndim_B-1);
  const int n_u_dofs = ((_varst_D[0] > 0) ? _el_dof_D[0] : 0);  // get_n_shapes(_ndim_B-1);

  // set 1 all the points for  bc (boundary condition) ------------------------
  for(int i1 = 0; i1 <  _top_offset_D; i1++) {_bc[0][i1] = 11; _bc[1][i1] = 11;  }

  // **************************************************************************
  // B) Reading  face_id vector (boundary zones) if the dataset exists
  // Open an existing file ----------------------------------------------------
  std::ostringstream file_bc; file_bc << _mgutils_D._inout_dir << _mgutils_D.get_file("INMESH");  //"/mesh.h5";
//   P_INFO(std::cout << " Reading bc_id from= " << file_bc.str() << std::endl;)
  hid_t file_id = H5Fopen(file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dims[2];
  // face id vector initial setting (set 0) -----------------------------------
  int* face_id_vect = new int[ _top_offset_D]; for(int i=0; i< _top_offset_D; i++) { face_id_vect[i] = 0; }
  // Getting dataset ----------------------------------------------------------
  hid_t dtset = H5Dopen(file_id, "NODES/BC", H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dtset); /* Get filespace handle first. */
  hid_t status = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status < 0) {   std::cerr << "GenIc::read dims not found"; abort();}
  else {  // reading (otherwise it stays 0)
    assert((int)dims[0] ==  _top_offset_D);
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &face_id_vect[0]);
  }
  H5Dclose(dtset);  H5Sclose(filespace);

  // **************************************************************************
  // C) Reading  mat_id (volume zones) if the dataset exists
  // mat id vector initialization ---------------------------------------------
  int* mat_id_vect = new int[ntot_el]; for(int i = 0; i < ntot_el; i++) {mat_id_vect[i] = 1; }
  // Getting dataset --------------------------------------------------------
  dtset = H5Dopen(file_id, ("ELEMS" + std::to_string(_mgeqnmap_D._mgmesh._GeomEl.n_q[0]) +"/MAT/MAT").c_str(), H5P_DEFAULT);
  filespace = H5Dget_space(dtset); /* Get filespace handle first. */
  status = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status < 0) {std::cerr << "GenIc::read mat dims not found"; abort();}
  else {  // reading if the dataset exists (otherwise it stays 1) ------------
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat_id_vect[0]);
  }  // end reading ---------------------------------------------------------

  // clean --------------------------------------------------------------------
  H5Sclose(filespace);  H5Dclose(dtset);  H5Fclose(file_id);

  // ----------------------------------------------------------------------------------------------
  // boundary condition from user
  GenBc_loop(0, _el_dof_D[0], n_u_dofs, n_l_dofs, n_k_dofs, face_id_vect, mat_id_vect);
  // ----------------------------------------------------------------------------------------------
  GenIc_loop(0, _el_dof_D[0],n_u_dofs, n_l_dofs, n_k_dofs, face_id_vect, mat_id_vect);
  // clean
  delete[] face_id_vect; delete[] mat_id_vect;

  return;
}
// ================================================================================================
void MGSolDA::GenIc_loop(
const int /*vvvb*/,          // 0-> volume   1-> boundary
const int /*ndof_femv*/,     // number of elment nodes
const int /*n_u_dofs*/,      // element quad dofs
const int /*n_l_dofs*/,      // element linear dofs
const int /*n_k_dofs*/,  // element konst dofs
int */*face_id_vect[]*/,      // group bc from gambit
int mat_id_vect[]        // bc from gambit
) {

  NumericVectorM& sol_top = *x[_NoLevels - 1];  // solution (top level)
  NumericVectorM& old_sol_top = *x_old[0][_NoLevels-1];
  NumericVectorM& nl_sol_top = *x_nonl[_NoLevels - 1];
  NumericVectorM& oold_sol_top = *x_old[1][_NoLevels - 1];
  const int* node_dof_top = _node_dof[_NoLevels - 1];
  const double* xyz_glb = _mgmesh._xyz;
  std::vector<double> u_value(_n_vars_D);    double xp[3];
  int off_proc = _iproc * _NoLevels; // loop reading
  int ndof_lev = 0;
  for(int pr = 0; pr < _mgmesh._iproc; pr++) {
    ndof_lev +=_mgmesh._off_el[0][pr*_NoLevels+_NoLevels-1+1] -_mgmesh._off_el[0][pr*_NoLevels + _NoLevels - 1];
  }

  for(int iel = 0; iel < _mgmesh._off_el[0][off_proc + _NoLevels] - _mgmesh._off_el[0][off_proc + _NoLevels - 1]; iel++) {
    int elem_gidx = (iel +  _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]) * _el_dof_D[0];
    int mat_id_elem = mat_id_vect[iel + ndof_lev];

    // the local nodes
    for(int i = 0; i < _el_dof_D[0]; i++) {
      int k = _mgmesh._el_map[0][elem_gidx + i];
      for(int idim = 0; idim < _ndim_B; idim++) xp[idim] = xyz_glb[k + idim * _top_offset_D];
       int face_id_node = _bc[0][k];
      // ===================================================
      ic_read_D(face_id_node, mat_id_elem, xp, iel + _mgmesh._off_el[0][off_proc + _NoLevels-1], u_value.data());
      // ===================================================

      // Set the quadratic and linear fields ----------------------------------------------------
      for(int ivar = 0; ivar < _varst_D[0]; ivar++) {
        int irrr = node_dof_top[k + ivar *  _top_offset_D];
        assert(irrr != -1);
        sol_top.set(irrr, u_value[ivar]);
      }  // -----------------------------------------------------------------------------------------

      // ----------Set linear vertex field --------------------------------------------
      if(_var_index_D[1] == 1 && i < _el_dof_D[1]) {
        for(int ivar = 0; ivar < _varst_D[1]; ivar++) {
          assert(node_dof_top[k + ivar *  _top_offset_D] != -1);
          sol_top.set(node_dof_top[k + (ivar+_varst_D[0])  *  _top_offset_D], u_value[ivar+_varst_D[0]]);
        } // -----------------------------------------------------------------------------------------
      }
      // --------------- Set linear mid-point field --------------------------------------
      if(_var_index_D[1] == 3 && (i >=_el_dof_D[0]-1-_mgmesh._GeomEl._n_sides[0] && i<_el_dof_D[0]-1)) {
        for(int ivar = 0; ivar < _varst_D[1]; ivar++) {
          assert(node_dof_top[k + ivar *  _top_offset_D] != -1);
          sol_top.set(node_dof_top[k + (ivar+_varst_D[0])  *  _top_offset_D], u_value[ivar+_varst_D[0]]);
        } // -----------------------------------------------------------------------------------------
      }
      //  --------------- Set discontinuous fields ---------------------------------
      if(i == _el_dof_D[0] - 1) {
        for(int ivar = 0; ivar < _varst_D[2]; ivar++) {
          assert(node_dof_top[k + (ivar + _varst_D[0] + _varst_D[1]) * _top_offset_D] != -1);
          sol_top.set(node_dof_top[k+(ivar+_varst_D[0]+_varst_D[1])*_top_offset_D] ,u_value[_varst_D[0]+_varst_D[1]+ivar]);
        }  // -----------------------------------------------------------------------------------------
      }
    } // i
  }
  // end of element loop

  sol_top.localize(nl_sol_top);
  sol_top.localize(oold_sol_top);
  old_sol_top.zero();

  // x[_NoLevels - 1]->print();

  return;
}


