// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ---------->
// #include "Printinfo_conf.h"

// algebric includes -----------
// #include "dense_matrixM.h"
// #include "dense_vectorM.h"
// #include "numeric_vectorM.h"
// #include "sparse_MmatrixM.h"
// #include "sparse_matrixM.h"

//  Femus class files ---------------->
#include "MGSolver_L2_DA.h"
// #include "MGSclass_conf.h"
// #include "MGFE_conf.h"
// #include "EquationSystemsExtendedM.h"
#include "MGGeom_L1_El.h"
// #include "MGSystem.h"
#include "MGMesh_L1_Extended.h"
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



// ==========================================================
//               MGSolDA functions
// ==========================================================

// ****************************************************************************
// Boundary condition and INITIAL CONDITION
// ****************************************************************************



// ========================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::GenBc_loop(
  const int vvvb,          // 0-> volume   1-> boundary
  const int ndof_femv,     // number of elment nodes
  const int /*n_u_dofs*/,      // element quad dofs
  const int /*n_l_dofs*/,      // element linear dofs
  const int /*n_k_dofs*/,  // element konst dofs
  int face_id_vect[],      // group bc from gambit
  int mat_id_vect[]        // bc from gambit

) {  // ======================================================
  TRACKING_FUN(printf(" --> MGSolDA_BCIC.C: MGSolDA::GenBc_loop \n  ");)
//   const int n_var = _varst_D[2] + _varst_D[1] + _varst_D[0];
  /// A) Set up: mesh,dof, bc
  // mesh
  const double* xyzgl = _mgmesh._xyz;
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  //   double normal[3];
  double xp[3];
  double xxb_qnds[27 * 3];
  double xx_qnds[27 * 3];
//   double normal[3];
  int el_neigh[27];
//   int el_flag[27];
  int el_conn[27]/*, elb_conn[9]*/;
  int sur_toply[9];  // boundary topology
  //   int  *bc_Neu  =new int[n_dofs]; // element bc
  //   int  *bc_value=new int[n_dofs]; // element bc
  int bc_Neu[3];    // new int[n_dofs]; // element bc
  int bc_value[3];  // new int[n_dofs]; // element bc
  int face_id_node = 0;
  int mat_id_elem = 0;
  int ndof_lev = 0;
  const int el_sides = _mgmesh._GeomEl._n_sides[0];

  // initial setting -----------------------------------------------------------------------------
  for(int isub = 0; isub < _mgmesh._n_subdom; ++isub) {
    int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
    int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
    int delta = ielf - iel0;
    for(int iel = 0; iel < delta; iel++) {   // element loop
      _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds, isub);
      _mgmesh.get_el_neighbor(el_sides, 0, _NoLevels - 1, iel, el_neigh, isub);
      for(int i = 0; i < ndof_femv; i++) {                               // node lement loop
        const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + i];  // global node
        // coordinates
        //         for(int idim = 0; idim < _ndim_B; idim++) { xp[idim] = xyzgl[k + idim * offset]; }
        if(_node_dof[_NoLevels - 1][k] > -1) { _bc[0][k] = -1000; }
      }
    }
  }  // --------------------------------------------------------------------------------------------

  /// B) Element Loop to set  bc[] (which is a node vector) ---------------------------------
  for(int isub = 0; isub < _mgmesh._n_subdom; ++isub) {
    int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
    int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
    int delta = ielf - iel0;

    for(int iel = 0; iel < delta; iel++) {   // element loop
      _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds, isub);
      _mgmesh.get_el_neighbor(el_sides, 0, _NoLevels - 1, iel, el_neigh, isub);
      mat_id_elem = 0;  // v==1 (boundari) no mat_id
      if(vvvb == 0) { mat_id_elem = mat_id_vect[iel + iel0]; }
      
      // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      // volume -----------------------------------------------------------------------------------
      for(int i = 0; i < ndof_femv; i++) {                               // node lement loop
        const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + i];  // global node

        if(_node_dof[_NoLevels - 1][k] > -1) {
          for(int idim = 0; idim < _ndim_B; idim++) {
            xp[idim] = xyzgl[k + idim * offset];
          }  // coordinates
          bc_value[0] = 1;
          bc_Neu[0] = 11;  // variable loop
          face_id_node = face_id_vect[k];
          // boundary (face_id_vect) and volume zones (mat_id_elem)
          bc_intern_read_D(face_id_node, mat_id_elem, xp, bc_Neu, bc_value);  // (idgroup, mat,x,.,.)
//           int old_val = _bc[0][_node_dof[_NoLevels - 1][k]];
           int old_val = _bc[0][k];
          if(old_val == -1000) { _bc[0][k] = bc_Neu[0]; }
          //           _bc[0][_node_dof[_NoLevels - 1][k]] = (old_val == -1000) ? bc_Neu[0] : (old_val);

          // sharing boundary nodes on the same element
//           for(int ivar = 0; ivar < n_dofs; ivar++) {
            // int kdof = _node_dof[_NoLevels - 1][k];  // kdof <-k
            // if(kdof > -1) {
            //   int number = _bc[1][k] / 10000;  // number of old count for double pt in an element
            //   // multiple corner
            //   if(abs(number) == 1) {
            //     _bc[1][k] = _bc[1][k] - _bc[1][k] * 10000 /
            //                    abs(_bc[1][k]);  // if the count is 1 set 0 if >1 leave
            //   }
            // }
//           }
        }
      }  // i-loop
        // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
     

      // -----------------------------  -----------------------------------------------------------------
      // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      // Boundary
      // _bc[0][all points] always assigned 
      for(int iside = 0; iside < el_sides; iside++) {
        if(el_neigh[iside] <0) {
          // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          int bc_side=abs(el_neigh[iside]);
        

//           int bc_min=1000; int face_min=1000;// max bc_flag  max bc_flag
//           // setup boundary element -> connectivity+coordinates
//           int bc_face=-1;
           for(int lbnode =_face_dof_D[0]-1; lbnode > -1; lbnode--) {
             int lnode = _mgmesh._GeomEl._surf_top(lbnode + _face_dof_D[0] * iside);  // local nodes
             sur_toply[lbnode] = lnode;                                          // lbnode -> lnode
             const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[lbnode]];
            //  if(k==11){
            //    int aa=1;
            // }
//            
//            
//             // read from user function -> bc_read_D(   ) ----------------------------------------------------
               for(int d=0;d<_ndim_B;d++) xp[d]/*=xxb_qnds[d*_face_dof_D[0]+lbnode]*/=xx_qnds[d*_el_dof_D[0]+lnode];              
//               bc_value[0] = 1;  bc_Neu[0] = 11;
//               face_id_node = face_id_vect[k];  // boundary (face_id_vect)
//  
                bc_read_D(bc_side, mat_id_elem, xp, bc_Neu, bc_value);
               int bc_face=bc_Neu[0];
//               bc_read_D(face_id_node, mat_id_elem, xp, bc_Neu, bc_value);
//          
//                if(lbnode ==_face_dof_D[0]-1)   bc_face=bc_Neu[0]; // midpoint
              int  tbc=_bc[0][k]; //old_value
              if(tbc<  bc_face)  _bc[0][k]=  bc_face; // rewrite the high bc-number
          }  // ----------------------------------------------------------------------------------------------
      
      
      
//       for(int iside = 0; iside < el_sides; iside++) {
//         std::cout<< el_neigh[iside] << " "; 
//         if(el_neigh[iside] <0) {
//  int bc_side=abs(el_neigh[iside]);
// //           int bc_min=1000; int face_min=1000;// max bc_flag  max bc_flag
//           // setup boundary element -> connectivity+coordinates
//           int bc_face=-1;
//           for(int lbnode =_face_dof_D[0]-1; lbnode > -1; lbnode--) {
//             int lnode = _mgmesh._GeomEl._surf_top(lbnode + _face_dof_D[0] * iside);  // local nodes
//             sur_toply[lbnode] = lnode;                                          // lbnode -> lnode
//             const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[lbnode]];
//            
// 
//             // read from user function -> bc_read_D(   ) ----------------------------------------------------
//               for(int d=0;d<_ndim_B;d++) xp[d]=xxb_qnds[d*_face_dof_D[0]+lbnode]=xx_qnds[d*_el_dof_D[0]+lnode];              
//               bc_value[0] = 1;  bc_Neu[0] = 11;
//               face_id_node = face_id_vect[k];  // boundary (face_id_vect)
//  
//               bc_read_D(face_id_node, mat_id_elem, xp, bc_Neu, bc_value);
// 
//               if(lbnode ==_face_dof_D[0]-1)   bc_face=bc_Neu[0]; // midpoint
//               int  tbc=_bc[0][k]; //old_value
//               if(tbc<  bc_face)  _bc[0][k]=  bc_face; // rewrite the high bc-number
//           }  // ----------------------------------------------------------------------------------------------
          
          
//           int face_nodes=_face_dof_D[0];
//           if(_varst_D[0]==0) { face_nodes=NDOF_PB; }
//           
// //           // we reassign bc -> face with bc higher gets the points
//           for(int lbnode = 0; lbnode < face_nodes; lbnode++) {
//             const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[lbnode]];
//             const int kdof=_node_dof[_NoLevels - 1][k];
//             if(kdof > -1) {
//               const int ibc_now=  _bc[0][k];
//               if(ibc_now <  bc_min) { _bc[0][k]=bc_min; }
//             }
// //                printf(" GenBc_loop: %d face %d  %d kdof %d k %d bc %d  bc_min %d \n",iel , iside,face_id_node,kdof,k,_bc[0][k],bc_min);
//           }  // ----------------------------------------------------------------------------------------------
// //           printf(" \n");
// 
// 
//           // normal
//           _fe_D[2]->normal_g(xxb_qnds, normal);
//           int dir_maxnormal = (fabs(normal[0]) > fabs(normal[1])) ? 0 : 1;
//           dir_maxnormal= (fabs(normal[dir_maxnormal])>fabs(normal[_ndim_B-1]))?dir_maxnormal:_ndim_B-1;
//           // computation of k_el and face_id_el
//           int k_el = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[_face_dof_D[0] - 1]];  // global node
//           int k_el_dof = _node_dof[_NoLevels - 1][k_el];
//           int bc_el = (k_el_dof > -1) ? (int)_bc[0][k_el] % 100 : 0;
// 
//           // ++++++++++++++++++++++++++++++++++++++++++++++++++++
//           const int k_face =
//             _mgmesh._el_map[0][(iel+iel0)*ndof_femv + sur_toply[_face_dof_D[0] - 1]];  // global node
//           int k_face_dof = _node_dof[_NoLevels - 1][k_face];
//           int face_id_mid_face = (k_face_dof > -1) ? face_id_vect[k_face] : 0;
//           int bc_face = (k_face_dof > -1) ? (int)_bc[0][k_face] : 0;
//           bc_face = bc_face % 100;
// 
//           for(int i = 0; i < _face_dof_D[0]; i++) {                                          // node lement loop
//             const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[i]];  // global node
//             int kdof = _node_dof[_NoLevels - 1][k];
// 
//             if(kdof > -1) {
//               face_id_node = face_id_vect[k];
//               int bc_id = (int)_bc[0][k];  // label surface pt 00
//               int mynormal = dir_maxnormal;    // dir_normal=geometrical normal
//               int sign = (bc_id == 0) ? 1 : (bc_id / (abs(bc_id)));
//               if((bc_id % 1000) / 100 > 0) { mynormal = (bc_id % 1000) / 100 - 1; }
//               if((bc_id % 1000) / 100 >3) {   //   pts label x00 ----------------------------------------------------------
// //                 std::cout << " ... single pressure pts .. \n";
// //                 if (i < NDOF_PB) { _bc[1][_node_dof[_NoLevels - 1][k + _ndim_B * offset]] = 4; }
//                 mynormal = abs((bc_id % 1000) / 100) % 4;  // force the normal
//               }
//               bc_id = bc_id % 100;  // pts label x00 forced normal
// 
//               //   pts label 00 ----------------------------------------------------------
//               // set the local boundary conditions into global vector bc[]
//               for(int ivar = 0; ivar < _varst_D[0]; ivar++) {   // quad el --
//                 int kdof = _node_dof[_NoLevels - 1][k + (ivar)*offset];
//                 int number = abs(_bc[1][k] / 10000) + 1;
//                 _bc[1][k] += sign * 10000;  // updating number of common nodes
//                 if(abs(_bc[1][k]) < 10000 || bc_id == bc_face) {
//                   _bc[1][k] = sign * (abs(bc_id) + (mynormal + 1) * 1000 + number * 10000);
// 
//                 }
// 
// 
//               }  // ----------------------------------------------------------------
//             }
//           }  // loop i +++++++++++++++++++++++++++++++++++++++++++++
         }    // iside -1
         
       }      // -----------------------------  End Boundary -------------------------------------
          // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
     }        // end of element loop
// 
     ndof_lev += delta;
   }  // i-sub

//   for(int inode=0;inode<25;inode++){
//    int kdof= _node_dof[_NoLevels - 1][inode];
//     if(kdof>-1) std::cout<< inode << " dof "<< kdof << " bc "  << _bc[1][inode] <<" \n";
//     else  std::cout<< inode << " dof "<< kdof << " bc "  << "  " <<" \n";
//   }
//   DEBUG(check_bc(0);)


  TRACKING_FUN(printf(" <--- MGSolDA_BCIC.C: MGSolDA::GenBc_loop \n  ");)
  return;
}

// ========================================
/// This function  defines the boundary conditions for the system:
void MGSolDA::bc_intern_read_D(
  int /*face_id_node*/,  ///<  face identity           (in)
  int /*mat_flag*/,      ///<  volume identity         (in)
  double /*xp*/[],       ///< xp[] node coordinates    (in)
  int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
  int bc_flag[]          ///< boundary condition flag  (out)
) {                        // ===================================
  /// Default: all Neumann
  //   for(int ivar=0; ivar<_n_vars; ivar++) {
  bc_flag[0] = 1;
  bc_Neum[0] = 11;

  //   }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void MGSolDA::get_el(
  const int Level,                   // level
  const int ivar0,                   // initial variable  <-
  const int nvars,                   // # of variables to get  <-
  const int el_nds,                  // # of element nodes for this variable  <-
  const int el_conn[],               // connectivity <-
  const int offset,                  // offset for connectivity <-
  std::vector<int>& el_dof_indices,  // element DOFs->
  int bc_dofs[][27],           // element boundary cond flags ->
  double uold[]                      // element node values ->
) const {                          // ==============================================================
  for(int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for(int ivar = ivar0; ivar < ivar0 + nvars; ivar++) {        // ivarq is like idim
      const int indx_loc = id + ivar * _el_dof_D[0];                 // local (element) index
      const int indx_glob = el_conn[id] + ivar * offset;         // global (mesh) index
      const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

      el_dof_indices[indx_loc] = _node_dof[Level][indx_glob];  // from mesh to dof
//       bc_dofs[0][indx_loc] = _bc[0][kdof_top];                 // element bc
//       bc_dofs[1][indx_loc] = _bc[1][kdof_top];                 // element bc
       bc_dofs[0][indx_loc] = _bc[0][ el_conn[id]];                 // element bc
      bc_dofs[1][indx_loc] = _bc[1][ el_conn[id]];                 // element bc
      uold[indx_loc] = (*x_old[0][_NoLevels - 1])(kdof_top);   // element sol
    }  // end quadratic ------------------------------------------------
  }

  return;
}

// ==============================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void MGSolDA::get_el_dof_bc(
  const int Level,  // level
  const int /*iel*/,    // element subdomain
  //   const int nvars[],       // # of variables to get  <-
//   const int el_nds[],                // # of element nodes for this variable  <-
  const int el_conn[],               // connectivity <-
  const int offset,                  // offset for connectivity <-
  std::vector<int>& el_dof_indices  // element connectivity ->
//   int bc_vol[],                      // element boundary cond flags ->
   // int bc_bd[]                        // element boundary cond flags ->
) const {                          // ==============================================================

  int offset_f=0;
  if(_var_index_D[1]==3)  offset_f=_el_dof_D[0]-1-_el_sides_D;
  for(int id = 0; id < _el_dof_D[0]; id++) {
    // quadratic -------------------------------------------------
    if(id < _el_dof_D[0])
      for(int ivar = 0; ivar < _varst_D[0]; ivar++) {   // ivarq is like idim
//         const int indx_loc = id + ivar * _el_dof_D[0];
        const int indx_loc_ql = id + ivar *  _el_dof_D[0];
        const int indx_glob = el_conn[id] + ivar * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
            assert( kdof_top !=-1);
        //         std::cout << el_dof_indices[indx_loc_ql] << " " << _node_dof[Level][indx_glob] << "\n";
        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
         assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
            // bc_bd[indx_loc_ql] = _bc[0][el_conn[id]];                         // element bc
//         bc_vol[indx_loc] = _bc[0][el_conn[id]];                        // element bc
      }  // end quadratic ------------------------------------------------

    //     // linear -----------------------------
    if(id <  _el_dof_D[1])
      for(int ivar = 0; ivar < _varst_D[1]; ivar++) {   // ivarq is like idim
        //        const int  indx_loc_l = id +ivar* _el_dof_D[1];
//         const int indx_loc = id + (ivar + _varst_D[0]) * _el_dof_D[0];
        const int indx_loc_ql = id + ivar *  _el_dof_D[1] + _varst_D[0] *  _el_dof_D[0];
        const int indx_glob = el_conn[id +offset_f] + (ivar + _varst_D[0]) * offset;
//         const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
//      assert( kdof_top !=-1);
        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
         assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                       // element bc
         // bc_bd[indx_loc_ql] = _bc[0][el_conn[id ] ];                         // element bc
//         bc_vol[indx_loc] = _bc[0][el_conn[id +offset_f] ];                       // element bc
        
        
      }  // end quadratic ------------------------------------------------

    //     // piecewise -----------------------------
    if(id == _el_dof_D[0]-1)
      for(int ivar = 0; ivar < _varst_D[2]; ivar++) {   // ivarq is like idim
        //        const int  indx_loc_l = id +ivar* el_nds[1];
//         const int indx_loc = 0 + (ivar + _varst_D[0] + _varst_D[1]) * _el_dof_D[0];
        const int indx_loc_ql = 0 + ivar * _el_dof_D[2] + _varst_D[0] * _el_dof_D[0] + _varst_D[1] * _el_dof_D[1];
        const int indx_glob = el_conn[_el_dof_D[0]-1] + (ivar + _varst_D[0] + _varst_D[1]) * offset;
//         const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
//         assert( kdof_top !=-1);

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
//         printf("loc %d  dof %d ind %d \n",indx_loc_ql,_node_dof[Level][indx_glob],indx_glob);
        assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
            // bc_bd[indx_loc_ql] = _bc[0][el_conn[_el_dof_D[0]-1]];                         // element bc
//             bc_vol[indx_loc] = _bc[0][el_conn[_el_dof_D[0]-1]];                        // element bc
        
        
      }  // end piecewise ------------------------------------------------
  }

  return;
}
// ==============================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void MGSolDA::get_el_dof_bc(
  const int Level,  // level
  const int /*iel*/,    // element subdomain
  //   const int nvars[],       // # of variables to get  <-
//   const int el_nds[],                // # of element nodes for this variable  <-
  const int el_conn[],               // connectivity <-
  const int offset,                  // offset for connectivity <-
  std::vector<int>& el_dof_indices,  // element connectivity ->
//   int bc_vol[],                      // element boundary cond flags ->
   int bc_bd[]                        // element boundary cond flags ->
) const {                          // ==============================================================

  int offset_f=0;
  if(_var_index_D[1]==3)  offset_f=_el_dof_D[0]-1-_el_sides_D;
  for(int id = 0; id < _el_dof_D[0]; id++) {
    // quadratic -------------------------------------------------
    if(id < _el_dof_D[0])
      for(int ivar = 0; ivar < _varst_D[0]; ivar++) {   // ivarq is like idim
//         const int indx_loc = id + ivar * _el_dof_D[0];
        const int indx_loc_ql = id + ivar *  _el_dof_D[0];
        const int indx_glob = el_conn[id] + ivar * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
            assert( kdof_top !=-1);
        //         std::cout << el_dof_indices[indx_loc_ql] << " " << _node_dof[Level][indx_glob] << "\n";
        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
         assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
             bc_bd[indx_loc_ql] = _bc[0][el_conn[id]];                         // element bc
//         bc_vol[indx_loc] = _bc[0][el_conn[id]];                        // element bc
      }  // end quadratic ------------------------------------------------

    //     // linear -----------------------------
    if(id <  _el_dof_D[1])
      for(int ivar = 0; ivar < _varst_D[1]; ivar++) {   // ivarq is like idim
        //        const int  indx_loc_l = id +ivar* _el_dof_D[1];
//         const int indx_loc = id + (ivar + _varst_D[0]) * _el_dof_D[0];
        const int indx_loc_ql = id + ivar *  _el_dof_D[1] + _varst_D[0] *  _el_dof_D[0];
        const int indx_glob = el_conn[id +offset_f] + (ivar + _varst_D[0]) * offset;
//         const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
//      assert( kdof_top !=-1);
        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
         assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                       // element bc
         bc_bd[indx_loc_ql] = _bc[0][el_conn[id ] ];                         // element bc
//         bc_vol[indx_loc] = _bc[0][el_conn[id +offset_f] ];                       // element bc
        
        
      }  // end quadratic ------------------------------------------------

    //     // piecewise -----------------------------
    if(id == _el_dof_D[0]-1)
      for(int ivar = 0; ivar < _varst_D[2]; ivar++) {   // ivarq is like idim
        //        const int  indx_loc_l = id +ivar* el_nds[1];
//         const int indx_loc = 0 + (ivar + _varst_D[0] + _varst_D[1]) * _el_dof_D[0];
        const int indx_loc_ql = 0 + ivar * _el_dof_D[2] + _varst_D[0] * _el_dof_D[0] + _varst_D[1] * _el_dof_D[1];
        const int indx_glob = el_conn[_el_dof_D[0]-1] + (ivar + _varst_D[0] + _varst_D[1]) * offset;
//         const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
//         assert( kdof_top !=-1);

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
//         printf("loc %d  dof %d ind %d \n",indx_loc_ql,_node_dof[Level][indx_glob],indx_glob);
        assert(_node_dof[Level][indx_glob] !=-1);
//         bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
//         bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
            bc_bd[indx_loc_ql] = _bc[0][el_conn[_el_dof_D[0]-1]];                         // element bc
//             bc_vol[indx_loc] = _bc[0][el_conn[_el_dof_D[0]-1]];                        // element bc
        
        
      }  // end piecewise ------------------------------------------------
  }

  return;
}
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void MGSolDA::set_el_dof_bc(
  const int Level,  // level
  const int /*iel*/,    // element subdomain
  //   const int nvars[],       // # of variables to get  <-
//   const int el_nds[],                // # of element nodes for this variable  <-
  const int el_conn[],               // connectivity <-
  const int offset,                  // offset for connectivity <-
  std::vector<int>& el_dof_indices,  // element connectivity ->
  int bc_vol[],                      // element boundary cond flags ->
  int bc_bd[]                        // element boundary cond flags ->
) {                          // ==============================================================

  for(int id = 0; id < _el_dof_D[0]; id++) {
    // quadratic -------------------------------------------------
    if(id < _el_dof_D[0])
      for(int ivar = 0; ivar < _varst_D[0]; ivar++) {   // ivarq is like idim
        const int indx_loc = id + ivar * _el_dof_D[0];
        const int indx_loc_ql = id + ivar * _el_dof_D[0];
        const int indx_glob = el_conn[id] + ivar * offset;
//         const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
//         _bc[1][kdof_top] = bc_bd[indx_loc];                         // element bc
//         _bc[0][kdof_top] = bc_vol[indx_loc];                        // element bc
              _bc[1][el_conn[id]] = bc_bd[indx_loc];                         // element bc
        _bc[0][el_conn[id]] = bc_vol[indx_loc];                        // element bc
        
      }  // end quadratic ------------------------------------------------
    //     // linear -----------------------------
    //       if ( id < el_nds[1] )    for ( int ivar = 0; ivar < _varst_D[1]; ivar++ ) { //ivarq is like idim
    // //        const int  indx_loc_l = id +ivar*el_nds[1];
    //             const int indx_loc = id + ( ivar + _varst_D[0] ) * _el_dof_D[0];
    //             const int indx_loc_ql = id + ivar * el_nds[1] + _varst_D[0] * el_nds[0];
    //             const int indx_glob = el_conn[id] + ( ivar + _varst_D[0] ) * offset;
    //             const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level
    //
    //             el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
    //             bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
    //             bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
    //             } // end quadratic ------------------------------------------------
    //
    //
    //       //     // piecewise -----------------------------
    //       if ( id < el_nds[0] )    for ( int ivar = 0; ivar < _varst_D[2]; ivar++ ) { //ivarq is like idim
    // //        const int  indx_loc_l = id +ivar*el_nds[1];
    //             const int indx_loc = id + ( ivar + _varst_D[0] + _varst_D[1] ) * _el_dof_D[0];
    //             const int indx_loc_ql = id + ivar * el_nds[0] + _varst_D[0] * el_nds[0] + _varst_D[1] *
    //             el_nds[1]; const int indx_glob = id + iel * el_nds[2] + ( ivar + _varst_D[0] + _varst_D[1] ) *
    //             offset; const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level
    //
    //             el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
    //             bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
    //             bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
    //             } // end piecewise ------------------------------------------------
  }

  return;
}

///---------USER SOURCE---------------
// ============================================================================
/// This function reads Boundary conditions  from function
// void MGSolDA::bc_read_D(
//     int /*face_id_node*/,  ///<  face identity          (in)
//     int /*mat_flag*/,      ///<  volume identity         (in)
//     double /*xp*/[],       ///< xp[] node coordinates    (in)
//     int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
//     int bc_flag[]          ///< boundary condition flag  (out)
// ) {                        // =========================================================================
//
//   bc_Neum[0] = 0;  bc_flag[0] = 0;
// }

// ===========================================================================================
// void MGSolDA::ic_read(
//     int face_id_node,
//     int mat_id_elem,  // If using GAMBIT_INTERFACE, mat_flag given in Gambit for this node
//     double xp[],      // Coordinates of the point
//     int iel,          // element
//     double ic[]           // Initial value of the [ivar] variable of the system
// ) {
//   for (int ivar = 0; ivar < _n_vars; ivar++) { ic[ivar] = 10.*xp[0]; }
// }

// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
