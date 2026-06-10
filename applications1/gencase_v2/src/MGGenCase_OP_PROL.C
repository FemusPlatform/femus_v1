// std libraries ------------------------------------------
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// LibMesh includes ---------------------------------------
// #include "MGFE_conf.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// configure includes -------------------------------------
// #include "Domain_conf.h"  //  domain dimensions
#include "gencase_conf.h"    //
#include "MGGenCase.h"
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"

// MED includes -------------------------------------------
#ifdef HAVE_MED
#include "MEDCoupling.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#endif

// local includes -----------------------------------------
#include "MGGenCase_conf.h"
#include "MGGeomEl.h"
#include "MGUtils.h"






// =====================================================================
/// This function computes the prologator operator
void MGGenCase::compute_prol(
    int* n_nodes_lev, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
    int* v_inv_el, int** g_indexL, int* min_nd[], int* max_nd[]) {
  // name
  std::string input_dir = _mgutils._inout_dir;
  std::string f_prol = _mgutils.get_file("F_PROL");
  
  

//   int nodes_row[3];
//   int nodes_cln[3];  // nodes row x cln
//   int levelFine[3];
//   int levelCoarse[3];

  for (int Level1 = 1; Level1 < _n_levels; Level1++) {
    // Set up --------------------------
    //  file  -------------------------------------------------------
    std::ostringstream name;
    name << input_dir << f_prol << Level1 - 1 << "_" << Level1 << ".h5";
    // Create file (data_in/Prol*.h5)
    hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // set up ----------------------------------------
    int dofType[4] = {3, 1, 0, 2}; // gindex_l: 0-> q, 1->l, 2->f, 3->c

    // matrix structures --------------------------
    int* Prol_q;       // compressed sparse prolongation
    double* values_q;  // compressed prolongation values
    int* len_q;        // # of diagonal row entries
    int* len_qoff;     // # of offdiagonal row entries

    // Prolongation Operator   (quad and linear)
    for (int iql = 0; iql <= 3; iql++) {   
        
        int levelFine = dofType[iql]*_n_levels + Level1;
        int levelCoarse = dofType[iql]*_n_levels + Level1 - 1;
        int nodes_row = n_nodes_lev[levelFine ];
        int nodes_cln = n_nodes_lev[ levelCoarse];
    
      //[0]=quad  [1]=linear
      int count_q = 0;
//       if (iql < 2) {
        len_q = new int[nodes_row + 1];
        len_qoff = new int[nodes_row + 1];
        Prol_q = new int[nodes_row * NDOF_FEM];
        values_q = new double[nodes_row * NDOF_FEM];

        // compute quadratic and linear prolongator
        count_q = compute_prol_qq(iql,
            Level1,levelFine, levelCoarse, nodes_row, off_el, v_el, elem_sto,
            v_inv_nd, g_indexL, min_nd[dofType[iql]],max_nd[dofType[iql]], values_q, Prol_q, len_q, len_qoff);
        
        assert(count_q <= nodes_row * NDOF_FEM);

//       } else {
//         len_q = new int[nodes_row[iql] + 1];
//         len_qoff = new int[nodes_row[iql] + 1];
//         Prol_q = new int[nodes_row[iql] /* *NDOF_K */];
//         values_q = new double[nodes_row[iql] /* *NDOF_K */];
// 
//         // compute piecewise prolongator
//         count_q = compute_prol_kk(
//             Level1,Level1+3*_n_levels, Level1+3*_n_levels-1, el_nodes[iql], nodes_row[iql], off_el, v_el, elem_sto,
//             v_inv_nd, v_inv_el, g_indexL, min_nd[3],max_nd[3], values_q, Prol_q, len_q, len_qoff);
//       }

      print_op_h5(
          name.str(), nodes_row, nodes_cln, count_q, Prol_q, values_q, len_q, len_qoff, iql);
      // clean
      delete[] Prol_q;
      delete[] values_q;
      delete[] len_q;
      delete[] len_qoff;
    }

    H5Fclose(fileP);

  }  // end Level1

  return;
}


// =======================================================
/// This function computes the quadratic matrix
int MGGenCase::compute_prol_qq(int name_label,
    int Level1,                             // level
    int Level_row,                          // row level
    int Level_clmn,                         // column level
    int n_nodes_f, int* off_el,             // offset elements vector
    std::vector<std::pair<int, int>> v_el,  // element ordering
    int** elem_sto,                         // element storage
    int* v_inv_nd,                          // node ordering
    int** g_indexL,                         // map nodes
    int* min_nd,                            // offset nodes vector
    int* max_nd,                            // offset nodes vector
    double values_q[],                      // prolongation nonzero values
    int Prol_q[],                           // prolongation compressed pos
    int len_q[],                            // # of diagonal nonzero row entries
    int len_qoff[]                          // # of offdiagonal nonzero row entries
) {
  // ======================================================

  // zeroing matrix
//   int *P_signf=new int[81]; int *P_signc=new int[81]; 
  for (int im = 0; im < n_nodes_f * NDOF_FEM; im++) { Prol_q[im] = -1;   values_q[im] = 0.;  }
//    for (int im = 0; im < n_nodes_f ; im++) { P_signf[im]=-1; P_signc[im]=-1; }
  
  // generation matrix entries
   int kk_max=NDOF_FEM; int kj_max=NDOF_FEM;   int kk_min=0;int kj_min=0;
   switch(name_label){
       case 0:      kk_min=NDOF_FEM-1; kk_max=NDOF_FEM;     break;
       case 1:      kk_min=0;   kk_max=NDOF_P;            break;
       case 2:      kk_min=0;   kk_max=NDOF_FEM;          break;
       case 3:     kk_min=NDOF_FEM-1-_n_local_faces;   kk_max=NDOF_FEM-1;       break;     
       default: abort();
   }
 
  
  // Prolongation Operator linear ---------------------------------------------------
  for (int pr = 0; pr < _n_subdomains; pr++) {
    for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
      int el_lib = v_el[iel].second;// element levelFine
      int par = elem_sto[el_lib][NDOF_FEM + 3];// processor
      int ch = 0; 
      int ch_id = 0;  //child
      while (el_lib != elem_sto[par][NDOF_FEM + 5 + ch]) { ch_id = ++ch; }
      for (int kk =0; kk < kk_max-kk_min; kk++) { 
        int el_n = v_inv_nd[elem_sto[el_lib][1 + kk + kk_min]];
          for (int kj =0; kj < kk_max-kk_min; kj++) { 
          int par_n = v_inv_nd[elem_sto[par][1 + kj + kk_min]];
          Prol_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = g_indexL[Level_clmn][par_n];
             switch(name_label){
       case 0:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] =1. ;    break;
       case 1:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_l[ch_id][kk][kj];           break;
       case 2:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_q[ch_id][kk][kj];          break;
#if ELTYPE == 27          
        case 3:   
              values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_f[ch_id][kk][kj];         
//             int el_f1=el_lib-4;
#endif 
//         if(el_lib > P_signf[el_n])  P_signf[el_n]=el_lib;
             
        break;     
       default: abort();
   }
          
//           if (ndof_row_ql == NDOF_P)
//             values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj] = _geomel._embedding_matrix_l[ch_id][kk][kj];
//           else
//             values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj] = _geomel._embedding_matrix_q[ch_id][kk][kj];
        }
      }
    }
  }
 
int _off_lev_fl[4]={1,0,0,1};
// Compressing prolongation ---------------------------
  
  int* offPc = new int[_n_subdomains + 1];
  offPc[0] = 0;
  for (int idom = 0; idom < _n_subdomains; idom++) {
    offPc[idom + 1] = offPc[idom] + (max_nd[idom * _n_levels + Level1-1] - min_nd[idom * _n_levels+_off_lev_fl[name_label]*(Level1 - 1)]);
  }
  offPc[_n_subdomains] += 1;

  
  // Compressing prolongation ---------------------------
  int count_q = 0;
  len_q[0] = count_q;
  int count_qoff = 0;
  len_qoff[0] = count_qoff;

  int ki = 0; int count_ki=0;
  for (int kdom = 0; kdom < _n_subdomains; kdom++) {
    for (int iki = 0; iki < max_nd[kdom * _n_levels + Level1] - min_nd[kdom * _n_levels + _off_lev_fl[name_label]*Level1]; iki++) {
      for (int kj = 0; kj < NDOF_FEM; kj++) {
        int kk = NDOF_FEM * ki + kj;
        if (fabs(values_q[kk]) > 1.e-8) {
          Prol_q[count_q] = Prol_q[kk];
          values_q[count_q] = values_q[kk];
          count_q++;
          if (Prol_q[kk] >= (int)offPc[kdom + 1] || Prol_q[kk] < (int)offPc[kdom]) { count_qoff++; }
        }
      }
//       len_q[ki + 1] = count_q;
//       len_qoff[ki + 1] = count_qoff;
      ki++;
      if(count_q-len_q[count_ki]>0) count_ki++;
      len_q[count_ki] = count_q;
      len_qoff[count_ki] = count_qoff;
    }
  }
//  for(int i=0;i<81;i++)    std::cout << i << "  " << P_signf[i] << std::endl;
  // clean -------------------------------------------------
  delete[] offPc;
  return count_q;
}

// // =======================================================
// /// This function computes the quadratic matrix
// int MGGenCase::compute_prol_kk(
//     int Level1,                             // level
//     int Level_row,                          // row level
//     int Level_clmn,                         // column level
//     int ndof_row_ql,                        // element dofs
//     int n_nodes_f, int* off_el,             // offset elements vector
//     std::vector<std::pair<int, int>> v_el,  // element ordering
//     int** elem_sto,                         // element storage
//     int* v_inv_nd,                          // node ordering
//     int* v_inv_el, int** g_indexL,          // map nodes
//     int* min_nd,                            // offset nodes vector
//      int* max_nd,                            // offset nodes vector
//    double values_q[],                      // prolongation nonzero values
//     int Prol_q[],                           // prolongation compressed pos
//     int len_q[],                            // # of diagonal nonzero row entries
//     int len_qoff[]                          // # of offdiagonal nonzero row entries
// ) {
//   // ======================================================
// 
//   // zeroing matrix
//   for (int im = 0; im < Level_row * 1; im++) {
//     Prol_q[im] = -1;
//     values_q[im] = 0.;
//   }
// 
//   // Simple !!! Prolongation Operator linear ---------------------------------------------------
//   int sum_lev = 0;
//   int count_q = 0;
//   len_q[0] = 0;
//   len_qoff[0] = 0;
//   for (int pr = 0; pr < _n_subdomains; pr++) {
//     int delta = off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
//     for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
//       int el_lib = v_el[iel].second;
//       int i_row = iel - off_el[pr * _n_levels + Level1] + sum_lev;
//       int par = elem_sto[el_lib][NDOF_FEM + 3];
//       int j_par = v_inv_el[par];
//       int sum_levm1 = 0;
//       for (int jpr = 0; jpr < pr; jpr++) {
//         sum_levm1 += off_el[jpr * _n_levels + Level1] - off_el[jpr * _n_levels + Level1 - 1];
//       }
//       for (int kk = 0; kk < ndof_row_ql; kk++) {
//         int ind_row = i_row * ndof_row_ql + kk;
//         int j_column = (j_par - off_el[pr * _n_levels + Level1 - 1] + sum_levm1) * ndof_row_ql + kk;
//         Prol_q[ind_row] = j_column;
//         values_q[ind_row] = 1.;
//         len_q[ind_row + 1] = count_q + 1;
//         len_qoff[ind_row + 1] = 0;
//         count_q++;
//       }
//     }
//     sum_lev += delta;
//   }
//   return count_q;
// }



// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
