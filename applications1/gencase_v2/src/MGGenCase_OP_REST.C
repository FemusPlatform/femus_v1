// std libraries ------------------------------------------
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// LibMesh includes ---------------------------------------
#include "gencase_conf.h"
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
/// This function computes the restrictor operator
void MGGenCase::compute_rest(
    int* n_nodes_lev,                      // node level vector
    int max_elnd,                          // max # of elements at each point
    int* off_el,                           // offset element vectors
    std::vector<std::pair<int, int>> v_el    , // element  ordering
    int** elem_sto,                 // storage element structure
    int* v_inv_nd,                  // node inverse ordering
    int* v_inv_el, int** g_indexL,  //  map nodes (level)
    int* min_nd[],                   //  node offset vector
    int* max_nd[]                   //  node offset vector
) {
  // =====================================================
  // hdf5 print Restrictor    data_in/Rest*.h5
  std::string input_dir = _mgutils._inout_dir;
  std::string f_rest = _mgutils.get_file("F_REST");

//   int level_row[3];
//   int level_cln[3];  // level row x cln
  int nodes_row[3];
  int nodes_cln[3];  // nodes row x cln
  int el_nodes[3];   // nodes for elements

  // level loop
  for (int Level1 = 0; Level1 < _n_levels - 1; Level1++) {
    // file name and creation  -----------------
    std::ostringstream name;
    name << input_dir << f_rest << Level1 + 1 << "_" << Level1 << ".h5";
    hid_t fileR = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    int dofType[4] = {3, 1, 0, 2}; // gindex_l: 0-> q, 1->l, 2->f, 3->c
//     // Set up --------------------------
//     nodes_row[0] = n_nodes_lev[Level1];
//     nodes_cln[0] = n_nodes_lev[Level1+1];  // quadratic
//     nodes_row[1] = n_nodes_lev[_n_levels+Level1];
//     nodes_cln[1] = n_nodes_lev[_n_levels+Level1+1];
//     nodes_row[2] = n_nodes_lev[3*_n_levels+Level1];
//     nodes_cln[2] = n_nodes_lev[3*_n_levels+Level1+1];
    
    
//     // levels
//     level_row[0] = Level1;
//     level_cln[0] = Level1 + 1;  // quadratic
//     level_row[1] = _n_levels;
//     if (Level1 > 0) { level_row[1] = Level1 - 1; }
//     level_cln[1] = Level1;  // linear
//     level_row[2] = Level1;
//     level_cln[2] = Level1 + 1;  // piecewise
//     // nodes
//     nodes_row[0] = n_nodes_lev[Level1];
//     nodes_cln[0] = n_nodes_lev[Level1 + 1];  // quadratic
//     nodes_row[1] = n_nodes_lev[2 * _n_levels];
//     if (Level1 > 0) {
//       nodes_row[1] = n_nodes_lev[Level1 - 1];  // linear
//     }
//     nodes_cln[1] = n_nodes_lev[Level1];
//     int sum_lev = 0;
//     int sum_levp1 = 0;
//     for (int jpr = 0; jpr < _n_subdomains; jpr++) {
//       sum_levp1 += off_el[jpr * _n_levels + Level1 + 2] - off_el[jpr * _n_levels + Level1 + 1];
//       sum_lev += off_el[jpr * _n_levels + Level1 + 1] - off_el[jpr * _n_levels + Level1];
//     }
//     nodes_row[2] = sum_lev;
//     nodes_cln[2] = sum_levp1;  // piecewise
    // element nodes
//     el_nodes[0] = NDOF_FEM;
//     el_nodes[1] = NDOF_P;
//     el_nodes[2] = NDOF_K;

    // matrix structures --------------------------
    int* Rest_q;       // pos
    double* values_q;  // values
    int* mem;          // mem row
    int* len_q;        // # of diagonal row entries
    int* len_qoff;     // # of offdiagonal row entries

    for (int iql = 0; iql <= 3; iql++) {
      // [0]=center [1]=vert [2]=quad [3] face
      int levelCoarse = dofType[iql]*_n_levels + Level1 ;  
      int levelFine = levelCoarse+1;
      int nodes_row = n_nodes_lev[levelCoarse];
      int nodes_cln = n_nodes_lev[levelFine];  
        
        
        
      int count_q = 0;
//       if (iql < 2) {
        // set up
        Rest_q = new int[nodes_row * NDOF_FEM * max_elnd];
        values_q = new double[nodes_row * NDOF_FEM * max_elnd];
        mem = new int[nodes_row * NDOF_FEM * max_elnd];
        len_qoff = new int[nodes_row + 1];
        len_q = new int[nodes_row + 1];

        //  Restrictor quadratic and linear--------------------------------------
//         count_q = compute_res_qq(iql,
//             Level1, Level1+iql*_n_levels, Level1+iql*_n_levels+1, el_nodes[iql], nodes_row[iql], nodes_cln[iql], max_elnd,
//             off_el, v_el, elem_sto, v_inv_nd, g_indexL, min_nd[iql],max_nd[iql], mem, values_q, Rest_q, len_q, len_qoff
//         );
        
           count_q = compute_res_qq(iql,
            Level1,  levelCoarse, levelFine, nodes_row, nodes_cln, max_elnd,
            off_el, v_el, elem_sto, v_inv_nd, g_indexL, min_nd[dofType[iql]],max_nd[dofType[iql]], mem, values_q, Rest_q, len_q, len_qoff
        ); 
           assert(count_q<= nodes_row * NDOF_FEM * max_elnd);
//       } else {
//                // compute quadratic and linear prolongator
//         count_q = compute_prol_qq(iql,
//             Level1,levelFine, levelCoarse, nodes_row, off_el, v_el, elem_sto,
//             v_inv_nd, g_indexL, min_nd[dofType[iql]],max_nd[dofType[iql]], values_q, Prol_q, len_q, len_qoff);
//         // set up
//         max_elnd =
//             elem_sto[0][NDOF_FEM + 4];  //<-pay attention to this instruction, do first quad and linear!!
//         Rest_q = new int[nodes_row[iql] * max_elnd];
//         values_q = new double[nodes_row[iql] * max_elnd];
//         mem = new int[nodes_row[iql] * max_elnd];
//         len_qoff = new int[nodes_row[iql] + 1];
//         len_q = new int[nodes_row[iql] + 1];
// 
//         //  Restrictor constant --------------------------------------
//         count_q = compute_res_kk(
//             Level1, Level1+3*_n_levels, Level1+3*_n_levels+1, el_nodes[iql], nodes_row[iql], nodes_cln[iql], max_elnd,
//             off_el, v_el, elem_sto, v_inv_nd, v_inv_el, g_indexL, min_nd[3],max_nd[3], mem, values_q, Rest_q, len_q,
//             len_qoff);
//       }
      // print ----------------------------------------------------
      print_op_h5(
          name.str(), nodes_row, nodes_cln, count_q, Rest_q, values_q, len_q, len_qoff, iql);
      // clean ----------------------------------------------------
      delete[] Rest_q;
      delete[] values_q;
      delete[] mem;  // Restrictor structures
      delete[] len_qoff;
      delete[] len_q;  // compressed row length
    }

    H5Fclose(fileR);

  }  // ----- end Level loop  ------

  return;
}

// =======================================================
/// This function computes the quadratic matrix
int MGGenCase::compute_res_qq(int name_label,
    int Level1,                             // level
    int Level_row,                          // row level
    int Level_cln,                          // colum level
    int n_nodes_row,                        // # of row nodes
    int n_nodes_cln,                        // # of column nodes
    int max_elnd,                           // max # of elements at each point
    int* off_el,                            // offset elements vector
    std::vector<std::pair<int, int>> v_el,  // element ordering
    int** elem_sto,                         // element storage
    int* v_inv_nd,                          // node ordering
    int** g_indexL,                         // map nodes
    int* min_nd,                            // offset nodes vector
    int* max_nd,                            // offset nodes vector
    int mem[],                              //  memory index
    double values_q[],                      // restriction nonzero values
    int Rest_q[],                           // restriction compressed pos
    int len_q[],                            // # of diagonal nonzero row entries
    int len_qoff[]                          // # of offdiagonal nonzero row entries
) {
  // ====================================================

  // zeroing the structures
  for (int i = 0; i < n_nodes_row * NDOF_FEM * max_elnd; i++) {
    Rest_q[i] = -1;
    values_q[i] = 0.;
  }
  for (int im = 0; im < n_nodes_row; im++) { mem[im] = 0; }
  
   int kk_max=NDOF_FEM;int kj_max=NDOF_FEM;   int kk_min=0;int kj_min=0;
   switch(name_label){
       case 0:      kk_min=NDOF_FEM-1; kk_max=NDOF_FEM;     break;
       case 1:      kk_min=0;   kk_max=NDOF_P;            break;
       case 2:      kk_min=0;   kk_max=NDOF_FEM;          break;
       case 3:     kk_min=NDOF_FEM-1-_n_local_faces;   kk_max=NDOF_FEM-1;       break;     
       default: abort();
   }
  // Restriction operator
  for (int pr = 0; pr < _n_subdomains; pr++) {
    for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
      int el_lib = v_el[iel].second;
      int n_childs = elem_sto[el_lib][NDOF_FEM + 4];  //(*it_n)->n_children();
      for (int i_ch = 0; i_ch < n_childs; i_ch++) {
        int ch_el = elem_sto[el_lib][NDOF_FEM + 5 + i_ch];  // Elem* child=(*it_n)->child(i_ch);
         for (int kk =0; kk < kk_max-kk_min; kk++) {
//         for (int kk = 0; kk < ndof_row_ql; kk++) {
          int op = v_inv_nd[elem_sto[el_lib][1 + kk+kk_min]];
          int irow = g_indexL[Level_row][op];
          
          for (int kj =0; kj < kk_max-kk_min; kj++) {  
//           for (int kj = 0; kj < ndof_row_ql; kj++) {
            int ch_op = v_inv_nd[elem_sto[ch_el][1 + kj+kk_min]];

// #ifdef REST_SIMPLE
//             if (op == ch_op) {
//               Rest_q[irow * NDOF_FEM * max_elnd] = g_indexL[Level_cln][ch_op];
//               values_q[irow * NDOF_FEM * max_elnd] = NDOF_P;
//               if (ndof_row_ql == NDOF_P) { values_q[irow * NDOF_FEM * max_elnd] = NDOF_P; }
//               mem[irow] = 1;
//             }
//                      switch(name_label){
//        case 0:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] =1. ;    break;
//        case 1:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_l[ch_id][kk][kj];           break;
//        case 2:     values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_q[ch_id][kk][kj]; ;          break;
//         case 3:    values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj + kk_min] = _geomel._embedding_matrix_f[ch_id][kk][kj];      break;     
//        default: abort();  
            
// #else
            double val_qq = 0;
            switch(name_label){
        case 0:     val_qq  =1. ;    break;
        case 1:     val_qq = _geomel._embedding_matrix_l[i_ch][kj][kk];           break;
        case 2:     val_qq = _geomel._embedding_matrix_q[i_ch][kj][kk]; ;          break;
#if  ELTYPE == 27      
        case 3:     val_qq = n_childs*0.5*_geomel._embedding_matrix_f[i_ch][kj][kk];      break;
#endif        
        default: abort(); 
            }            
            
            
            
//             if (ndof_row_ql == NDOF_P) {
//               val_qq = _geomel._embedding_matrix_l[i_ch][kj][kk];
//             } else {
//               val_qq = _geomel._embedding_matrix_q[i_ch][kj][kk];
//             }
            if (fabs(val_qq) > 1.e-6) {
              int jcln = g_indexL[Level_cln][ch_op];
              int indx0 = irow * NDOF_FEM * max_elnd;
              int indx_mem = mem[irow];
              int flag = 0;
              for (int ks = 0; ks < indx_mem; ks++)  if (Rest_q[indx0 + ks] == jcln) { flag = 1; 
               if( name_label==3 && fabs(val_qq)>0.3*n_childs*0.5) values_q[indx0 + ks] = val_qq;
            }
// 
               if (flag == 0)  if( name_label!=3 || fabs(val_qq)>0.3*n_childs*0.5) {
                 Rest_q[indx0 + indx_mem] = jcln;
                 values_q[indx0 + indx_mem] = val_qq;
                 mem[irow]++;
                 //printf("irow %d  jcln %d  val_qq %g ",irow,  jcln ,  val_qq );
               }
            }
// #endif
          }
        }
      }
       //printf(" \n"); std::cout << std::flush;
    }
     
  }
int _off_lev_fl[4]={1,0,0,1};
  // compressing zeros
  int* offPc = new int[_n_subdomains + 1];
  offPc[0] = 0;
  for (int idom = 0; idom < _n_subdomains; idom++) {
    offPc[idom + 1] = offPc[idom] + (max_nd[idom * _n_levels + Level1 + 1] - min_nd[idom * _n_levels+_off_lev_fl[name_label]*(Level1 + 1)]);
  }
  offPc[_n_subdomains] += 1;
  //  Compressing quadratic restrictor
  int count_q = 0;
  len_q[0] = count_q;
  int count_qoff = 0;
  len_qoff[0] = count_qoff;
  int ki = 0;
  for (int kdom = 0; kdom < _n_subdomains; kdom++) {
    for (int iki = 0; iki < max_nd[kdom * _n_levels + Level1] - min_nd[kdom * _n_levels+_off_lev_fl[name_label]*Level1]; iki++) {
      // computing various flags and lengths
      for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
        int kk = NDOF_FEM * max_elnd * ki + kj;
        if (fabs(values_q[kk]) > 1.e-8) {
          Rest_q[count_q] = Rest_q[kk];
          values_q[count_q] = values_q[kk];
          count_q++;
          if (Rest_q[kk] >= (int)offPc[kdom + 1] || Rest_q[kk] < (int)offPc[kdom]) { count_qoff++; }
        }
      }
      len_q[ki + 1] = count_q;
      len_qoff[ki + 1] = count_qoff;
      ki++;
    }
  }

  delete[] offPc;
  return count_q;
}

// =======================================================
/// This function computes the quadratic matrix
int MGGenCase::compute_res_kk(
    int Level1,                             // level
    int Level_row,                          // row level
    int Level_cln,                          // colum level
    int ndof_row_ql,                        // element dofs
    int n_nodes_row,                        // # of row nodes
    int n_nodes_cln,                        // # of column nodes
    int max_elnd,                           // max # of elements at each point
    int* off_el,                            // offset elements vector
    std::vector<std::pair<int, int>> v_el,  // element ordering
    int** elem_sto,                         // element storage
    int* v_inv_nd,                          // node ordering
    int* v_inv_el, int** g_indexL,          // map nodes
    int* min_nd,                            // offset nodes vector
    int* max_nd,                            // offset nodes vector
    int mem[],                              //  memory index
    double values_q[],                      // restriction nonzero values
    int Rest_q[],                           // restriction compressed pos
    int len_q[],                            // # of diagonal nonzero row entries
    int len_qoff[]                          // # of offdiagonal nonzero row entries
) {
  // ====================================================

  // zeroing the structures
  for (int i = 0; i < n_nodes_row * max_elnd; i++) {
    Rest_q[i] = -1;
    values_q[i] = 0.;
  }
  for (int im = 0; im < n_nodes_row; im++) { mem[im] = 0; }

  // Restriction operator
  int sum_lev = 0;
  int count_q = 0;
  len_q[0] = 0;
  len_qoff[0] = 0;
  for (int pr = 0; pr < _n_subdomains; pr++) {
    int delta = off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
    for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
      int el_lib = v_el[iel].second;
      int n_childs = elem_sto[el_lib][NDOF_FEM + 4];  //(*it_n)->n_children();
      for (int kk = 0; kk < ndof_row_ql; kk++) {
        int i_row = kk + (iel - off_el[pr * _n_levels + Level1] + sum_lev) * ndof_row_ql;
        int ind_row = i_row * n_childs;
        for (int i_ch = 0; i_ch < n_childs; i_ch++) {
          int ch_el = elem_sto[el_lib][NDOF_FEM + 5 + i_ch];  // Elem* child=(*it_n)->child(i_ch);
          int j_par = v_inv_el[ch_el];
          int sum_levp1 = 0;
          for (int jpr = 0; jpr < pr; jpr++) {
            sum_levp1 += off_el[jpr * _n_levels + Level1 + 2] - off_el[jpr * _n_levels + Level1 + 1];
          }
          int j_column = kk + (j_par - off_el[pr * _n_levels + Level1 + 1] + sum_levp1) * ndof_row_ql;
          double val_qq = 1.;
          Rest_q[ind_row + i_ch] = j_column;
          values_q[ind_row + i_ch] = val_qq;
          count_q++;
        }
        len_q[i_row + 1] = count_q;
        len_qoff[i_row + 1] = 0;
      }
    }
    sum_lev += delta;
  }
  return count_q;
}



// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
