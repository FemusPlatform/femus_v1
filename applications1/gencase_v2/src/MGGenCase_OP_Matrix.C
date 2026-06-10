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






//-----------------------------------------------------------------
/// This function computes the sparse matrix patern
void MGGenCase::compute_matrix(
    int* n_nodes_lev,                       // nodes for level
    int max_elnd,                           // max # of elements at each point
    int* off_el,                            // element offset
    std::vector<std::pair<int, int>> v_el,  // element ordering
    int** elem_sto,                         // element structure
    int* v_inv_nd,                          // node inverse ordering
    int** g_indexL,                         // node map
    int* min_nd[],                           // node offset
    int* max_nd[],
    int** elem_conn                         // neighbours  
) {
  // ====================================================

  std::cout << " Computing matrix :\n";

  // structures
  int n_nodes = n_nodes_lev[_n_levels - 1];   // quadratic nodes
  int n_nodesl = n_nodes_lev[_n_levels - 1];  // linear nodes
  int* len_q = new int[n_nodes + 1];          // row length
  int* len_qoff = new int[n_nodes + 1];       // row offdiagonal length
  int* mem = new int[n_nodes];                // row memorization index

  int* Mat_q = new int[n_nodes * NDOF_FEM * max_elnd];
  int* aux = new int[NDOF_FEM * max_elnd * n_nodes];

  // directory
  std::string input_dir = _mgutils._inout_dir;
  std::string f_matrix = _mgutils.get_file("F_MATRIX");

  int nodes_top[4];  // nodes row x cln
  

  for (int Level1 = 0; Level1 < _n_levels; Level1++) {
    std::ostringstream name;
    name << input_dir << f_matrix << Level1 << ".h5";
    hid_t fileM = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    std::cout << " Level =" << Level1 << " \n ";       // quad level print
    n_nodes = n_nodes_lev[Level1];  //  nodes (quad)
    // quad matrix dimensions
    // Set up --------------------------
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t start_time = std::clock();
#endif
   
int* offset_sub = new int[_n_subdomains + 1]; offset_sub[0] = 0;
      for (int idom = 0; idom < _n_subdomains; idom++) {
        offset_sub[idom + 1] = offset_sub[idom] + (max_nd[0][idom * _n_levels + Level1] - min_nd[0][idom * _n_levels]);
      }
 
#ifdef PRINT_INFO  // ---------------------------------------
    printf(
        " \n Off processor limits quad and linear for level = %d and proc %d \n ---------------------------  "
        "\n",
        Level1, _n_subdomains);
    for (int idom = 0; idom <= _n_subdomains; idom++) {
      printf("  %d %d \n", idom, offset_sub[idom]);
    }
#endif
    offset_sub[_n_subdomains] += 1;// for interval bound

    // nodes [0]=quad  [1]=linear
    nodes_top[0] = n_nodes_lev[3*_n_levels+ Level1];
    nodes_top[1] = n_nodes_lev[_n_levels+ Level1];  //linear
    nodes_top[2] = n_nodes_lev[ Level1];   // quadratic
    nodes_top[3] = n_nodes_lev[2*_n_levels+ Level1];

    // Generating matrix operators ----------------------------------------

    // matrices qq-ql-lq-ll  ===============================================
    for (int iql = 0; iql <=3; iql++) {
      for (int jql = 0; jql <=3; jql++) {
        int name_label = iql * 10 +  jql;
        int count_q;

        if ((iql == 1) && (jql == 1)) {
          int  max_elnd_ghost = max_elnd + max_elnd*elem_conn[0][0];
          int* Mat_q_ghost = new int[n_nodes * NDOF_FEM * max_elnd_ghost];
          int* aux_ghost = new int[NDOF_FEM * max_elnd_ghost * n_nodes];

          count_q = compute_mat_qq_ghost(
            name_label,Level1, n_nodes, max_elnd_ghost, 
            off_el,
            v_el, elem_sto, v_inv_nd, g_indexL, min_nd[0],max_nd[0], offset_sub, aux_ghost, mem, Mat_q_ghost, len_q, len_qoff,
            elem_conn);

          print_op_h5(
            name.str(), nodes_top[iql], nodes_top[jql], count_q, Mat_q_ghost, NULL, len_q, len_qoff, name_label);
        }
        else {
          count_q = compute_mat_qq(
            name_label,Level1, n_nodes, max_elnd, 
            off_el,
            v_el, elem_sto, v_inv_nd, g_indexL, min_nd[0],max_nd[0], offset_sub, aux, mem, Mat_q, len_q, len_qoff);

          print_op_h5(
            name.str(), nodes_top[iql], nodes_top[jql], count_q, Mat_q, NULL, len_q, len_qoff, name_label);
        }
           
        
      }
    }

#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t end_time = std::clock();
    std::cout << " quad Matrix compute time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif
    delete[] offset_sub;
    H5Fclose(fileM);
  }

  // cleaning -------------------------------------------------
  delete[] len_q;  delete[] len_qoff;
  delete[] mem;  delete[] Mat_q;
  delete[] aux;

  return;
}

// =======================================================
/// This function computes the quadratic matrix
int MGGenCase::compute_mat_qq(int name_label,
    int Level1,        // level
    int n_nodes,       // top nodes
    int max_elnd, 
    int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, 
    int* v_inv_nd, int** g_indexL, int* min_nd, int* max_nd, 
    int offset_sub[], int aux[], int mem[],  // number elemets row
    int Mat_q[], int len_q[], int len_qoff[]
) {
  // ================================================

  // zeroing matrix
  for (int im = 0; im < n_nodes; im++) { mem[im] = 0; }
  for (int im = 0; im < n_nodes * NDOF_FEM * max_elnd; im++) { Mat_q[im] = -1; }

  // generation matrix entries
   int kk_max=NDOF_FEM;int kj_max=NDOF_FEM;
   int kk_min=0;int kj_min=0;
   switch(name_label/10){
       case 0:
           kk_max=NDOF_FEM; kk_min=NDOF_FEM-1;
           break;
       case 1:
           kk_min=0;   kk_max=NDOF_P;
            break;
       case 2:
           kk_min=0;   kk_max=NDOF_FEM;
            break;
        case 3:
           kk_min=NDOF_FEM-1-_n_local_faces;   kk_max=NDOF_FEM-1;
            break;     
       default: abort();
   }
      switch(name_label%10){
       case 0:
           kj_max=NDOF_FEM; kj_min=NDOF_FEM-1;
           break;
       case 1:
           kj_min=0;   kj_max=NDOF_P;
            break;
       case 2:
           kj_min=0;   kj_max=NDOF_FEM;
            break;
        case 3:
           kj_min=NDOF_FEM-1-_n_local_faces;   kj_max=NDOF_FEM-1;
            break;     
       default: abort();
   }
 

  for (int pr = 0; pr < _n_subdomains; pr++) {
    for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
      int el_lib = v_el[iel].second;
              for (int kk = kk_min; kk <kk_max; kk++) {
        int op = v_inv_nd[elem_sto[el_lib][1 + kk]];
        int i_row = g_indexL[Level1][op];
       
        for (int kj = kj_min; kj <kj_max ; kj++)  {
          int opj = v_inv_nd[elem_sto[el_lib][1 + kj]];
          int ind_row = i_row * NDOF_FEM * max_elnd;
     
             Mat_q[ind_row + mem[i_row] + kj] = g_indexL[Level1][opj];
        }
        mem[i_row] += NDOF_FEM;
      }
    }
  }  
  // compressing zeros and counting -------------------------------------
  // counters
  int count_q = 0;     len_q[0] = 0;
  int count_qoff = 0;  len_qoff[0] = 0;
  int ki = 0; int count_ki=0;
  // zeroing aux
  for (int i = 0; i < NDOF_FEM * max_elnd * n_nodes; i++) { aux[i] = 0; }
  // compressing
  for (int kdom = 0; kdom < _n_subdomains; kdom++) {
    for (int iki = 0; iki < (int)(max_nd[kdom * _n_levels + Level1] - min_nd[kdom * _n_levels + Level1]); iki++) {
      // setting aux
      for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
        int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
        if (knod_q >= 0) { aux[knod_q] = 1; }
      }
      // computing various flags and lengths
      for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
        int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
        if (knod_q >= 0 && aux[knod_q] == 1) {
          aux[knod_q] = 0;
          Mat_q[count_q] = knod_q;
          count_q++;
          if (knod_q < (int)offset_sub[kdom] || knod_q >= (int)offset_sub[kdom + 1]) { count_qoff++; }
        }
      }
 
      ki++;
      if(count_q-len_q[count_ki]>0) count_ki++;
      len_q[count_ki] = count_q;
      len_qoff[count_ki] = count_qoff;
      
    }
  }

  return count_q;
}

// =======================================================
/// This function computes the quadratic matrix
int MGGenCase::compute_mat_qq_ghost(int name_label,
    int Level1,        // level
    int n_nodes,       // top nodes
    int max_elnd, 
    int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, 
    int* v_inv_nd, int** g_indexL, int* min_nd, int* max_nd, 
    int offset_sub[], int aux[], int mem[],  // number elemets row
    int Mat_q[], int len_q[], int len_qoff[],
    int** elem_conn
) {
  // ================================================

  // zeroing matrix
  for (int im = 0; im < n_nodes; im++) { mem[im] = 0; }
  for (int im = 0; im < n_nodes * NDOF_FEM * max_elnd; im++) { Mat_q[im] = -1; }

  // generation matrix entries
   int kk_max=NDOF_FEM;int kj_max=NDOF_FEM;
   int kk_min=0;int kj_min=0;
   switch(name_label/10){
       case 0:
           kk_max=NDOF_FEM; kk_min=NDOF_FEM-1;
           break;
       case 1:
           kk_min=0;   kk_max=NDOF_P;
            break;
       case 2:
           kk_min=0;   kk_max=NDOF_FEM;
            break;
        case 3:
           kk_min=NDOF_FEM-1-_n_local_faces;   kk_max=NDOF_FEM-1;
            break;     
       default: abort();
   }
      switch(name_label%10){
       case 0:
           kj_max=NDOF_FEM; kj_min=NDOF_FEM-1;
           break;
       case 1:
           kj_min=0;   kj_max=NDOF_P;
            break;
       case 2:
           kj_min=0;   kj_max=NDOF_FEM;
            break;
        case 3:
           kj_min=NDOF_FEM-1-_n_local_faces;   kj_max=NDOF_FEM-1;
            break;     
       default: abort();
   }

   int n_faces = elem_conn[0][0];
 
  for (int pr = 0; pr < _n_subdomains; pr++) {
    for (int iel = off_el[pr * _n_levels + Level1]; iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
      int el_lib = v_el[iel].second;
      // std::cout<<"iel "<<el_lib<<std::endl;

      for (int kk = kk_min; kk <kk_max; kk++) {
        int op = v_inv_nd[elem_sto[el_lib][1 + kk]];
        int i_row = g_indexL[Level1][op];
       
        for (int kj = kj_min; kj <kj_max ; kj++)  {
          int opj = v_inv_nd[elem_sto[el_lib][1 + kj]];
          int ind_row = i_row * NDOF_FEM * max_elnd;
     
          Mat_q[ind_row + mem[i_row] + kj] = g_indexL[Level1][opj];
          // if (i_row == 531) std::cout<<"iel "<<iel<<" "<<g_indexL[Level1][opj]<<std::endl;
        }
        mem[i_row] += NDOF_FEM;

        for (int f = 0; f < n_faces; f ++){
          int el_side = elem_conn[el_lib][f + 1];
          if (el_side == -1 ) continue;

          for (int kj = kj_min; kj <kj_max ; kj++)  {
            int opj = v_inv_nd[elem_sto[el_side][1 + kj]];
            int ind_row = i_row * NDOF_FEM * max_elnd;
            Mat_q[ind_row + mem[i_row] + kj] = g_indexL[Level1][opj];
            // if (i_row == 531) std::cout<<"iel "<<iel<<" neigh "<<el_side<<" "<<g_indexL[Level1][opj]<<std::endl;
          }

          mem[i_row] += NDOF_FEM;
        }

        // const int stride = NDOF_FEM * max_elnd;
        // const int max_print = stride;

        // std::cout << "Row " << i_row << " (first " << max_print << " slots): ";
        // for (int j = 0; j < max_print; ++j) {
        //   int ind_row = i_row * NDOF_FEM * max_elnd;
        //   int v = Mat_q[ind_row + j];
        //   if (v < 0) std::cout << ".\t";
        //   else       std::cout << v << "\t";
        // }
        // std::cout << "\n";
          
        
      }

    }
    
  }  
  


  
  // compressing zeros and counting -------------------------------------
  // counters
  int count_q = 0;     len_q[0] = 0;
  int count_qoff = 0;  len_qoff[0] = 0;
  int ki = 0; int count_ki=0;
  // zeroing aux
  for (int i = 0; i < NDOF_FEM * max_elnd * n_nodes; i++) { aux[i] = 0; }
  // compressing
  for (int kdom = 0; kdom < _n_subdomains; kdom++) {
    for (int iki = 0; iki < (int)(max_nd[kdom * _n_levels + Level1] - min_nd[kdom * _n_levels + Level1]); iki++) {
      // setting aux
      for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
        int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
        if (knod_q >= 0) { aux[knod_q] = 1; }
      }
      // computing various flags and lengths
      for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
        int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
        if (knod_q >= 0 && aux[knod_q] == 1) {
          aux[knod_q] = 0;
          Mat_q[count_q] = knod_q;
          count_q++;
          if (knod_q < (int)offset_sub[kdom] || knod_q >= (int)offset_sub[kdom + 1]) { count_qoff++; }
        }
      }
 
      ki++;
      if(count_q-len_q[count_ki]>0) count_ki++;
      len_q[count_ki] = count_q;
      len_qoff[count_ki] = count_qoff;
      
    }
  }

  return count_q;
}



// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
