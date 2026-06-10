// std libraries ------------------------------------------
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// LibMesh includes ---------------------------------------
// #include "gencase_conf.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// // configure includes -------------------------------------
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

// #include <libmesh>
// ==============================================
/// This function print the mesh
// ----------------------------------------------
void MGGenCase::printMesh(
    BoundaryMesh& bd_msht,   // boundary mesh
    Mesh& msht,              // Volume mesh -> refined one
    BoundaryMesh& bd_msh0,   // boundary mesh
    Mesh& msh0,              // coarse mesh
    const int map_mesh_in[]  // mesh map
) {
  // ============================================

  // set up
  int n_nodes = msht.n_nodes();                     // from mesh
  int n_elements = msht.n_elem();                   // from mesh
  int n_elements_b = bd_msht.n_elem() * _n_levels;  // from mesh
  int n_nodes_b = bd_msht.n_nodes();                // from mesh
  std::cout << " \n bnodes   " << n_nodes_b << "\n";
  Mesh::const_element_iterator it_t00 = msh0.elements_begin();
  const Mesh::const_element_iterator end_t00 = msh0.elements_end();
  Mesh::const_element_iterator it_tr = msht.elements_begin();
  const Mesh::const_element_iterator end_tr = msht.elements_end();
  // // //

  // // // // active elemnt iterator
  //   Mesh::const_element_iterator        it_t00 = msh0.active_elements_begin();
  //   const Mesh::const_element_iterator  end_t00 = msh0.active_elements_end();
  //   Mesh::const_element_iterator           it_tr = msht.active_elements_begin();
  //   const Mesh::const_element_iterator   end_tr = msht.active_elements_end();

  const int n_faces = (*it_t00)->n_sides();
  _n_local_faces=n_faces ;
  //   const int  n_nodes_lev0=msh0.n_nodes();
  int n_variables = 1;
  // counting
  int* n_elements_lev = new int[2 * _n_levels];
  int* n_nodes_lev = new int[4 * _n_levels]; // 0: quad, 1: vertices, 2: face, 3: central 
  int* bd_n_elements_sl = new int[_n_subdomains * _n_levels + 1];
  int* n_nodes_sl_q = new int[_n_subdomains * _n_levels + 1];
  int* n_nodes_sl_v = new int[_n_subdomains * _n_levels + 1];
  int* n_elements_sl = new int[_n_subdomains * _n_levels + 1];
  int* n_nodes_sl_f = new int[_n_subdomains * _n_levels + 1];
  int* n_nodes_sl_c = new int[_n_subdomains * _n_levels + 1];
  
  int* elxnode = new int[n_nodes];

  // structures with original mesh -------
  double* nod_val;
  nod_val = new double[n_nodes * 3];
  int* nod_flag;
  nod_flag = new int[n_nodes];
  int* mat_flag;
  mat_flag = new int[n_elements];
  int* map;
  map = new int[n_elements * NDOF_FEM];
  int** elem_sto;
  elem_sto = new int*[n_elements];
  int** elem_conn;
  elem_conn = new int*[n_elements];
  int** bd_elem_sto;
  bd_elem_sto = new int*[n_elements_b];
  int** nod_sto;
  nod_sto = new int*[n_nodes];
  // order vectors -------------------
  std::vector<std::pair<int, int>> v_el(n_elements);
  int* v_inv_el = new int[n_elements];
  std::vector<std::pair<int, int>> v(n_nodes);
  int* v_inv_nd = new int[n_nodes];
  //  local bc and material for each element
  int bc_tmp[NDOF_FEM];
  int maptmp[NDOF_FEM];
  // zeroes ------------------
  for (int i = 0; i < n_nodes; i++) {
    nod_flag[i] = 0;  // bc condition
  }
  for (int i = 0; i < n_nodes * 3; i++) {
    nod_val[i] = 0.;  // nodes
  }
  for (int i = 0; i < n_elements * NDOF_FEM; i++) {
    map[i] = -1;  // nodes
  }
  for (int i = 0; i < _n_levels * _n_subdomains; i++) { n_elements_sl[i] = 0; }
  for (int i = 0; i < _n_levels * _n_subdomains; i++) { bd_n_elements_sl[i] = 0; }
  for (int ilev = 0; ilev < _n_levels * 2; ilev++) {
    n_elements_lev[ilev] = 0;  // nodes
  }
  for (int ilev = 0; ilev < 4 * _n_levels; ilev++) {
    n_nodes_lev[ilev] = 0;  // nodes
  }
  for (int is = 0; is <= _n_subdomains * _n_levels; is++) { n_nodes_sl_q[is] = 0;n_nodes_sl_v[is] = 0; n_nodes_sl_f[is] = 0; n_nodes_sl_c[is] = 0; }
  for (int i = 0; i < n_elements; i++) {
    // elements
    elem_sto[i] = new int[_dcl_nel];
    mat_flag[i] = 1;  // material flag (default->1)
    for (int k = 0; k < _dcl_nel; k++) { elem_sto[i][k] = -1; }
    elem_conn[i] = new int[n_faces + 1];
    elem_conn[i][0] = n_faces;
    for (int k = 1; k < n_faces + 1; k++) { elem_conn[i][k] = -1; }
  }
  for (int i = 0; i < n_elements_b; i++) {
    // bd elements
    bd_elem_sto[i] = new int[_dclb_nel];
    for (int j = 0; j < _dclb_nel; j++) { bd_elem_sto[i][j] = 0; }
  }
  //   ================================================
  //    ELEMENTS (in elem_sto and bd_elem_sto)
  // ---------------------------------------------
  //    elem_sto[Id(0), nodes (NDOF_FEM),
  //           lev(NDOF_FEM+1),pr(NDOF_FEM+2),parent (NDOF_FEM+3),
  //           Nchildren(NDOF_FEM+4), CHILDREN]
  // --------------------------------------------------
  //    bd_elem_sto[Id(0), elem_id (1),
  //           lev(2),side(3), n_nodes (4),
  //           side-nodes]
  // ================================================
  // Storing the mesh information
  int count_e = 0;
  int count_eb = 0;
  int n_nodes_el;

  //   Mesh::const_element_iterator       it_tr = msht->elements_begin();
  //   const Mesh::const_element_iterator end_tr= msht->elements_end();
  for (; it_tr != end_tr; ++it_tr) {
    Elem* elem = *it_tr;
    // element id
    int id_el = elem->id();
    elem_sto[count_e][0] = id_el;

    // element nodes
    n_nodes_el = elem->n_nodes();
    //     int a= NDOF_FEM;
    for (int inode = 0; inode < n_nodes_el; inode++) {
      int knode = elem->node_id(inode);
      elem_sto[count_e][1 + inode] = knode;
      // coordinates storage
      for (int idim = 0; idim < DIMENSION; idim++) {
        double xyz = msht.point(knode)(idim);
        nod_val[knode + idim * n_nodes] = xyz;
      }
//       printf("libmesh %d %d %g %g %g\n", inode, knode, nod_val[knode + 0 * n_nodes], nod_val[knode + 1 * n_nodes], nod_val[knode + 2 * n_nodes]);
    }
    // level and subdomain flag
    int lev = elem->level();
    elem_sto[count_e][NDOF_FEM + 1] = lev;
    int subdom_id = elem->processor_id();
    elem_sto[count_e][NDOF_FEM + 2] = subdom_id;
    // parent flag
    Elem* parent = elem->parent();
    if (parent != NULL) { elem_sto[count_e][NDOF_FEM + 3] = parent->id(); }
    // children
    int n_childs = elem->n_children();
    if (elem->has_children()) {
      elem_sto[count_e][NDOF_FEM + 4] = n_childs;
      // children in the element
      for (int i_ch = 0; i_ch < n_childs; i_ch++) {
        Elem* child = (*it_tr)->child_ptr(i_ch);
        elem_sto[count_e][NDOF_FEM + 5 + i_ch] = child->id();
      }
    }

    //  boundary ------------------------------------------
    for (int s = 0; s < (int)elem->n_sides(); s++) {
      if (elem->neighbor_ptr(s) == NULL) {
        bd_elem_sto[count_eb][0] = count_eb;
        bd_elem_sto[count_eb][NDOF_FEMB + 1] = elem->id();
        bd_elem_sto[count_eb][NDOF_FEMB + 2] = lev;
        n_elements_lev[lev + _n_levels]++;
        bd_elem_sto[count_eb][NDOF_FEMB + 3] = s;
        std::unique_ptr<Elem> side(elem->build_side_ptr(s));
        bd_elem_sto[count_eb][NDOF_FEMB + 4] = (int)side->n_nodes();
        for (int ns = 0; ns < (int)side->n_nodes(); ns++) {
          bd_elem_sto[count_eb][1 + ns] = side->node_id(ns);
        }
        count_eb++;
      } else {
        elem_conn[count_e][s + 1] = (elem->neighbor_ptr(s))->id();
      }
    }
    count_e++;
  }

  //    ================================================================

  int bd_n_elements = count_eb;
  assert(n_elements == count_e);
  int count_lev0 = 0;
  // ================================================
  // ELEMENT GROUPING  (lev(0)>lev(n))
  // ================================================
  // mesh0 sets subdomain (prodcessor id)
  //   Mesh::const_element_iterator       it_t00 = msh0->elements_begin();
  //   const Mesh::const_element_iterator end_t00= msh0->elements_end();

  for (; it_t00 != end_t00; ++it_t00) {
    // mesh00 processors
    Elem* elem = *it_t00;                                  // element pt
    int el_id = elem->id();                                // el identity
    elem_sto[el_id][NDOF_FEM + 2] = elem->processor_id();  // processor
    elem_sto[el_id][NDOF_FEM + 1] = 0;                     // father (coarse=0)
    count_lev0++;                                          // count
  }

  int n_groups_names = 1;
  int* group_id_names = new int[1];
#ifdef MATBC_INTERFACE

  const int libmesh_gen = (int)_mgutils._geometry["libmesh_gen"];  // gen param flag
  if (libmesh_gen == 0) {
    // ****************************************************************************
    // Boundary condition and material from file
    // file with mat bc info ------------------------------------------------------
    hid_t status = 0;   // read error flag
    hsize_t dims10[2];  // dimension vector for dset (two index)
    // file with bc and mat cond --------------------------------------------------
    std::ostringstream order_bc;  // file name to read
    order_bc << _mgutils._inout_dir << _mgutils.get_file("BASEBC");
    hid_t file_id = H5Fopen(order_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // Reading bc_flag ----------------------------------------------------------
    hid_t dtset = H5Dopen(
        file_id, "/bc_flag"
#if HDF5_VERSIONM != 1808
        ,
        H5P_DEFAULT
#endif
    );                                                            // get dset
    hid_t filespace = H5Dget_space(dtset);                        // get filespace
    status = H5Sget_simple_extent_dims(filespace, dims10, NULL);  // get dim
    int n_nodes_lev00 = dims10[0];                                // level zero original element -> 00
    int* bc_condition_lev0 = new int[dims10[0]];                  // coarse bc
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, bc_condition_lev0);  // read
    assert(status == 0);                                                                        // check

    // Reading mat_condition ----------------------------------------------------
    dtset = H5Dopen(
        file_id, "/material"
#if HDF5_VERSIONM != 1808
        ,
        H5P_DEFAULT
#endif
    );                                                            // get dset
    filespace = H5Dget_space(dtset);                              // get filespace
    status = H5Sget_simple_extent_dims(filespace, dims10, NULL);  // get dim
    int n_elements_lev0 = dims10[0];                              // number of elem c level
    int* mat_condition_lev0 = new int[n_elements_lev0];           // coarse mat c
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat_condition_lev0);  // read
    assert(status == 0);                                                                         // check

    //   int n_groups_names=1;
    //   int * group_id_names =new int[1];

#ifdef HAVE_GROUP
    // Reading group name  ----------------------------------------------------------
    dtset = H5Dopen(file_id, "/group_names", H5P_DEFAULT);                                                            // get dset
    filespace = H5Dget_space(dtset);                              // get filespace
    status = H5Sget_simple_extent_dims(filespace, dims10, NULL);  // get dim
    n_groups_names = dims10[0];                                   // level zero original element -> 00
    delete[] group_id_names;
    group_id_names = new int[dims10[0]];                                                     // coarse bc
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_id_names);  // read
    assert(status == 0);
#endif  // end HAVE_GROU  ----------------------------------------------------------
        // creating map file-mesh  to libmesh-mesh ----------------------------------
    for (int i = 0; i < n_nodes_lev00; i++) { nod_flag[i] = bc_condition_lev0[i]; }
    const int order_quad = (int)_mgutils._geometry["second_order"];  // param
    for (int ielem = 0; ielem < n_elements_lev0; ielem++) {
      mat_flag[ielem] = mat_condition_lev0[ielem];
      if (order_quad != 0) {
        for (int inode = 0; inode < NDOF_FEM; inode++) { maptmp[inode] = elem_sto[ielem][1 + inode]; }

        // writing bc conditions in the elements (min wins)
        for (int i = 0; i < NDOF_FEM; i++) {
          // loop over all child nodes
          int min = 100000000;
          int count = 0;
          for (int j = 0; j < NDOF_P; j++) {
            // loop over father nodes
            int bc_tmp = nod_flag[maptmp[j]];

            if (MGGeomEl::Prol[j + i * NDOF_P] != 0) {
              // nodes with proj
              if (min > bc_tmp) {
                min = bc_tmp;  // min wins
              }
              count++;
            }
          }
          if (count < 5) {
            nod_flag[maptmp[i]] = min;  // store the value
          }
        }
        // +++++++++++++++++++++++++++++

      }  // ibc_flag
    }
    //     // bc and material *****************************************

    // clean
    delete[] bc_condition_lev0;
    delete[] mat_condition_lev0;
  }
#endif

  // ****************************************************************************
  //  From coarse level to level n (ordering and filling)
  // subdomain ordering from level 0
  for (int ilev = 0; ilev < _n_levels - 1; ilev++) {
    // level loop
    for (int ielem = 0; ielem < n_elements; ielem++) {
      if (elem_sto[ielem][NDOF_FEM + 1] == ilev) {
        // only el with same level
        for (int ich = 0; ich < elem_sto[ielem][NDOF_FEM + 4]; ich++) {
          // child loop

          elem_sto[elem_sto[ielem][NDOF_FEM + 5 + ich]][NDOF_FEM + 2] =
              elem_sto[ielem][NDOF_FEM + 2];  // writing into my child the rank

#ifdef MATBC_INTERFACE
          // MATBC_INTERFACE ++++++++++++++++++++++++++++++++++++++++++++++++++++
          // child material flag --------------------
          mat_flag[elem_sto[ielem][NDOF_FEM + 5 + ich]] = mat_flag[ielem];

          // child node flag -----------------------------
          // storing  the  NDOF_FEM  nod_flag nodes of the child (ich)
          for (int inode = 0; inode < NDOF_FEM; inode++) {
            int icnode = elem_sto[elem_sto[ielem][NDOF_FEM + 5 + ich]][1 + inode];
            bc_tmp[inode] = nod_flag[icnode];  // child node flag
            maptmp[inode] = icnode;            // child nodes
          }

          // writing bc conditions in the elements (min wins)
          for (int i = 0; i < NDOF_FEM; i++) {
            // loop over all child nodes
            int min = 100000000;
            for (int j = 0; j < NDOF_P; j++) {
              // loop over father nodes
              if (MGGeomEl::Prol[j + i * NDOF_P] != 0) {
                // nodes with proj
                if (min > bc_tmp[j]) {
                  min = bc_tmp[j];  // min wins
                }
              }
            }
            nod_flag[maptmp[i]] = min;  // store the value in the child node
          }
#endif
          // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }  // end for ich
      }    //  ==ilev
    }      // end for ielem
  }        // end for (ilev)
  // *************************************************************

  // ================================================
  // REORDERING ELEMENTS (lev>subdomain)
  // ================================================

  // v_el[i].first  : level at which cell i belongs
  // v_el[i].second : cell id, with global numbering
  // level 0: nl0 cells, level 1: nl1 cells, level 2: nl2 cells
  // then global cells numbering is from 0 to (nl0 + nl1 + nl2)

  std::vector<int> ElementsPerLevel(_n_levels);
  for (int k = 0; k < _n_levels; k++) { ElementsPerLevel[k] = 0; }
  for (int kel = 0; kel < n_elements; kel++) {
    // counting (level and subdom)
    int proc_id = elem_sto[kel][NDOF_FEM + 2];
    int level_id = elem_sto[kel][NDOF_FEM + 1];
    n_elements_lev[level_id]++;
    n_elements_sl[level_id + proc_id * _n_levels]++;
    // setup ordering vector
    v_el[kel].first = level_id + _n_levels * proc_id;
    v_el[kel].second = kel;

    _ElPerProcPerLevel[proc_id][level_id]++;
    ElementsPerLevel[level_id]++;
  }

  for (int k = 0; k < _n_levels; k++) {
    std::cout << "Elements in level " << k << ": " << ElementsPerLevel[k] << std::endl;
    std::cout << "Parallel subdivision: ";
    for (int pp = 0; pp < _n_subdomains; pp++)
      std::cout << " proc: " << pp << " n elem " << _ElPerProcPerLevel[pp][k];
    std::cout << std::endl;
  }

  int count2 = 0;
  for (int pp = 0; pp < _n_subdomains; pp++) {
    for (int k = 0; k < _n_levels; k++) {
      std::cout << " p" << pp << "l" << k << ": from " << count2 << " <-> "
                << count2 + _ElPerProcPerLevel[pp][k] - 1 << " |";
      count2 += _ElPerProcPerLevel[pp][k];
    }
  }
  std::cout << "Total Elements " << n_elements << std::endl;

  // std reordering on pairs NDOF_FEM
  std::sort(v_el.begin(), v_el.end());
  for (int i = 0; i < n_elements; i++) {
    int kel = v_el[i].second;
    v_inv_el[kel] = i;
  }

  // boundary ---------------------------------
  std::vector<std::pair<int, int>> v_elb(bd_n_elements);
  int* v_inv_elb = new int[bd_n_elements];
  for (int kel = 0; kel < bd_n_elements; kel++) {
    int isubdom = elem_sto[bd_elem_sto[kel][NDOF_FEMB + 1]][NDOF_FEM + 2];
    int ilev = bd_elem_sto[kel][NDOF_FEMB + 2];

    bd_n_elements_sl[ilev + isubdom * _n_levels]++;
    v_elb[kel].first = ilev + _n_levels * isubdom;
    v_elb[kel].second = kel;
  }
  // std reordering on pairs
  std::sort(v_elb.begin(), v_elb.end());
  //  print
  for (int i = 0; i < bd_n_elements; i++) {
    int kel = v_elb[i].second;
    v_inv_elb[i] = kel;
  }

  // ================================================
  //    NODES (in nd_sto and nod_val)
  // ---------------------------------------------
  //    nd_sto[Index(0), pr(1),lev (2), vertex(3), edge(6), face(4), center(5)]
  // --------------------------------------------------
  //    nod_val[x,y,z]
  // ================================================

  // Set up -------------------------------------
  for (int i = 0; i < n_nodes; i++) {
    nod_sto[i] = new int[_N_NDV + 3];
    for (int k = 0; k < _N_NDV + 3; k++) { nod_sto[i][k] = -1; }
    nod_sto[i][2] = 1000000;
    nod_sto[i][3] = _n_levels;
      nod_sto[i][4] = _n_levels;
        nod_sto[i][5] = _n_levels;
        nod_sto[i][6] = _n_levels;
  }

  // Storing the node information--------------
  for (int i = 0; i < n_elements; i++) {
    for (int k = 0; k < n_nodes_el; k++) {
      int knode = elem_sto[i][k + 1];
      nod_sto[knode][0] = knode;
      if (elem_sto[i][NDOF_FEM + 2] > nod_sto[knode][1]) {
        nod_sto[knode][1] = elem_sto[i][NDOF_FEM + 2];  // subdomain
      }
      if (elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][2]) {
        nod_sto[knode][2] = elem_sto[i][NDOF_FEM + 1];  // level
      }
      if (k < NDOF_P && elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][3]) {
        nod_sto[knode][3] = elem_sto[i][NDOF_FEM + 1];
      }
        if (k> NDOF_P-1 && k< NDOF_FEM-1-n_faces && elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][6]) {
        nod_sto[knode][6] = elem_sto[i][NDOF_FEM + 1];
      }
        if (k< NDOF_FEM-1 && k> NDOF_FEM-2-n_faces && elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][4]) {
        nod_sto[knode][4] = elem_sto[i][NDOF_FEM + 1];
      }
         if (k== NDOF_FEM-1 && elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][5]) {
        nod_sto[knode][5] = elem_sto[i][NDOF_FEM + 1];
      }
      
      // nodes=2(vertex);3(face);4(center)
//       if(elem_sto[i][NDOF_FEM + 1] == 0) {
//       if(nod_sto[knode][4] == -1) nod_sto[knode][4] = 1;
//       if(k< NDOF_P) nod_sto[knode][4] = std::max(2, nod_sto[knode][4]);
//        if(k< NDOF_FEM-1 && k> NDOF_FEM-2-n_faces) nod_sto[knode][4] = std::max(3, nod_sto[knode][4]);
//        if(k== NDOF_FEM-1) nod_sto[knode][4] = std::max(4, nod_sto[knode][4]);
//       }
    }
  }

  // -----------------------------------------------
  // REORDERING NODES (subdom>lev>levP>var)
  for (int knode = 0; knode < n_nodes; knode++) {
    n_nodes_sl_q[_n_levels * nod_sto[knode][1] + nod_sto[knode][2]]++;
    if (nod_sto[knode][3] < _n_levels) { n_nodes_sl_v[_n_levels * nod_sto[knode][1] + nod_sto[knode][3]]++; }
      if (nod_sto[knode][4] < _n_levels) { n_nodes_sl_f[_n_levels * nod_sto[knode][1] + nod_sto[knode][4]]++; }
      if (nod_sto[knode][5] < _n_levels) { n_nodes_sl_c[_n_levels * nod_sto[knode][1] + nod_sto[knode][5]]++; }
    for (int iflag = 2; iflag< 6; iflag++)
        if (nod_sto[knode][iflag] < _n_levels)
           n_nodes_lev[(iflag-2) * _n_levels + nod_sto[knode][iflag]]++;
    
    // for (int klev = nod_sto[knode][2]; klev < _n_levels; klev++) { n_nodes_lev[klev]++; }
//   if (nod_sto[knode][2] == 0) 
//               n_nodes_lev[(nod_sto[knode][4] - 1) * _n_levels]++;
    
    
    // old !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      v[knode].first = nod_sto[knode][4] +
//                       n_variables * (nod_sto[knode][3] +
//                                      (_n_levels + 1) * (nod_sto[knode][2] + _n_levels * nod_sto[knode][1]));
   // new  v[knode].first=(proc)(level(2))(pressure(1))(node pos(1))
   // fun =102 1 1 -> p=1 l=2 pr=1 nodes=1(vertex)
      v[knode].first = nod_sto[knode][5] +
      10*(nod_sto[knode][4] +
      10* (nod_sto[knode][6] +
                      (10) * (nod_sto[knode][3] +
                                     (10) * (nod_sto[knode][2] +  _n_levels * nod_sto[knode][1]))));                
                     
    v[knode].second = knode;
  }

  for (int ilev=1; ilev<_n_levels; ilev++)
      for (int iflag = 2; iflag< 4; iflag++)
       n_nodes_lev[(iflag-2) * _n_levels + ilev] += n_nodes_lev[(iflag-2) * _n_levels + ilev - 1];
  
  sort(v.begin(), v.end());
  // order increasing functional fun  v[knode].first
  //  print
  for (int i = 0; i < (int)v.size(); i++) { v_inv_nd[v[i].second] = i; }

  for (int ind = 0; ind < n_nodes; ind++) { delete[] nod_sto[ind]; }
  delete[] nod_sto;

  // ============================================
  // OFFSET
  // ============================================
 
  int* off_el = new int[_n_subdomains * _n_levels + 1];   off_el[0] = 0;  int sum = 0;
  int* bd_off_el = new int[_n_subdomains * _n_levels + 1]; bd_off_el[0] = 0;  int sumb = 0;
//   int* off_nd[4];  
  int* min_nd[4]; int* max_nd[4];
//   for(int i_type=0;i_type<4;i_type++) {off_nd[i_type] = new int[_n_subdomains * _n_levels + 1]; }
  for(int i_type=0;i_type<4;i_type++) {min_nd[i_type] = new int[_n_subdomains * _n_levels]; max_nd[i_type] = new int[_n_subdomains * _n_levels];}
  
  int start_q = 0; int start_v = 0;
  for(int s=0;s<_n_subdomains; s++) {
    int offset_q = 0; int offset_v = 0; 
    for(int l=0; l<_n_levels; l++) {
        offset_q += n_nodes_sl_q[s*_n_levels + l];
        offset_v += n_nodes_sl_v[s*_n_levels + l];
        min_nd[0][s * _n_levels + l] = start_q;
        max_nd[0][s * _n_levels + l] = start_q + offset_q;
        min_nd[1][s * _n_levels + l] = start_v;
        max_nd[1][s * _n_levels + l] = start_v + offset_v;
        min_nd[2][s * _n_levels + l] = start_q + offset_q - n_nodes_sl_c[s * _n_levels + l] - n_nodes_sl_f[s * _n_levels + l];
        max_nd[2][s * _n_levels + l] = start_q + offset_q - n_nodes_sl_c[s * _n_levels + l];
        min_nd[3][s * _n_levels + l] = start_q + offset_q - n_nodes_sl_c[s * _n_levels + l];
        max_nd[3][s * _n_levels + l] = start_q + offset_q;
    }
    start_q = max_nd[0][s * _n_levels + _n_levels - 1];
    start_v = max_nd[0][s * _n_levels + _n_levels - 1];
  }
  
// off_nd[0][0] = 0;off_nd[1][0] = 0;
//   off_nd[2][0] = n_nodes_sl_q[0] - n_nodes_sl_c[0] - n_nodes_sl_f[0];
//   off_nd[3][0] = n_nodes_sl_q[0] - n_nodes_sl_c[0];
//   int sum_nd_q = 0; int sum_nd_v = 0;
// 
// 
for (int isub = 0; isub < _n_subdomains * _n_levels; isub++) {
     sum += n_elements_sl[isub];       sumb += bd_n_elements_sl[isub];
//     sum_nd_q += n_nodes_sl_q[isub];    sum_nd_v += n_nodes_sl_v[isub];
     off_el[isub + 1] = sum;                      bd_off_el[isub + 1] = sumb;
//     off_nd[0][isub + 1] = sum_nd_q;              off_nd[1][isub + 1] = sum_nd_v;
//     off_nd[2][isub + 1] = off_nd[2][isub] + n_nodes_sl_f[isub];    off_nd[3][isub + 1] = off_nd[3][isub] + n_nodes_sl_c[isub];
    
}
  // max number of element on each node
  for (int inode = 0; inode < n_nodes; inode++) { elxnode[inode] = 0; }
  for (int i = 0; i < n_elements; i++)  // ALERT
    for (int k = 0; k < n_nodes_el; k++) {
      if (elem_sto[i][1 + NDOF_FEM] == _n_levels - 1) { elxnode[elem_sto[i][1 + k]]++; }
    }

  int maxelxnode = 0;
  for (int inode = 0; inode < n_nodes; inode++)
    if (elxnode[inode] > maxelxnode) { maxelxnode = elxnode[inode]; }
  std::cout << " Max number of element on each node = " << maxelxnode << " \n";

  // clean
  delete[] elxnode;
  delete[] bd_n_elements_sl; delete[] n_elements_sl;
  delete[] n_nodes_sl_v; delete[] n_nodes_sl_f; delete[] n_nodes_sl_c;  delete[] n_nodes_sl_q;
  
 


  // =================
  // node map -------------------------------------------
  // packaging
  int** g_indexL = new int*[4*_n_levels]; int  i_count[4];
  
  for (int ilev = 0; ilev < _n_levels; ilev++) {
    for(int i_type=0;i_type<4;i_type++){
        g_indexL[ilev+i_type*_n_levels] = new int[n_nodes];
        i_count[i_type]=0;
        for (int inode = 0; inode < n_nodes; inode++)  g_indexL[ilev+i_type*_n_levels][inode] = -1;    
        for (int iproc = 0; iproc < _n_subdomains; iproc++) {
//           for (int inode = off_nd[i_type][iproc * _n_levels]; inode < off_nd[i_type][iproc * _n_levels + ilev + 1];inode++) 
//                                                                                  g_indexL[ilev+i_type*_n_levels][inode] = i_count[i_type]++;
            int base = (i_type < 2)? 0 : ilev;
          for (int inode = min_nd[i_type][iproc * _n_levels + base]; inode < max_nd[i_type][iproc * _n_levels + ilev];inode++) 
                                                                                 g_indexL[ilev+i_type*_n_levels][inode] = i_count[i_type]++;
       }
    }
  }
  
//   int count = 0;
//   for (int iproc = 0; iproc < _n_subdomains; iproc++) {
//     for (int inode = 0; inode < off_nd[1][iproc * _n_levels + 1] - off_nd[1][iproc * _n_levels]; inode++) {
//       g_indexL[_n_levels][off_nd[0][iproc * _n_levels] + inode] = count++;
//     }
//   }

  //   print mesh
  print_mesh_h5(
      n_nodes_lev, map_mesh_in, n_nodes, nod_val,  min_nd, max_nd, n_elements, n_elements_lev, bd_n_elements,
      n_nodes_b, off_el, v_inv_nd, v_inv_el, elem_sto, elem_conn, v_el, bd_elem_sto, v_elb, bd_off_el,
      g_indexL, v, nod_flag, mat_flag);

  for (int iel = 0; iel < n_elements_b; iel++) { delete[] bd_elem_sto[iel]; }
  delete[] bd_elem_sto;
  delete[] nod_val;
  delete[] bd_off_el;

  // data_in/mesh.xmf
  print_multimesh(n_elements_lev, n_nodes_lev);
  // print med

#ifdef HAVE_MED
  const int ibcflag = (int)_mgutils._geometry["ibc_gen"];  /// parameter defined in parameters.in

  if (ibcflag != 0) {
    // name file med-mesh
    if (libmesh_gen == 0) {
      std::string mesh_name = _mgutils.get_file("F_MESH_READ");
      unsigned pos = mesh_name.find(".");           // position of "live" in str
      std::string str3 = mesh_name.substr(0, pos);  // get from "live" to the end
      print_med(
          _n_levels - 1, bd_msht, msht, n_groups_names, group_id_names, v, v_elb, nod_flag, mat_flag,
          str3.c_str());
      print_MedToMg(
          _n_levels - 1, bd_msht, msht, n_groups_names, group_id_names, v, v_inv_nd, v_elb, nod_flag,
          mat_flag, off_el, v_inv_el, v_el, ElementsPerLevel, elem_sto, str3.c_str());
    } else {
      print_med(_n_levels - 1, bd_msht, msht, n_groups_names, group_id_names, v, v_elb, nod_flag, mat_flag);
      print_MedToMg(
          _n_levels - 1, bd_msht, msht, n_groups_names, group_id_names, v, v_inv_nd, v_elb, nod_flag,
          mat_flag, off_el, v_inv_el, v_el, ElementsPerLevel, elem_sto);
    }
  }
#endif

  delete[] n_elements_lev;  //

  // *********************************
  //          MGOperators
  // *********************************
  // #if EL_TYPE==18
  // return;
  // #else
  const int mgops_gen = (int)_mgutils._geometry["mgops_gen"];

  if (mgops_gen) {
    compute_and_print_MGOps(
        n_nodes_lev, maxelxnode, off_el, v_el, elem_sto, v_inv_nd, v_inv_el, g_indexL, min_nd,max_nd,elem_conn);
  }
  // clean
  delete[] v_inv_elb;
  delete[] v_inv_el;
  delete[] v_inv_nd;
  delete[] off_el;

  for (int iel = 0; iel < n_elements; iel++) { delete[] elem_sto[iel]; }
  delete[] elem_sto;
  for (int ind = 0; ind < 4*_n_levels; ind++) { delete[] g_indexL[ind]; } 
  delete[] g_indexL;
  for (int ind = 0; ind < 4; ind++)   {/*delete[] off_nd[ind];*/ delete[] min_nd[ind]; delete[] max_nd[ind];}
  delete[] n_nodes_lev;
// #endif
#ifdef PRINT_INFO
  std::cout << " MGGenCase::printMesh: Operators  printed \n";
#endif
  return;
}  // end print_mesh

// ====================================================
/// This function computes and prints MGOps
void MGGenCase::compute_and_print_MGOps(
    int* n_nodes_lev,                       // nodes for level
    int maxelxnode,                         // max # of elements for node
    int* off_el,                            // offset element
    std::vector<std::pair<int, int>> v_el,  // element ordering
    int** elem_sto,                         // element storage
    int* v_inv_nd,                          // node ordering
    int* v_inv_el, int** g_indexL,          // node maps
    int* min_nd[],                           // node subdomain  offset
    int *max_nd[],
    int** elem_conn
) {
  int max_elnd = maxelxnode;
#ifdef PRINT_INFO
  std::cout << " MGGenCase::compute_and_print_MGOps:  start \n";
#endif

  std::cout << " MGGenCase::compute_and_print_MGOps:  start \n";
  compute_matrix(n_nodes_lev, max_elnd, off_el, v_el, elem_sto, v_inv_nd, g_indexL, min_nd,max_nd,elem_conn);
#ifdef PRINT_INFO
  std::cout << " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
#endif
  compute_prol(n_nodes_lev, off_el, v_el, elem_sto, v_inv_nd, v_inv_el, g_indexL, min_nd,max_nd);
#ifdef PRINT_INFO
  std::cout << " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
#endif
  compute_rest(n_nodes_lev, max_elnd, off_el, v_el, elem_sto, v_inv_nd, v_inv_el, g_indexL, min_nd,max_nd);
#ifdef PRINT_INFO
  std::cout << " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
#endif

  // //===========================================
  // #ifdef Q2Q0                  //quadratic-piecewise
  //     compute_matrix_Q2Q0(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
  // #endif
  //     compute_prol_Q2Q0(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd, v_inv_el,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
  // #endif
  //     compute_rest_Q2Q0(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,v_inv_el,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
  // #endif
  //
  // #endif
  //
  // //==========================================
  // #ifdef Q2P1                 //quadratic-3piecewise
  //     compute_matrix_Q2P1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
  // #endif
  //     compute_prol_Q2P1(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd, v_inv_el,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
  // #endif
  //     compute_rest_Q2P1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,v_inv_el,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
  // #endif
  // #endif
  // //=======================================
  //
  // #ifdef Q2Q1                 //quadratic-linear
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps:  start \n";
  // #endif
  //     compute_matrix_Q2Q1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
  // #endif
  //     compute_prol_Q2Q1(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
  // #endif
  //     compute_rest_Q2Q1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,g_indexL,off_nd);
  // #ifdef PRINT_INFO
  //     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
  // #endif
  // #endif
  // //==========================

  // end switch

  return;
}

// ==============================================================
/// This function prints the mesh structure elem_sto of the MGGenCase class
void MGGenCase::print_lib_mesh(MGGeomEl& /*Eltype */, int** elem_sto, int n_elements, int dcl_nel) {
  //  const int NDOF_FEM = Eltype.n_q[0]; const int max_elnd = Eltype.n_l[0];
  std::cerr << " N : Id ,  NDOF_P , NDOF_FEM; (lev,pr,parent) ; "
            << " Nchildren: CHILDREN \n ";
  for (int i = 0; i < n_elements; i++) {
    std::cerr << i << " : " << elem_sto[i][0] << ", ";
    for (int k = 1; k < 1 + NDOF_P; k++) {
      std::cerr << elem_sto[i][k] << " ";  // nodes
    }
    std::cerr << ", ";
    for (int k = NDOF_P + 1; k < 1 + NDOF_FEM; k++) { std::cerr << elem_sto[i][k] << " "; }
    std::cerr << "; ";
    std::cerr << "(" << elem_sto[i][1 + NDOF_FEM] << ",";     //  level
    std::cerr << elem_sto[i][NDOF_FEM + 2] << ",";            //  lproc
    std::cerr << elem_sto[i][NDOF_FEM + 3] << ");  ";         //  parent
    std::cerr << " " << elem_sto[i][NDOF_FEM + 4] << " :  ";  //  n children
    for (int k = 1 + NDOF_FEM + 4; k < dcl_nel; k++) { std::cerr << elem_sto[i][k] << " "; }
    std::cerr << " \n ";
  }
  return;
}

// ================================================================
/// This function prints the mesh structure nod_sto of the MGGenCase class
void MGGenCase::print_lib_node(int** nod_sto, int n_nodes, int /*n_ndv */) {
  std::cerr << " N : Id ,  pr,lev (P); var \n ";
  for (int i = 0; i < n_nodes; i++) {
    std::cerr << i << " : " << nod_sto[i][0] << "; ";  // id
    std::cerr << nod_sto[i][1] << ",";                 //  lproc
    std::cerr << nod_sto[i][2] << " (";                //  level
    std::cerr << nod_sto[i][3] << ");  ";              //  level P
    std::cerr << nod_sto[i][4] << " \n ";              //  variable
  }
  return;
}

void MGGenCase::print_Mat(int* Mat, int nrow, int ncln) {
  std::cerr << " Matrix \n ";
  for (int i = 0; i < nrow; i++) {
    std::cerr << i << " - ";
    for (int j = 0; j < ncln; j++) { std::cerr << Mat[i * ncln + j] << " "; }
    std::cerr << "\n ";
  }
  return;
}

// ======================================================
/// This function print the mesh xdmf format
void MGGenCase::print_multimesh(
    const int* n_elements_lev,  // elements for level <-
    const int* n_nodes_lev      // nodes for level <-
) {
  // ==================================================

  std::string multimesh = _mgutils.get_file("MULTIMESH");
  std::string basemesh = _mgutils.get_file("BASEMESH");
  std::ostringstream inmesh_xmf;
  inmesh_xmf << _mgutils._inout_dir << multimesh << ".xmf";
  std::ofstream out(inmesh_xmf.str().c_str());

  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
  //   out << "[ <!ENTITY HeavyData \"mesh.h5\"> ] ";
  out << "> \n"
      << " \n";
  out << "<Xdmf> \n";
  out << "<Domain> \n";
  // Volume mesh ----------------------------------------------
  for (int ilev = 0; ilev < _n_levels; ilev++) {
    out << "<Grid Name=\"Mesh" << ilev << "\"> \n";
    out << "<Topology Type=\"" << _geomel.name[0] << "\" Dimensions=\"" << n_elements_lev[ilev] << "\"> \n";
    out << "<DataStructure DataType=\"Int\" Dimensions=\"" << n_elements_lev[ilev] << " " << NDOF_FEM
        << "\" Format=\"HDF\">  \n";
    out << basemesh << ".h5"  //"mesh.h5"
        << ":/ELEMS" << _geomel.n_q[0] << "/CONN/CONN_L" << ilev << "\n";
    out << "</DataStructure> \n";
    out << "</Topology> \n";
    out << "<Geometry Type=\"X_Y_Z\"> \n";
    for (int idim = 0; idim < 3; idim++) {
      out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_nodes_lev[_n_levels - 1]
          << " 1 \" Format=\"HDF\">  \n";
      out << basemesh << ".h5"  //"mesh.h5"
          << ":/NODES/COORD/X" << idim + 1 << " \n";
      out << "</DataStructure> \n";
    }
    out << " </Geometry>\n";
    // #ifdef MATBC_INTERFACE
    out << " <Attribute Name=\"Material\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    out << "<DataItem DataType=\"Int\" Precision=\"8\" Dimensions=\"" << n_elements_lev[ilev]
        << " 1 \" Format=\"HDF\">  \n";
    out << basemesh << ".h5"
        << ":/ELEMS" << _geomel.n_q[0] << "/MAT/MAT_L" << ilev << "\n";
    out << "</DataItem> \n";
    out << " </Attribute> \n";
    out << " <Attribute Name=\"PROCID\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    out << "<DataItem DataType=\"Int\" Precision=\"8\" Dimensions=\"" << n_elements_lev[ilev]
        << " 1 \" Format=\"HDF\">  \n";
    out << basemesh << ".h5"
        << ":/ELEMS" << _geomel.n_q[0] << "/PROCID_L" << ilev << "\n";
    out << "</DataItem> \n";
    out << " </Attribute> \n";
    out << " <Attribute Name=\"BC\" > \n";
    out << "<DataStructure DataType=\"Int\" Precision=\"8\" Dimensions=\"" << n_nodes_lev[_n_levels - 1]
        << " 1 \" Format=\"HDF\">  \n";
    out << basemesh << ".h5"
        << ":/NODES/BC \n";
    out << "</DataStructure> \n";
    out << " </Attribute> \n";
    // end MATBC_INTERFACE
    // #endif
    out << "</Grid> \n";
  }  // ---------- end volume mesh -----------------------------
  // Boundary mesh ----------------------------------------------
//   for (int ilev = 0; ilev < _n_levels; ilev++) {
//     out << "<Grid Name=\"Meshb" << ilev << "\"> \n";
//     out << "<Topology Type=\"" << _geomel.name[1] << "\" Dimensions=\"" << n_elements_lev[_n_levels + ilev]
//         << "\"> \n";
//     out << "<DataStructure DataType=\"Int\" Dimensions=\"" << n_elements_lev[_n_levels + ilev] << " "
//         << NDOF_FEMB << "\" Format=\"HDF\">  \n";
//     out << basemesh << ".h5"  //"mesh.h5"
//         << ":/ELEMS" << _geomel.n_q[0] << "/FEM1/MSH" << ilev << "\n";
//     out << "</DataStructure> \n";
//     out << "</Topology> \n";
//     out << "<Geometry Type=\"X_Y_Z\"> \n";
//     for (int idim = 0; idim < 3; idim++) {
//       out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_nodes_lev[_n_levels - 1]
//           << " 1 \" Format=\"HDF\">  \n";
//       out << basemesh << ".h5"  //"mesh.h5"
//           << ":/NODES/COORD/X" << idim + 1 << " \n";
//       out << "</DataStructure> \n";
//     }
//     out << " </Geometry>\n";
//     out << "</Grid> \n";
//   }
  // ---------- end boundary mesh -----------------------------

  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();

  return;
}

// =========================================================
/// This function print the mesh in hdf5 format (hdf5 version  1.8.8)
void MGGenCase::print_mesh_h5(
    int* n_nodes_lev, const int* map_mesh_in, int n_nodes, double* nod_val, int* min_nd[], int* max_nd[], int n_elements,
    int* n_elements_lev, int bd_n_elements, int bd_n_nodes, int* off_el, int* v_inv_nd, int* v_inv_el,
    int** elem_sto, int** elem_conn, std::vector<std::pair<int, int>> v_el, int** bd_elem_sto,
    std::vector<std::pair<int, int>> v_elb, int* bd_off_el, int** g_indexL,
    std::vector<std::pair<int, int>> v, int* nod_flag, int* mat_flag) {
  // ===================================================

  hid_t status = 0;
  // file h5 to print ------------------------------------
  std::ostringstream name;
  std::string basemesh = _mgutils.get_file("BASEMESH");
  std::string appl_dir = _mgutils.get_file("APPL_DIR");

  std::ostringstream inmesh;

  inmesh << _mgutils._inout_dir << basemesh << ".h5";
  hid_t file = H5Fcreate(inmesh.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  int n_meshes = 2 * _n_levels;

  //==================================
  // ROOT
  // =====================
  // _dim,_NoFamFEM,_NoLevels,Subdom ------------------
  // packaging data
  int* tdata;
  tdata = new int[4];
  tdata[0] = DIMENSION;  tdata[1] = 1;
  tdata[2] = _n_levels;  tdata[3] = _n_subdomains;
  // DFLS vector (hdf5 sorage)
  hsize_t dimsf[2];  dimsf[0] = 4;  dimsf[1] = 1;
  status = _mgutils.print_Ihdf5(file, "DFLS", dimsf, tdata);
  // clean
  assert(status == 0);
  delete[] tdata;
  //  _type_FEM -------------------------------
  //  data packaging
  tdata = new int[1]; tdata[0] = _geomel.n_q[0];
  // hdf5 sorage
  dimsf[0] = 1;  dimsf[1] = 1;
  status = _mgutils.print_Ihdf5(file, "GEOM_FAM", dimsf, tdata);
  // clean
  assert(status == 0);
  delete[] tdata;
  // // ===========================================
  // //  NODE COORDINATES  (COORD)
  // // ===========================================
  // ++++++++++++++++++++++++++++++++++++++++++++++++++
  hid_t group_id = H5Gcreate(file, "/NODES", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // packaging data
  int* n_point_b;
  n_point_b = new int[1]; n_point_b[0] = n_nodes;
  // hdf5 print
  dimsf[0] = 1; dimsf[1] = 1;
  name.str("");
  name << "/NODES/SIZE";
  status = _mgutils.print_Ihdf5(file, name.str(), dimsf, n_point_b);
  delete[] n_point_b;
  //      H5Gclose ( subgroup99_id );

  //  boundary conditions -----------------------------------------------------
  int* bc_con2;
  bc_con2 = new int[n_nodes];
  for (int inode = 0; inode < n_nodes; inode++) { bc_con2[inode] = nod_flag[v[inode].second]; }

  dimsf[0] = n_nodes;
  dimsf[1] = 1;
  name.str("");
  name << "/NODES/BC";
  status = _mgutils.print_Ihdf5(file, name.str(), dimsf, bc_con2);
  assert(status == 0);
  delete[] bc_con2;

  // node coordinates
  hid_t subgroup_id = H5Gcreate(file, "/NODES/COORD", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // packaging data
  // _xyz ---------------------------------------------------------------------
  // packaging data
  double* xcoord;
  xcoord = new double[n_nodes];
  for (int kc = 0; kc < 3; kc++) {
    for (int inode = 0; inode < n_nodes; inode++) { xcoord[inode] = nod_val[v[inode].second + kc * n_nodes]; }
    // hdf5 print
    dimsf[0] = n_nodes;
    dimsf[1] = 1;
    name.str("");
    name << "/NODES/COORD/X" << kc + 1;
    status = _mgutils.print_Dhdf5(file, name.str(), dimsf, xcoord);
  }
  assert(status == 0);


  // clean --------------------------------------------------------------------
  delete[] xcoord;
  H5Gclose(subgroup_id);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++
  //  /NODES/MAP
  subgroup_id = H5Gcreate(file, "/NODES/MAP", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // +++++++++++++++++++++++++++++++++++++++++++++++++
   std::cout << "\n ------ OFFSETV_LS ----- \n ";
   for (int i = 0; i < _n_subdomains; i++) {
    std::cout << " Proc " << i << "\n ";
    for (int j = 0; j < _n_levels; j++) {
      std::cout << " Level Nodes " << j << " from " << min_nd[0][j + i * _n_levels] << " to "
                << max_nd[0][j + i * _n_levels] << "\n ";
      std::cout << " Level Cells " << j << " from " << off_el[j + i * _n_levels] << " to "
                << off_el[j + 1 + i * _n_levels] << "\n ";
    }
  }
  // number of nodes on levels (Q(node) +V(vertex)+F(face)+C(center))
   dimsf[0] = 4 * _n_levels;  dimsf[1] = 1;
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/SIZE_QVFC", dimsf, n_nodes_lev);
  // node (quad) offset
  dimsf[0] = _n_subdomains * _n_levels;  dimsf[1] = 1;
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MINQ_SL", dimsf, min_nd[0]);
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MINV_SL", dimsf, min_nd[1]);
   status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MINF_SL", dimsf, min_nd[2]);
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MINC_SL", dimsf, min_nd[3]);
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MAXQ_SL", dimsf, max_nd[0]);
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MAXV_SL", dimsf, max_nd[1]);
   status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MAXF_SL", dimsf, max_nd[2]);
  status = _mgutils.print_Ihdf5(file, "/NODES/MAP/MAXC_SL", dimsf, max_nd[3]);
 std::cout << " \n ";
  // vertex (linear) offset 
 
  // node map -------------------------------------------------------------------
  // level map (hdf5 print)
  dimsf[0] = n_nodes;  dimsf[1] = 1; // dim dataset
  for (int ilev = 0; ilev < _n_levels; ilev++) {
    name.str("");  name << "/NODES/MAP/MAPQ_L" << ilev; // name dataset
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev]);
    name.str("");  name << "/NODES/MAP/MAPV_L" << ilev; // name dataset
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev + _n_levels]);
    name.str("");  name << "/NODES/MAP/MAPF_L" << ilev; // name dataset
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev + 2*_n_levels]);
    name.str("");  name << "/NODES/MAP/MAPC_L" << ilev; // name dataset
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev + 3*_n_levels]);
  }
  // vertex (linear) coarse level map (hdf5 print)
//   name.str("");  name << "/NODES/MAP/MAPV_L" << 0;
//   status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[_n_levels]);
//    for (int ilev = 0; ilev < _n_levels-1; ilev++) {
//     name.str("");  name << "/NODES/MAP/MAPV_L" << ilev+1; // name dataset
//     status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev]);
//   }
  
  // clean group -------------------------------------------
  H5Gclose(subgroup_id);
  H5Gclose(group_id);
  
  
 /* 
  
  
  
  dimsf[0] =  n_nodes;
  dimsf[1] = 1;
   for (int ilev = 0; ilev < _n_levels; ilev++) {
    name.str("");
    // status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL);
    hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
    hid_t dataset = H5Dcreate(file, name.str().c_str(), H5T_NATIVE_INT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, g_indexL);
    H5Sclose(dataspace);
    H5Dclose(dataset);*/
//   // }
//   // linear coarse level map (hdf5 print)
//   name.str("");
//   name << "/NODES/MAP/MAPV_L" << _n_levels;
//   status = _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[_n_levels]);
//   H5Gclose(subgroup_id);
//   H5Gclose(group_id);

  // ===========================================
  //   /ELEMS (CONNECTIVITY)
  // ===========================================
  std::string group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]);
  group_id = H5Gcreate(file, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // +++++++++++++++++++++++++++++++++++
  //  /ELEMS27/FEM1  (volume mesh)
  // +++++++++++++++++++++++++++++++++++
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/CONN"; 
  subgroup_id = H5Gcreate(file, group_name.c_str(), H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/MAT";
  hid_t subgroup_id2 = H5Gcreate(file, group_name.c_str(), H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  dimsf[0] = 4;  dimsf[1] = 1;
  int ndofm[4];  ndofm[0] =  _geomel.n_q[0];  ndofm[1] = _geomel.n_l[0];  ndofm[2] =  _geomel._n_sides[0]; ndofm[3] = 1;
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/NDOF_QVFC";
  status = _mgutils.print_Ihdf5(file, group_name.c_str(), dimsf, ndofm);

  // NoElements ------------------------------------
  dimsf[0] = _n_levels;  dimsf[1] = 1;
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/SIZE_L";
  status = _mgutils.print_Ihdf5(file, group_name.c_str(), dimsf, &n_elements_lev[0]);

  // offset
  dimsf[0] = _n_subdomains * _n_levels + 1;
  dimsf[1] = 1;
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/OFFSET_SL";
  status = _mgutils.print_Ihdf5(file, group_name.c_str(), dimsf, off_el);

  // Connectivity element-node all levels  packaging data  (volume) ------
  int* tempconn;
  tempconn = new int[n_elements * NDOF_FEM];  // coonectivity all mesh
  int* temp_material = new int[n_elements];
  for (int ielem = 0; ielem < n_elements; ielem++) {
    for (int inode = 0; inode < NDOF_FEM; inode++) {
      tempconn[inode + ielem * NDOF_FEM] =
          v_inv_nd[elem_sto[v_el[ielem].second][inode + 1]];  // swap for proc and lev
         // printf(" conn %d  %d %d %d \n  ", inode,v_inv_nd[elem_sto[v_el[ielem].second][inode + 1]], v_el[ielem].second, elem_sto[v_el[ielem].second][inode + 1]    );
    }
    temp_material[ielem] = mat_flag[v_el[ielem].second];
  }
  dimsf[0] = n_elements * NDOF_FEM;
  dimsf[1] = 1;  // global mesh hdf5 storage
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/CONN/CONN_ALL";
  status = _mgutils.print_Ihdf5(file, group_name.c_str(), dimsf, tempconn);

  // material condition on lev (ilev)    ------------------------------------
  dimsf[0] = n_elements;
  dimsf[1] = 1;
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/MAT/MAT";
  status = _mgutils.print_Ihdf5(file, group_name, dimsf, temp_material);

  int* temp_procid = new int[n_elements];
  for (int ielem = 0; ielem < n_elements; ielem++) {
    temp_procid[v_inv_el[ielem]] = elem_sto[ielem][NDOF_FEM + 2];
  }
  dimsf[0] = n_elements;
  dimsf[1] = 1;
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/PROCID";
  status = _mgutils.print_Ihdf5(file, group_name, dimsf, temp_procid);


  // -----------------------------------------------------------------------------

  // Connectivity element-element (neighbrs) at all levels (volume) -------
  const int n_faces = elem_conn[0][0];
  int* tempconn_el;
  tempconn_el = new int[n_elements * n_faces];
  for (int ielem = 0; ielem < n_elements; ielem++) {
    for (int iside = 0; iside < n_faces; iside++) {
      int el_side = elem_conn[v_el[ielem].second][iside + 1];
      // internal element ->  v_inv_el[el_side]; boundary element -1
        // std::cout << "ielem new: " << ielem << "        ielem old: "
        //     <<  v_el[ielem].second << "        el_side: " << el_side <<  "       v_inv_el:  "
        //     <<   v_inv_el[el_side] << std::endl;
      if (el_side != -1) {
        tempconn_el[iside + ielem * n_faces] = v_inv_el[el_side];
      } else {
        tempconn_el[iside + ielem * n_faces] = -1;  // boundary element -1
      }
    }
  }

  dimsf[0] = n_elements * n_faces;
  dimsf[1] = 1;  // global mesh hdf5 storage
  group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/EL_NEIG";
  status = _mgutils.print_Ihdf5(file, group_name, dimsf, tempconn_el);
  int icount_m1=0;
   for (int iface = 0; iface < n_elements * n_faces; iface++) {
        if(tempconn_el[iface]==-1) icount_m1++;
   }
   _n_glob_faces=(n_elements * n_faces-icount_m1)/2+icount_m1;
  
  
  // --------------------------------------------------------------------------
  // level connectivity (for each level) (start from zero)
#if ELTYPE == 18
  int wedge_18[18] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17};
#endif

  for (int ilev = 0; ilev < _n_levels; ilev++) {
    // mesh level vector for connectivity and material
    int* tempconnf = new int[n_elements_lev[ilev] * NDOF_FEM];
    int* temp_material2 = new int[n_elements_lev[ilev]];
    int ltot = 0;
    for (int iproc = 0; iproc < _n_subdomains; iproc++) {
      for (int iel = off_el[iproc * _n_levels + ilev]; iel < off_el[iproc * _n_levels + ilev + 1]; iel++) {
        temp_material2[ltot] = temp_material[iel];
        for (int inode = 0; inode < NDOF_FEM; inode++) {
          tempconnf[ltot * NDOF_FEM + inode] = tempconn
              [iel * NDOF_FEM +
#if ELTYPE == 18
               wedge_18[inode]
#else
               inode
#endif
          ];
        }  // end for inode
        ltot++;
      }  // end for iel
    }    // end for iproc
    // mesh connectivity on lev (ilev)  ---------------------------------------
    dimsf[0] = n_elements_lev[ilev] * NDOF_FEM;
    dimsf[1] = 1;
    group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/CONN/CONN_L" + std::to_string(ilev);
    status = _mgutils.print_Ihdf5(file, group_name, dimsf, tempconnf);

    // material condition on lev (ilev)    ------------------------------------
    dimsf[0] = n_elements_lev[ilev];
    dimsf[1] = 1;
    group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/MAT/MAT_L" + std::to_string(ilev);
    status = _mgutils.print_Ihdf5(file, group_name, dimsf, temp_material2);

    int* temp_procid2 = new int[n_elements_lev[ilev]];
    int ltot_procid = 0;
    for (int iproc = 0; iproc < _n_subdomains; iproc++) {
        std::cout <<  off_el[iproc * _n_levels + ilev] << "  " <<  off_el[iproc * _n_levels + ilev + 1] <<  std::endl;
      for (int iel = off_el[iproc * _n_levels + ilev]; iel < off_el[iproc * _n_levels + ilev + 1]; iel++) {
        temp_procid2[ltot_procid] = iproc;
        ltot_procid++;
      }  // end for iel
    }    // end for iproc
    // material condition on lev (ilev)    ------------------------------------
    dimsf[0] = n_elements_lev[ilev];
    dimsf[1] = 1;
    group_name = "/ELEMS" + std::to_string( _geomel.n_q[0]) + "/PROCID_L" + std::to_string(ilev);
    status = _mgutils.print_Ihdf5(file, group_name, dimsf, temp_procid2);

    // clean lev structure
    delete[] temp_material2;
    delete[] temp_procid2;
    delete[] tempconnf;

  }  // end for ilev

  // clean
  delete[] tempconn;
  delete[] tempconn_el;
  delete[] temp_material;
  delete[] temp_procid;
  H5Gclose(subgroup_id2);
  H5Gclose(subgroup_id);
  // =+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  CONN/FEM2 --> BOUNDARY MESH
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   subgroup_id = H5Gcreate(
//       file, "/ELEMS27/FEM1", H5P_DEFAULT
// #if HDF5_VERSIONM != 1808
//       ,
//       H5P_DEFAULT, H5P_DEFAULT
// #endif
//   );
//   dimsf[0] = 4;
//   dimsf[1] = 1;
//     ndofm[0] =  _geomel.n_q[1];  ndofm[1] = _geomel.n_l[1];  ndofm[2] =  _geomel._n_sides[1]; ndofm[3] = 1;
// //   ndofm[0] = NDOF_FEMB;
// //   ndofm[1] = 0;
// //   ndofm[2] = NDOF_PB;
//   status = _mgutils.print_Ihdf5(file, "/ELEMS27/FEM1/NDOFEM", dimsf, ndofm);
//   // No elements
//   dimsf[0] = _n_levels;
//   dimsf[1] = 1;
//   status = _mgutils.print_Ihdf5(file, "/ELEMS27/FEM1/NExLEV", dimsf, &n_elements_lev[_n_levels]);
// 
//   // Global Surface mesh at all levels  packaging data  (volume) ----------
//   tempconn = new int[bd_n_elements * NDOF_FEMB];
//   for (int ielem = 0; ielem < bd_n_elements; ielem++) {
//     for (int inode = 0; inode < NDOF_FEMB; inode++) {
//       tempconn[inode + ielem * NDOF_FEMB] = v_inv_nd[bd_elem_sto[v_elb[ielem].second][inode + 1]];
//     }
//   }
//   // global boundary mesh hdf5 sorage ---------------------
//   dimsf[0] = bd_n_elements * NDOF_FEMB;
//   dimsf[1] = 1;
//   status = _mgutils.print_Ihdf5(file, "/ELEMS27/FEM1/MSH", dimsf, tempconn);
// 
//   dimsf[0] = _n_subdomains * _n_levels + 1;
//   dimsf[1] = 1;
//   status = _mgutils.print_Ihdf5(file, "/ELEMS27/FEM1/OFF_EL", dimsf, bd_off_el);
//   //  Level surface mesh connectivity ---------------------
//   int tot_el = 0;
//   for (int ilev = 0; ilev < _n_levels; ilev++) {
//     tot_el += n_elements_lev[ilev + _n_levels];
//     int* tempconnf = new int[tot_el * NDOF_FEMB];
//     int ltot = 0;
//     for (int isubdom = 0; isubdom < _n_subdomains; isubdom++) {
//       for (int iel = bd_off_el[ilev + isubdom * _n_levels]; iel < bd_off_el[ilev + 1 + isubdom * _n_levels];
//            iel++) {
//         for (int inode = 0; inode < NDOF_FEMB; inode++) {
//           tempconnf[ltot * NDOF_FEMB + inode] = tempconn[iel * NDOF_FEMB + inode];
//         }
//         ltot++;
//       }
//     }
//     // partial print -------------------------
//     dimsf[0] = ltot * NDOF_FEMB;
//     dimsf[1] = 1;
//     name.str("");
//     name << "/ELEMS27/FEM1/MSH" << ilev;
//     status = _mgutils.print_Ihdf5(file, name.str(), dimsf, tempconnf);
//     // clean
//     delete[] tempconnf;
//   }
//   delete[] tempconn;
//   H5Gclose(subgroup_id);
  H5Gclose(group_id);
  // file closure
  H5Fclose(file);

  return;
}

// =============================================================
/// This function prints (hdf5 format) the prolongation matrix structures
void MGGenCase::print_op_h5(  // hdf5 version  1.8.8
    std::string name,         // file name
    int n_nodes_row,          // # of fine grid nodes
    int n_nodes_cln,          // # of coarse grid nodes
    int count_q,              // # of entries in Prol_q
    int* Res_q,               // compressed pos for prolongation Operator
    double* values_q,         // compressed prolongation Operator
    int* len_q,               // row length
    int* len_qoff,            // offset row len
    int qq                    // mode linear-quadratic
) {
  // ====================================================

  std::ostringstream mode;
  if (qq <10) {
    mode << "0" << qq;
  } else {
    mode << qq;
  }

  hid_t status = 0;
  // file
  hid_t fileR = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  // Prolongation dimensions
  std::ostringstream name0;
  name0 << "DIM" << mode.str().c_str();
  hsize_t dimsf[2];
  dimsf[0] = 2;
  dimsf[1] = 1;
  int rowcln[2];
  rowcln[0] = n_nodes_row;
  rowcln[1] = n_nodes_cln;
  status = _mgutils.print_Ihdf5(fileR, name0.str().c_str(), dimsf, rowcln);
  assert(status == 0);
  // node offset
  if (Res_q != NULL) {
    // print matrix position of P
    std::ostringstream name1;
    name1 << "POS" << mode.str().c_str();
    dimsf[0] = count_q;
    status = _mgutils.print_Ihdf5(fileR, name1.str().c_str(), dimsf, Res_q);
  }
  if (values_q != NULL) {
    // print values
    std::ostringstream name2;
    name2 << "VAL" << mode.str().c_str();
    status = _mgutils.print_Dhdf5(fileR, name2.str().c_str(), dimsf, values_q);
  }
  // print row length of R
  std::ostringstream name3;
  name3 << "LEN" << mode.str().c_str();
  dimsf[0] = n_nodes_row + 1;
  status = _mgutils.print_Ihdf5(fileR, name3.str().c_str(), dimsf, len_q);
  // print row length of R
  std::ostringstream name4;
  name4 << "OFFLEN" << mode.str().c_str();
  status = _mgutils.print_Ihdf5(fileR, name4.str().c_str(), dimsf, len_qoff);
  //   clean
  H5Fclose(fileR);

  return;
}






// ====================================================
// Print mesh in med format:
// conversion form libmesh to med format
void MGGenCase::print_med(
    int Level,                                                              // Level
    BoundaryMesh& bd_msh0,                                                  // boundary coarse mesh
    Mesh& msh0,                                                             // coarse mesh
    int n_groups, int* group_id_names, std::vector<std::pair<int, int>> v,  // node order
    std::vector<std::pair<int, int>> v_elb, int* nod_flag,                  // bc node
    int* mat_flag,                                                          // mat element
    std::string filename                                                    // filename
) {
// ==================================================
#ifdef HAVE_MED
  // name file directory
  std::string input_dir = _mgutils._mesh_dir;
  std::ostringstream name;
  name << input_dir << filename << "_gen.med";

  // coarse  mesh
  int n_nodes = msh0.n_nodes();    // from mesh
  int n_elements = msh0.n_elem();  // from mesh
  //  boundary coarse mesh
  int bd_n_nodes = bd_msh0.n_nodes();    // from boundary mesh
  int bd_n_elements = bd_msh0.n_elem();  // from boundary mesh

  // cordinates
  double* coord;
  coord = new double[n_nodes * _dim];
  // connectivity (vol+boundary)
  mcIdType* conn;
  conn = new mcIdType[n_elements * NDOF_FEM];
  mcIdType* conn_bd;
  conn_bd = new mcIdType[bd_n_elements * NDOF_FEMB];
  int* elem_bd_id2;
  elem_bd_id2 = new int[bd_n_elements];
  // element (nodes for element)
  int n_nodes_el;
  int bd_n_nodes_el;
  int count_eb = 0;

// med->libmesh map  second order only
#if ELTYPE == 27
#if DIMENSION == 3
  const unsigned int nodesinv[] = {4, 7,  3,  0,  5,  6, 2,  1,  19, 15, 11, 12, 17, 14,
                                   9, 13, 16, 18, 10, 8, 24, 25, 23, 20, 21, 22, 26};
  const unsigned int nodesinvbd[] = {3, 0, 1, 2, 7, 4, 5, 6, 8};
#endif
#if DIMENSION == 2
  const unsigned int nodesinv[] = {3, 0, 1, 2, 7, 4, 5, 6, 8};

  const unsigned int nodesinvbd[] = {0, 2, 1};
#endif
#if DIMENSION == 1
  const unsigned int nodesinv[] = {0, 2, 1};

  const unsigned int nodesinvbd[] = {0};
#endif
#endif
#if ELTYPE == 10
#if DIMENSION == 3
  const unsigned int nodesinv[] = {2, 3, 1, 0, 9, 8, 5, 6, 7, 4};
  const unsigned int nodesinvbd[] = {1, 2, 0, 4, 5, 3};
#endif
#if DIMENSION == 2
  const unsigned int nodesinv[] = {1, 2, 0, 4, 5, 3};
  const unsigned int nodesinvbd[] = {0, 2, 1};
#endif
#if DIMENSION == 1
  const unsigned int nodesinv[] = { 0, 2, 1 } const unsigned int nodesinvbd[] = {0};
#endif
#endif
  //  mesh (volume) -----------------------------------------------------
  int* map_elem = new int[n_elements];
  int n_element_top = 0;
  Mesh::const_element_iterator it_tr = msh0.elements_begin();
  const Mesh::const_element_iterator end_tr = msh0.elements_end();
  for (; it_tr != end_tr; ++it_tr) {
    Elem* elem = *it_tr;  // element
    int lev = elem->level();
    if (lev == Level) {
      int id_el = n_element_top;
      map_elem[id_el] = elem->id();  // element id
      n_nodes_el = elem->n_nodes();  // number of element nodes
      for (int inode = 0; inode < n_nodes_el; inode++) {
        int knode = elem->node_id(nodesinv[inode]);  // global node through map
        conn[id_el * n_nodes_el + inode] = knode;    // connectivity

        // coordinates storage
        for (int idim = 0; idim < _dim; idim++) {
          double xyz = msh0.point(knode)(idim);
          coord[knode * _dim + idim] = xyz;  // med cordinate tuple
        }
      }
      n_element_top++;
    }
    //  mesh (boundary) ----------------------------------------
    if (lev == Level) {
      for (int s = 0; s < (int)elem->n_sides(); s++) {
        if (elem->neighbor_ptr(s) == NULL) {
          std::unique_ptr<Elem> side(elem->build_side_ptr(s));  // face element
          bd_n_nodes_el = (int)side->n_nodes();           // face nodes
          int min = 100000;
          for (int ns = 0; ns < bd_n_nodes_el; ns++) {
            conn_bd[count_eb * bd_n_nodes_el + ns] = side->node_id(nodesinvbd[ns]);
            if (min > nod_flag[conn_bd[count_eb * bd_n_nodes_el + ns]]) {
              min = nod_flag[conn_bd[count_eb * bd_n_nodes_el + ns]];
            }
          }
          elem_bd_id2[count_eb] = min;
          count_eb++;  // counter
        }
      }
    }
  }

  bd_n_elements = count_eb;
  assert(bd_n_elements == count_eb);  // check
  std::cout << " Printing med file " << name.str().c_str() << "\n for bc and mat: nodes =" << n_nodes
            << "; elements =" << n_elements << "; n_nodes_el= " << n_nodes_el << std::endl;

  // MED mesh *************************************************

  // MEDCouplingUMesh mesh connectivity (volume)
  MEDCoupling::MEDCouplingUMesh* mesh1 = MEDCoupling::MEDCouplingUMesh::New("Mesh_1", _dim);
  mesh1->allocateCells(n_element_top);
  for (int i = 0; i < n_element_top; i++) {
    //         if (MED_EL_TYPE!=n_nodes_el) {
    //             std::cout<<"Attention!!! MED_EL_TYPE!=n_nodes_el in MGGenCase::print_med"<<std::endl;
    //             std::cout<<"Usually this means that your mesh is not linear or bi-quadratic"<<std::endl;
    //             abort();
    //         }
    mesh1->insertNextCell(MED_EL_TYPE, n_nodes_el, conn + i * n_nodes_el);
  }
  mesh1->finishInsertingCells();

  // MEDCouplingUMesh Mesh connectivity (boundary)
  MEDCoupling::MEDCouplingUMesh* mesh2 = MEDCoupling::MEDCouplingUMesh::New("Mesh_1", DIMENSION - 1);
  mesh2->allocateCells(bd_n_elements);
  for (int i = 0; i < bd_n_elements; i++) {
    mesh2->insertNextCell(MED_EL_BDTYPE, bd_n_nodes_el, conn_bd + i * bd_n_nodes_el);
  }
  mesh2->finishInsertingCells();

  // coord (same node set for both meshes)
  MEDCoupling::DataArrayDouble* coordarr = MEDCoupling::DataArrayDouble::New();
  coordarr->alloc(n_nodes, _dim);
  std::copy(coord, coord + n_nodes * _dim, coordarr->getPointer());
  mesh1->setCoords(coordarr);
  mesh2->setCoords(coordarr);

  // Setting MEDCouplingUMesh into MEDFileUMesh
  MEDCoupling::MEDFileUMesh* mm = MEDCoupling::MEDFileUMesh::New();
  mm->setName("Mesh_1");  // name needed to be non empty
  mm->setDescription("Description Mesh_1");
  mm->setCoords(mesh1->getCoords());
  mm->setMeshAtLevel(0, mesh1, false);
  mm->setMeshAtLevel(-1, mesh2, false);

  // Volume Groups

  std::map<int, std::vector<int>> vol_group_nodes;
  std::map<int, int> vol_group;
  for (int i = 0; i < n_element_top; i++) {
    vol_group[mat_flag[map_elem[i]]]++;
    vol_group_nodes[mat_flag[map_elem[i]]].push_back(i);
  }

  int n_vol_group = vol_group.size();
  std::vector<const MEDCoupling::DataArrayIdType*> gr_vol(n_vol_group);
  MEDCoupling::DataArrayIdType** g_vol = new MEDCoupling::DataArrayIdType*[n_vol_group];

  int js = 0;
  // defining the vol group data to store
  std::map<int, std::vector<int>>::iterator it_vol;
  for (it_vol = vol_group_nodes.begin(); it_vol != vol_group_nodes.end(); ++it_vol) {
    int igroup = it_vol->first;
    int is = it_vol->second.size();
    g_vol[js] = MEDCoupling::DataArrayIdType::New();
    g_vol[js]->alloc(is, 1);
    std::ostringstream name_p;
    name_p << igroup;
    g_vol[js]->setName(name_p.str().c_str());
    int* val1 = new int[is];
    for (int iv = 0; iv < is; iv++) { val1[iv] = it_vol->second[iv]; }
    std::copy(val1, val1 + is, g_vol[js]->getPointer());
    delete[] val1;
    gr_vol[js] = g_vol[js];
    js++;
  }
  // inserting  the volume groups into the med-file
  mm->setGroupsAtLevel(0, gr_vol, false);

  // Boundary group ******************************************
  // Finding bd_group and  bd_group_nodes map
  std::map<int, std::vector<int>> bd_group_nodes;
  std::map<int, int> bd_group;
  for (int i = 0; i < bd_n_elements; i++) {
    bd_group[elem_bd_id2[i]]++;
    bd_group_nodes[elem_bd_id2[i]].push_back(i);
  }
  int n_bd_group = bd_group.size();
  // group vector
  std::vector<const MEDCoupling::DataArrayIdType*> gr_bd(n_bd_group);
  MEDCoupling::DataArrayIdType** g_bd = new MEDCoupling::DataArrayIdType*[n_bd_group];

  js = 0;
  // defining the  group data to store
  std::map<int, std::vector<int>>::iterator it;
  for (it = bd_group_nodes.begin(); it != bd_group_nodes.end(); ++it) {
    int igroup = it->first;
    int is = it->second.size();
    //     std::cout << igroup << "\n";
    g_bd[js] = MEDCoupling::DataArrayIdType::New();
    g_bd[js]->alloc(is, 1);
    std::ostringstream name_p;
    name_p << igroup;
    g_bd[js]->setName(name_p.str().c_str());
    int* valb1 = new int[is]; /*int icount=0; */
    for (int iv = 0; iv < is; iv++) { valb1[iv] = it->second[iv]; }
    std::copy(valb1, valb1 + is, g_bd[js]->getPointer());
    delete[] valb1;
    gr_bd[js] = g_bd[js];
    js++;
  }
  // insert the boundary groups into the med-file
  mm->setGroupsAtLevel(-1, gr_bd, false);

  // Printing Group and Family
  std::cout << "\n \n =================================== ";
  std::cout << "\n Group Names -> Families : \n";
  std::map<std::string, std::vector<std::string>> a = (mm->getGroupInfo());
  std::map<std::string, std::vector<std::string>>::iterator ita;
  for (ita = a.begin(); ita != a.end(); ++ita) {
    std::string igroup = ita->first;
    int is = ita->second.size();
    std::cout << "\n " << igroup << "-> ";
    for (int i = 0; i < is; i++) {
      std::string a2 = ita->second[i];
      std::cout << a2 << "  ";
    }
  }
  std::cout << "\n \n =================================== ";
  std::cout << " \n Families -> Groups id: \n";
  std::map<std::string, mcIdType> fama = mm->getFamilyInfo();
  std::map<std::string, mcIdType>::iterator itfa;
  for (itfa = fama.begin(); itfa != fama.end(); ++itfa) {
    std::string igroup = itfa->first;
    std::cout << "\n " << igroup << "-> ";
    mcIdType a2 = itfa->second;
    std::cout << a2 << "  ";
  }
  std::cout << "\n\n ";
  mm->write(name.str().c_str(), 2);

  // clean
  delete[] map_elem;
  delete[] coord;
  fama.clear();
  a.clear();
  delete[] conn;
  coordarr->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
  mm->decrRef();

#else
  std::cout << "\n \n MGGenCase::Print_med you don't have MED library included in your project\n\n";
#endif

  return;
}

void MGGenCase::print_MedToMg(
    int Level,                                                              // Level
    BoundaryMesh& bd_msh0,                                                  // boundary coarse mesh
    Mesh& msh0,                                                             // coarse mesh
    int n_groups, int* group_id_names, std::vector<std::pair<int, int>> v,  // node order
    int* v_inv_nd, std::vector<std::pair<int, int>> v_elb,
    int* nod_flag,  // bc node
    int* mat_flag,  // mat element
    int* off_el, int* v_inv_el, std::vector<std::pair<int, int>> v_el, std::vector<int> ElementsPerLevel,
    int** elem_sto,
    std::string filename  // filename
) {
  // ==================================================
  // THIS ROUTINE IS USED TO PRINT MED MESH FILES WITH FIELDS THAT HELP CREATING INTERFACES
  // IN PARTICULAR WE SAVE THE PROC BY PROC MESH PARTITIONING AND THE NODE RENUMBERING FOR PARALLEL
  // DISTRIBUTION

  std::cout << "\033[1;31m\n............................................... \
                        \n PRINTING THE MED FILE FOR COUPLING INTERFACES \
                        \n...............................................\n\033[0m\n";

#ifdef HAVE_MED
  // name file directory
  std::string input_dir = _mgutils._mesh_dir;
  std::string info_name = input_dir + filename + "_MedToMg" + ".med";

  //  boundary coarse mesh
  int bd_n_nodes = bd_msh0.n_nodes();    // from boundary mesh
  int bd_n_elements = bd_msh0.n_elem();  // from boundary mesh
  mcIdType* conn_bd;
  conn_bd = new mcIdType[bd_n_elements * NDOF_FEMB];
  int* elem_bd_id2;
  elem_bd_id2 = new int[bd_n_elements];
  // element (nodes for element)
  int n_nodes_el;
  int NumBdNodes;
  int count_eb = 0;

// MAP BETWEEN MED AND LIBMESH ELEMENT NODE NUMBERING
#if ELTYPE == 27
#if DIMENSION == 3
  const unsigned int LibToMed[] = {7, 4,  5,  6,  3,  0,  1,  2,  19, 16, 17, 18, 11, 8,
                                   9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};

  const unsigned int LibToMed_bd[] = {3, 0, 1, 2, 7, 4, 5, 6, 8};
#endif
#if DIMENSION == 2
  const unsigned int LibToMed[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const unsigned int LibToMed_bd[] = {0, 1, 2};
#endif
#if DIMENSION == 1
  const unsigned int LibToMed[] = {0, 1, 2};
  const unsigned int LibToMed_bd[] = {0};
#endif
#endif
#if ELTYPE == 10
#if DIMENSION == 3
  const unsigned int LibToMed[] = {2, 3, 1, 0, 9, 8, 5, 6, 7, 4};
  const unsigned int LibToMed_bd[] = {1, 2, 0, 4, 5, 3};
#endif
#if DIMENSION == 2
  const unsigned int LibToMed[] = {1, 2, 0, 4, 5, 3};
  const unsigned int LibToMed_bd[] = {0, 2, 1};
#endif
#if DIMENSION == 1
  const unsigned int LibToMed[] = {0, 2, 1};
  const unsigned int LibToMed_bd[] = {0};
#endif
#endif

  /*  ELEMENT SORTING FOR MULTI PROC AND MULTILEVEL DOMAIN DECOMPOSITION
   *          PROC 0              PROC 1                  PROC N
   *  | .. :::: ======== | .. :::: ======== [...] | .. :::: ======== |
   *    L0  L1     LN      L0  L1     LN            L0  L1     LN
   *  --------------------------------------------------------------->
   *                INCREASING ELEMENT NUMBERING ORDER
   *
   *  | ------------------ A ---------------------| -- B -- |-- C --|
   *
   */

  for (int LEVEL = 0; LEVEL <= Level; LEVEL++) {
    Mesh::const_element_iterator it_tr = msh0.level_elements_begin(LEVEL);
    const Mesh::const_element_iterator end_tr = msh0.level_elements_end(LEVEL);

    int n_elements = msh0.n_elem();  // from mesh
    int* map_elem = new int[n_elements];

    int n_element_top = 0;

    MEDCoupling::DataArrayDouble* Volcells = MEDCoupling::DataArrayDouble::New();
    Volcells->alloc(ElementsPerLevel[LEVEL], 1);
    double* CellMap = const_cast<double*>(Volcells->getPointer());
    MEDCoupling::DataArrayDouble* ProcArray = MEDCoupling::DataArrayDouble::New();
    ProcArray->alloc(ElementsPerLevel[LEVEL], 1);
    double* ProcMap = const_cast<double*>(ProcArray->getPointer());

    mcIdType* conn = new mcIdType[ElementsPerLevel[LEVEL] * NDOF_FEM];
    std::vector<int> Nodes;

    for (; it_tr != end_tr; ++it_tr) {
      // LOOP OVER THE MESH ELEMENTS ------------------------
      Elem* elem = *it_tr;      // element
      int lev = elem->level();  // actual level of the mesh element
      int id_el = n_element_top;

      n_nodes_el = elem->n_nodes();  // number of element nodes

      if (LEVEL == Level) map_elem[id_el] = elem->id();
      int orig_el_id = v_inv_el[elem->id()];
      CellMap[id_el] = orig_el_id;
      ProcMap[id_el] = elem_sto[orig_el_id][NDOF_FEM + 2];

      // VOLUME
      for (int inode = 0; inode < n_nodes_el; inode++) {
        int knode = elem->node_id(LibToMed[inode]);  // global node through map
        conn[id_el * n_nodes_el + inode] = knode;    // element connectivity
        Nodes.push_back(knode);
      }
      n_element_top++;

      // BOUNDARY
      for (int s = 0; s < (int)elem->n_sides(); s++) {
        // loop over element sides
        if (elem->neighbor_ptr(s) == NULL) {
          // if neighbor_ptr == null then the side is on boundary
          std::unique_ptr<Elem> side(elem->build_side_ptr(s));  // face element
          NumBdNodes = (int)side->n_nodes();              // face nodes
          int min = 100000;
          for (int ns = 0; ns < NumBdNodes; ns++) {
            conn_bd[count_eb * NumBdNodes + ns] =
                side->node_id(LibToMed_bd[ns]);  // connectivity of boundary elements
            if (min > nod_flag[conn_bd[count_eb * NumBdNodes + ns]]) {
              min = nod_flag[conn_bd[count_eb * NumBdNodes + ns]];  // we set the boundary flag of the node ->
                                                                    // group id
            }
          }
          elem_bd_id2[count_eb] =
              min;     // boundary flag of the element -> it is equal to the lowest bd node flag of the side
          count_eb++;  // counter
        }
      }
    }  // END LOOP OVER MESH ELEMENTS -------------------------------------------------------
    bd_n_elements = count_eb;           // total number of boundary elements (sides)
    assert(bd_n_elements == count_eb);  // check

    // MED mesh *************************************************
    MEDCoupling::MEDCouplingUMesh* VolumeMesh =
        MEDCoupling::MEDCouplingUMesh::New("Mesh_Lev_" + std::to_string(LEVEL), _dim);
    VolumeMesh->allocateCells(n_element_top);
    for (int i = 0; i < n_element_top; i++) {
      VolumeMesh->insertNextCell(MED_EL_TYPE, n_nodes_el, conn + i * n_nodes_el);
    }
    VolumeMesh->finishInsertingCells();
    // MEDCouplingUMesh Mesh connectivity (boundary)
    MEDCoupling::MEDCouplingUMesh* BoundaryMesh =
        MEDCoupling::MEDCouplingUMesh::New("Mesh_Lev_" + std::to_string(LEVEL), _dim - 1);
    BoundaryMesh->allocateCells(bd_n_elements);
    for (int i = 0; i < bd_n_elements; i++) {
      BoundaryMesh->insertNextCell(MED_EL_BDTYPE, NumBdNodes, conn_bd + i * NumBdNodes);
    }
    BoundaryMesh->finishInsertingCells();

    std::set<int> s(Nodes.begin(), Nodes.end());
    Nodes.clear();
    int LevelMeshNodes = s.size();
    double* Coords = new double[_dim * LevelMeshNodes];

    MEDCoupling::DataArrayDouble* FinerConn = MEDCoupling::DataArrayDouble::New();
    FinerConn->alloc(LevelMeshNodes, 1);
    double* MapArray = const_cast<double*>(FinerConn->getPointer());

    for (int node = 0; node < LevelMeshNodes; node++) {
      MapArray[node] = v_inv_nd[node];
      for (int dim = 0; dim < _dim; dim++) { Coords[node * _dim + dim] = msh0.point(node)(dim); }
    }

    MEDCoupling::DataArrayDouble* coordarr = MEDCoupling::DataArrayDouble::New();
    coordarr->alloc(LevelMeshNodes, _dim);
    std::copy(Coords, Coords + LevelMeshNodes * _dim, coordarr->getPointer());
    VolumeMesh->setCoords(coordarr);
    BoundaryMesh->setCoords(coordarr);

    // Setting MEDCouplingUMesh into MEDFileUMesh
    MEDCoupling::MEDFileUMesh* mm = MEDCoupling::MEDFileUMesh::New();
    mm->setName("Mesh_Lev_" + std::to_string(LEVEL));  // name needed to be non empty
    mm->setDescription("Description Mesh_1");
    mm->setCoords(VolumeMesh->getCoords());
    mm->setMeshAtLevel(0, VolumeMesh, false);
    mm->setMeshAtLevel(-1, BoundaryMesh, false);

    if (LEVEL == Level) {
      // boundary and volume groups available only at finer level
      std::map<int, std::vector<int>> vol_group_nodes;
      std::map<int, int> vol_group;
      for (int i = 0; i < n_element_top; i++) {
        vol_group[mat_flag[map_elem[i]]]++;
        vol_group_nodes[mat_flag[map_elem[i]]].push_back(i);
      }

      int n_vol_group = vol_group.size();
      std::vector<const MEDCoupling::DataArrayIdType*> gr_vol(n_vol_group);
      MEDCoupling::DataArrayIdType** g_vol = new MEDCoupling::DataArrayIdType*[n_vol_group];

      int js = 0;
      // defining the vol group data to store
      std::map<int, std::vector<int>>::iterator it_vol;
      for (it_vol = vol_group_nodes.begin(); it_vol != vol_group_nodes.end(); ++it_vol) {
        int igroup = it_vol->first;
        int is = it_vol->second.size();
        g_vol[js] = MEDCoupling::DataArrayIdType::New();
        g_vol[js]->alloc(is, 1);
        std::ostringstream name_p;
        name_p << igroup;
        g_vol[js]->setName(name_p.str().c_str());
        int* val1 = new int[is];
        for (int iv = 0; iv < is; iv++) { val1[iv] = it_vol->second[iv]; }
        std::copy(val1, val1 + is, g_vol[js]->getPointer());
        delete[] val1;
        gr_vol[js] = g_vol[js];
        js++;
      }
      // inserting  the volume groups into the med-file
      mm->setGroupsAtLevel(0, gr_vol, false);

      // Boundary group ******************************************
      // Finding bd_group and  bd_group_nodes map
      std::map<int, std::vector<int>> bd_group_nodes;
      std::map<int, int> bd_group;
      for (int i = 0; i < bd_n_elements; i++) {
        bd_group[elem_bd_id2[i]]++;
        bd_group_nodes[elem_bd_id2[i]].push_back(i);
      }
      int n_bd_group = bd_group.size();
      // group vector
      std::vector<const MEDCoupling::DataArrayIdType*> gr_bd(n_bd_group);
      MEDCoupling::DataArrayIdType** g_bd = new MEDCoupling::DataArrayIdType*[n_bd_group];

      js = 0;
      // defining the  group data to store
      std::map<int, std::vector<int>>::iterator it;
      for (it = bd_group_nodes.begin(); it != bd_group_nodes.end(); ++it) {
        int igroup = it->first;
        int is = it->second.size();
        g_bd[js] = MEDCoupling::DataArrayIdType::New();
        g_bd[js]->alloc(is, 1);
        std::ostringstream name_p;
        name_p << igroup;
        g_bd[js]->setName(name_p.str().c_str());
        int* valb1 = new int[is]; /*int icount=0; */
        for (int iv = 0; iv < is; iv++) { valb1[iv] = it->second[iv]; }
        std::copy(valb1, valb1 + is, g_bd[js]->getPointer());
        delete[] valb1;
        gr_bd[js] = g_bd[js];
        js++;
      }
      // insert the boundary groups into the med-file
      mm->setGroupsAtLevel(-1, gr_bd, false);
    }

    if (LEVEL == 0) {
      mm->write(info_name, 2);
    } else {
      mm->write(info_name, 1);
    }

    // CELL FIELD STORING THE MED-TO-MG MESH CELL NUMBERING
    MEDCoupling::MEDCouplingFieldDouble* CellField =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
    CellField->setMesh(VolumeMesh);
    CellField->setArray(Volcells);
    CellField->setName("MG_cell_id_Lev_" + std::to_string(LEVEL));
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, CellField);

    // PROC FIELD STORING PROC ID
    MEDCoupling::MEDCouplingFieldDouble* ProcField =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
    ProcField->setMesh(VolumeMesh);
    ProcField->setArray(ProcArray);
    ProcField->setName("Proc_Lev_" + std::to_string(LEVEL));
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, ProcField);

    // MAP TO MG NODE NUMBERING
    MEDCoupling::MEDCouplingFieldDouble* FinerLevelIDS =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
    FinerLevelIDS->setMesh(VolumeMesh);
    FinerLevelIDS->setArray(FinerConn);
    FinerLevelIDS->setName("FinerLevelNodeIDS_Lev_" + std::to_string(LEVEL));
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, FinerLevelIDS);

    CellField->decrRef();
    ProcField->decrRef();
    FinerLevelIDS->decrRef();
    ProcArray->decrRef();
    FinerConn->decrRef();
    Volcells->decrRef();

    delete[] Coords;
    delete[] conn;
    delete[] map_elem;
    coordarr->decrRef();
    VolumeMesh->decrRef();
    BoundaryMesh->decrRef();
    mm->decrRef();
  }

  delete[] elem_bd_id2;
#else
  std::cout << "\n \n MGGenCase::Print_med you don't have MED library included in your project\n\n";
#endif
  std::cout << "\033[1;32m "
               "---------------------------------------------------------------------------------------------"
               "---------- \n";
  std::cout << " Mesh for interface creation printed in file " << info_name << std::endl;
  std::cout << " The domain is defined on " << _n_levels << " levels \n";
  std::cout << " For each level the following maps are printed: \n\
    Proc_Lev_<level>               ->  cell-wise field showing the processor id for each cell \n\
    MG_cell_id_Lev_<level>         ->  cell-wise field the cell-id within the global cell numbering (from level 0 to highest level)\n\
    FinerLevelNodeIDS_Lev_<level>  ->  point-wise field: it's a map from local (level) node numbering to highest level node numbering\n";
  std::cout << " --------------------------------------------------------------------------------------------"
               "----------- \033[0m\n";

  return;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;

