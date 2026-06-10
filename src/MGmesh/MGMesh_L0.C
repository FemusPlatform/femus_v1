// std lib ------------------>
#include <cmath>
#include <iomanip>
#include <sstream>

// conf includes ------------------
#include "Printinfo_conf.h"  // print config file
#include "Solverlib_conf.h"  // petsc conf

// class includes ---------------------
#include "MGGeom_L1_El.h"
#include "MGUtils_L0.h"
#include "MGMesh_L1_Extended.h"



// // #include "parallelM.h"
#ifdef HAVE_PETSCM
#include <mpi.h>  // This is needed in the constructor
#endif

#define NUMGEOMELS (1)


// ====================================================
/// This function is the mesh constructor
MGMesh::MGMesh(
    const ParallelObjectM& comm,
    MGUtils& mgutils_in,  // input file name class
    MGGeomElBase& geomel_in   // element geometry class input
    //   const double Lref    // reference length input
    )
    :  _GeomEl(geomel_in),     //  element geometry class
      _mgutils(mgutils_in),   // file name class
     _n_GeomEl(NUMGEOMELS),  // element geometry class
      _comm(comm) {           // communicator
TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::MGMesh (Constructor *****) \n  ");)      
// ===============================================================================================
  // read parameter from GeometrySetting.in
  // ==================================

// geometry parameters from  GeometrySetting.in  for rotation during mesh reading 
// // _theta[2]={0.,0.}; double _Coxyz[3]={0.,0.,0.};
// int err=_mgutils._dict_geom.size();
//  _theta = _mgutils._dict_geom["theta"];
//  _Coxyz[0]=  _mgutils._dict_geom["Ox"];
//  _Coxyz[1]=  _mgutils._dict_geom["Oy"];
//  _Coxyz[2]=  _mgutils._dict_geom["Oz"];
//  // _label[0,1,2] defined by  VART0,1,2 
//  _type_nonodes[0]=(int)(_mgutils._dict_geom["VART0"]);
//  _type_nonodes[1]=(int)(_mgutils._dict_geom["VART1"]);
//  _type_nonodes[2]=(int)(_mgutils._dict_geom["VART2"]);
//  for(int ktype=0;ktype<3;ktype++)_label[ktype]= _dict_label[_type_nonodes[ktype]];
//  // restart
//  _restart = stoi(_mgutils._dict_config["restart"]);
// if(_mgutils._dict_geom.size()-err) {
//   std::cout<< "MGMesh::MGMesh: Missing parameter in GeometrySetting.in"; std::exit(EXIT_FAILURE);}
//   
 
  
  
  // reference length
 _ord_FEM = new int[27];  for (int ik = 0; ik < 27; ik++) _ord_FEM[ik] = ik;
  if(_GeomEl.n_q[0]==27) {  int tmporder[27] ={ 0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,24,22,21,23,20,25,26};
       for (int ik = 0; ik < 27; ik++)  _ord_FEM[ik]=tmporder[ik]; 
 }
// 0 1 2 3 4 5 6 7 8 9 10 11 16 17 18 19 12 13 14 15 24 22 21 23 20 25 26
// _ord_FEM
  // set subdomain number
  int i;
#ifdef HAVE_PETSCM
  MPI_Comm_rank(MPI_COMM_WORLD, &i);
  _iproc = static_cast<int>(i);
P_INFO(std::cout << " ===== Mesh: ===== \n communicator: " << MPI_COMM_WORLD << std::endl;printf(" \n %d proc=: ", _iproc);)
#else
  _iproc = 0;
#endif
  this->read_c();  // reading
  // _Oxyz_pts[1]; //  normal at points
TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::MGMesh (Constructor *****) \n  ");)
  return;
}


// ================================================================================================
/// This function is the mesh destructor
MGMesh::~MGMesh() {  // ===========================================================================
TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::~MGMesh (Destructor *****) \n  ");)
// ================================================================================================
  clear();
  delete[] _NoElements;  delete[] _xyz;
  delete[] _xyzo;  delete[] _dxdydz;
  delete[] _NoNodes;

  delete[] _type_FEM;  delete[] _el_map;  //  delete[]_conn_map;
  delete[] _off_el;
  //   delete[] _off_nd;
  //   delete[] _off_nd1;
  delete[] _ord_FEM;  delete[] _dist;
//   delete[] _VolFrac;
//   delete[] _ctrl_dom;
  delete[] _el_neighbor;  delete[] _node_map;
  delete[] _NodeBDgroup;  delete[] _normalb;
TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::~MGMesh (Destructor *****) \n  ");)
return;
}

// ================================================================================================
/// This function destructs the mesh substructures
void MGMesh::clear() {  // ============================================
TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::clear \n  ");)
// ================================================================================================
  for (int imesh = 0; imesh < 3*_NoLevels; imesh++) delete[] _node_map[imesh];

  //   for(int  imesh =0; imesh<_NoLevels; imesh++)   delete []_dist; //[imesh];

 for(int type=0; type<3;type++) { delete[] _min_nd[type]; delete[] _max_nd[type];}

  for (int imesh = 0; imesh < _NoFamFEM; imesh++) {
    delete[] _el_map[imesh];
    delete[] _off_el[imesh];
    delete[] _NoElements[imesh];
  }

  for (int imesh = 0; imesh < _NoFamFEM; imesh++) {
    delete[] _el_neighbor[imesh];  // now _el_neighbor is allocated only for the volume family
  }
TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::clear \n  ");)
  return;
}


// ================================================================================================
/// This function computes the element center
void MGMesh::get_el_ctr(
    const int el_nnodes,   // element nodes  (<- input)=_GeomEl.n_q[bdry]
    const int /*bdry*/,        // zone vol/bd    (<- input)=(1)/(0) flag
    const double* xx_nds,  // element coords (-> output)
    double* el_xm          // element center (-> output)
    ) const {              // =======================================
TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::get_el_ctr  \n  ");)
// ================================================================================================
  for (int idim = 0; idim < _dim; idim++) {
    el_xm[idim] = 0.;  // zero

    for (int eln = 0; eln < el_nnodes; eln++) el_xm[idim] += xx_nds[eln + idim * el_nnodes];

    el_xm[idim] = el_xm[idim] / el_nnodes;  // normalization
  }
TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::get_el_ctr  \n  ");)
  return;
}

// ================================================================================================
/// This function gets the global node numbers for that element
/// and their coordinates
void MGMesh::get_el_nodes(
    const int el_nnodes,  // element nodes(<- input)= _GeomEl.n_q[bdry];
    const int bdry,       // zone vol/bd   (<- input)=(1)/(0) flag
    const int Level,      // level         (<- input)
    const int iel,        // element       (<- input)
    double xx[]           // element coords(-> output)
    ) const {             // =======================================
// ================================================================================================
  // element  definition
  for (int n = 0; n < el_nnodes; n++) {
    // get the global node number
    int el_conn = _el_map[bdry][(iel + _off_el[bdry][Level + _NoLevels * _iproc]) * el_nnodes + n];

    // get the element coordinades
    for (int idim = 0; idim < _dim; idim++) {
      const int indxn = n + idim * el_nnodes;
      xx[indxn] = _xyz[el_conn + idim * _NoNodes[_NoLevels - 1]];
    }
  }

  return;
}


// ================================================================================================
/// This function return the element connectivity and  coordinates
void MGMesh::get_el_nod_disp(
    const int bdry,   // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
    const int Level,  // level                (<- input)
    const int iel,    // element              (<- input)
                      //   int  el_conn[],   // element connectivity (<- intput)
    double xx[]       // element coordinates  (-> output)
    ) const {         // =================================
// ================================================================================================
  const int el_nnodes = _GeomEl.n_q[bdry];  // element nodes
                                            //   const int  offset   =_NoNodes[_NoLevels-1];
                                            //   // element  definition
                                            //   for(int  inode=0; inode<el_nnodes; inode++)    {
                                            //     //get the global node number
  // //     el_conn[inode] = _el_map[bdry][(iel+_off_el[bdry][Level+_NoLevels*_iproc])*el_nnodes+inode];
  //     // get the element coordinades
  //     for(int  idim=0; idim<_dim; idim++) {
  //       xx[inode+idim*el_nnodes] = _dxdydz[el_conn[inode]+idim*offset];
  // //       xx[inode+idim*el_nnodes] = _xyz[el_conn[inode]+idim*offset]+_dxdydz[el_conn[inode]+idim*offset];
  //
  // //        xx[inode+idim*el_nnodes] = _dxdydz[el_conn[inode]+idim*offset];
  //     }
  //   }

  // element  definition
  for (int n = 0; n < el_nnodes; n++) {
    // get the global node number
    int el_conn = _el_map[bdry][(iel + _off_el[bdry][Level + _NoLevels * _iproc]) * el_nnodes + n];

    // get the element coordinades
    for (int idim = 0; idim < _dim; idim++) {
      const int indxn = n + idim * el_nnodes;
      xx[indxn] = _dxdydz[el_conn + idim * _NoNodes[_NoLevels - 1]];
    }
  }

  return;
}

// ================================================================================================
/// This function return the element connectivity and  coordinates
void MGMesh::get_el_nod_conn(
    const int bdry,     // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
    const int Level,    // level                (<- input)
    const int iel,      // element              (<- input)
    int el_conn[],      // element connectivity (-> output)
    double xx[],        // element coordinates (-> output)
    int iproc) const {  // =================================
// ================================================================================================
  const int el_nnodes = _GeomEl.n_q[bdry];  // element nodes
  const int offset = _NoNodes[_NoLevels - 1];

  // element  definition
  for (int inode = 0; inode < el_nnodes; inode++) {
    // get the global node number
    el_conn[inode] = _el_map[bdry][(iel + _off_el[bdry][Level + _NoLevels * iproc]) * el_nnodes + inode];

    // get the element coordinades
    for (int idim = 0; idim < _dim; idim++) {
      xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
    }
  }

  return;
}

// ================================================================================================
/// This function returns the element connectivity and coordinates
void MGMesh::get_el_nod_conn(
    const int bdry,   // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
    const int Level,  // level                (<- input)
    const int iel,    // element              (<- input)
    int el_conn[],    // element connectivity (-> output)
    double xx[]       // element coordinates  (-> output)
    ) const {         // =================================
// ================================================================================================
  const int el_nnodes = _GeomEl.n_q[bdry];  // element nodes
  const int offset = _NoNodes[_NoLevels - 1];

  // element  definition
  for (int inode = 0; inode < el_nnodes; inode++) {
    // get the global node number
    el_conn[inode] = _el_map[bdry][(iel + _off_el[bdry][Level + _NoLevels * _iproc]) * el_nnodes + inode];

    // get the element coordinates
    for (int idim = 0; idim < _dim; idim++) {
      xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
    }
  }

  return;
}

// ================================================================================================
/// This function returns the element connectivity and coordinates
void MGMesh::get_el_nod_conn(
    const int bdry,   // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
    const int iel,    // element              (<- input)      -> this function gets the exact element position in the _el_map vector
    int el_conn[],    // element connectivity (-> output)
    double xx[]       // element coordinates  (-> output)
    ) const {         // =================================
// ================================================================================================
  const int el_nnodes = _GeomEl.n_q[bdry];  // element nodes
  const int offset = _NoNodes[_NoLevels - 1];

  // element  definition
  for (int inode = 0; inode < el_nnodes; inode++) {
    // get the global node number
    el_conn[inode] = _el_map[bdry][iel * el_nnodes + inode];

    // get the element coordinates
    for (int idim = 0; idim < _dim; idim++) {
      xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
    }
  }

  return;
}

// ================================================================================================
/// This function returns the element connectivity and coordinates
void MGMesh::get_el_nod_conn_global(
    const int bdry,   // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
    const int Level,  // level                (<- input)
    const int iel,    // element              (<- input)
    int el_conn[],    // element connectivity (-> output)
    double xx[]       // element coordinates  (-> output)
    ) const {         // =================================
// ================================================================================================
  const int el_nnodes = _GeomEl.n_q[bdry];  // element nodes
  const int offset = _NoNodes[_NoLevels - 1];

  // element  definition
  for (int inode = 0; inode < el_nnodes; inode++) {
    // get the global node number
    el_conn[inode] = _el_map[bdry][(iel) * el_nnodes + inode];

    // get the element coordinates
    for (int idim = 0; idim < _dim; idim++) {
      xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
    }
  }

  return;
}

// ================================================================================================
/// This function return the element neigbours
void MGMesh::get_el_neighbor(
    const int el_sides,  // element sides        (<- input) _GeomEl._n_sides[bdry];
    const int bdry,      // bd(1)/vol(0) flag    (<- input) (1)/(0) flag
    const int Level,     // level                (<- input)
    const int iel,       // element              (<- input)
    int el_neigh[],      // element connectivity (-> output)
    int iproc) const {   // =======================================================================
// ================================================================================================
  //   const int  el_sides= _GeomEl._n_sides[bdry]; // element nodes
  for (int iside = 0; iside < el_sides; iside++) {
    // get the global node number
    el_neigh[iside] = _el_neighbor[bdry][(iel + _off_el[bdry][Level + _NoLevels * iproc]) * el_sides + iside];
  }

  return;
}

// ================================================================================================
/// This function return the element neigbours
void MGMesh::get_el_neighbor(
    const int el_sides,  // element sides        (<- input) _GeomEl._n_sides[bdry];
    const int bdry,      // bd(1)/vol(0) flag    (<- input) (1)/(0) flag
    const int Level,     // level                (<- input)
    const int iel,       // element              (<- input)
    int el_neigh[]       // element connectivity (-> output)
    ) const {            // ======================================================================
// ================================================================================================
  //   const int  el_sides= _GeomEl._n_sides[bdry]; // element nodes
  for (int iside = 0; iside < el_sides; iside++) {
    // get the global node number
    el_neigh[iside] =
        _el_neighbor[bdry][(iel + _off_el[bdry][Level + _NoLevels * _iproc]) * el_sides + iside];
  }
  return;
}

// ================================================================================================
/// This function gets the element connectivity
void MGMesh::get_el_conn(
    const int el_nnodes,  // element nodes    (<- input) = _GeomEl.n_q[bdry];
    const int bdry,       // vol/surf         (<- input)=(1)/(0) flag
    const int Level,      // Level            (<- input)
    const int iel,        // element          (<- input)
    int el_conn[]         // connectivity map (-> output)
    ) const {             // =====================================================================
// ================================================================================================
  // get the global node number
  for (int n = 0; n < el_nnodes; n++) {
    el_conn[n] = _el_map[bdry][(iel + _off_el[bdry][Level + _NoLevels * _iproc]) * el_nnodes + n];
  }
  return;
}

// ================================================================================================
/// This function returns the coordinate vector
void MGMesh::xcoord(
    double cvect[], const int n_nodes,
    const int offset) const {  // ================================================================
// ================================================================================================     
  // computation xyz vector in variable order
  for (int i = 0; i < n_nodes; i++) cvect[i] = _xyz[i + offset];
}

// ================================================================================================
void MGMesh::nodesxyz(
    double xyzvect[],  // coord vector ->
    const int n_nodes  // node number  <-
    ) const {          // =========================================================================
// ================================================================================================      
  // computation xyz vector in node order
  for (int i = 0; i < n_nodes; i++) {
    xyzvect[i * 3] = _xyz[i];
    xyzvect[i * 3 + 1] = 0.;
    xyzvect[i * 3 + 2] = 0.;
// #if DIMENSION == 2
if(_dim==2)     xyzvect[i * 3 + 1] = _xyz[i + n_nodes];
// #endif
if(_dim==3){
    xyzvect[i * 3 + 1] = _xyz[i + n_nodes];
    xyzvect[i * 3 + 2] = _xyz[i + 2 * n_nodes];
}
// #endif
  }
}



// ================================================================================================
/// This function computes the subconnectivity
int MGMesh::sub_conn(
    int gl_conn[],         // global sub-connectivity map
    const int ifem,        // fem type
    const int ilev,        // Level
    const int n_points_el  // point number
    ) const {  // =================================================================================           
// ================================================================================================
  int n_elements = _NoElements[ifem][ilev];
  int type_family_in = _type_FEM[ilev + ifem * _NoLevels];

  int icount = 0;

  for (int ik = 0; ik < n_elements; ik++) {
    for (int inode = 0; inode < n_points_el; inode++) {
      gl_conn[icount] = _el_map[ifem][inode + ik * type_family_in];
      icount++;
    }
  }

P_INFO(  std::cout << icount << std::endl;)
  return icount;
}

// =============================================
/// This function computes the connectivity
void MGMesh::conn(
    int gl_conn[],       // global sub-connectivity map
    const int ifem,      // fem type
    const int indx_mesh  // mesh
    ) const {            // =============================================
  for (int ik = 0; ik < _NoElements[ifem][indx_mesh] * _type_FEM[indx_mesh]; ik++) {
    gl_conn[ik] = _el_map[ifem][ik];
  }

  return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGMesh::Translate(const int dir, NumericVectorM& x_old) {  // ------------------------------------

  const int flag_moving_mesh = (int)_mgutils._dict_geom["moving_mesh"];

  if (flag_moving_mesh) {
    // #ifdef FINE_MOVEMENT
    /// E) mesh update
    double ww = 1.;
    const int n_nodes = _NoNodes[_NoLevels - 1];
    const int offsetp = dir * n_nodes;

    for (int inode = 0; inode < n_nodes; inode++) {
      _xyz[inode + offsetp] = _xyzo[inode + offsetp] + ww * x_old(inode);
    }

    // #endif
    // #ifdef COARSE_MOVEMENT
    //     MoveMesh(_NoLevels-1);
    // // disp[_NoLevels-1]->zero();
    //     for(int inode=0; inode<n_nodes; inode++) {
    //       _mgmesh._xyz[inode+offsetp] += (*disp[_NoLevels-1])(inode);
    //       _mgmesh._dxdydz[inode+offsetp]= (*disp[_NoLevels-1])(inode);
    //     }
    //     MoveMesh(_NoLevels-1);
    // #endif
  }

  // ==============================================================
  return;
}





// ================================================================================================
// This function sets the geometric normals based on the Group regions
void MGMesh::set_node_normal() 
// ================================================================================================
{TRACKING_FUN( printf(" --> MGMesh.C: MGMesh::set_node_normal  \n  "); )
// ================================================================================================  
  // -------------------------- Geometry -------------------------------------------------------------
  // top level mesh info
  const int Level = _NoLevels - 1;             // top Level always
  
//   const int offset = _NoNodes[_NoLevels - 1];  // mesh nodes
  //  element info
  double x_m[3];                  // element central point coords
  double xx_qnds[27 * 3];   // volume point coords
  double xxb_qnds[27 * 3];  // boundary point coords
  int el_conn[27];                  // node volume/boundary connectivity
  int el_neigh[27];                 // element connectivity
  // boundary
  const int el_sides = _GeomEl._n_sides[0];  // boundary  element sides
  int elb_conn[9];                   // node volume/boundary connectivity
  double normal[3 * 3];      // normal+tg to the boundary
  const int el_nnodes = _GeomEl.n_q[0]; 
  const int el_nnodes_bd = _GeomEl.n_q[1]; 
  // loop over procs+elements at top level Level(_NoLevels-1)
  for (int pr = 0; pr < _n_subdom; pr++) {
    const int nel_e = _off_el[0][Level + _NoLevels * pr + 1];  // start element
    const int nel_b = _off_el[0][Level + _NoLevels * pr];      // stop element
    for (int iel = 0; iel < (nel_e - nel_b); iel++) {
      // ---------------------------------------------------------------------------------
      // Element Connectivity (el_conn)  and coordinates (xx_qnds)
      get_el_nod_conn(0, Level, iel, el_conn, xx_qnds, pr);
      get_el_neighbor(el_sides, 0, Level, iel, el_neigh, pr);

      // central point for check extern normal
      for (int idim = 0; idim < _dim; idim++) {  // quad loop entities (vector)
        x_m[idim] = 0.;
        for (int d = 0; d < el_nnodes; d++) x_m[idim] += xx_qnds[idim * el_nnodes + d] / el_nnodes;
      }

      // Loop over the boundary
      for (int iside = 0; iside < el_sides; iside++) {
        if (el_neigh[iside] <0) {
          // setup boundary element  ----------------------------------------------------------------
          for (int lbnode = 0; lbnode <el_nnodes_bd; lbnode++) {          // quad quantities
            int lnode = _GeomEl._surf_top(lbnode + el_nnodes_bd * iside);  // local nodes
            //             sur_toply[lbnode] = lnode;        // lbnode -> lnode
            elb_conn[lbnode] = el_conn[lnode];  // connctivity el_conn->elb_conn

            for (int idim = 0; idim < _dim; idim++) {  // coordinates    xxb_qnds[x_EDGE;y_EDGE;z_EDGE]
              xxb_qnds[idim * el_nnodes_bd + lbnode] = xx_qnds[idim * el_nnodes + lnode];
            }
          }
          get_normal_b(xxb_qnds, x_m, normal);  // normal,tg1,tg2
          //         int dir_maxnormal = (fabs(normal[0]) > fabs(normal[1])) ? 0 : 1 ;
          //         dir_maxnormal = (fabs(normal[dir_maxnormal]) > fabs(normal[DIMENSION - 1])) ?
          //         dir_maxnormal : DIMENSION - 1; std::cout<<"
          //         .........................................................\n Element "<<iel<<" side "<<
          //         iside <<std::endl;
          //         set_bc_matrix(dir_maxnormal, sur_toply, el_ndof, elb_ndof, elb_ngauss, normal, el_conn);

          for (int lbnode = 0; lbnode < el_nnodes_bd; lbnode++) {
            int in = elb_conn[lbnode];
            int ibc = _NodeBDgroup[in];
            int in_mid = elb_conn[el_nnodes_bd - 1];
            int ibc_mid = _NodeBDgroup[in_mid];
            // Note; ibc == ibc_mid geometric -> normals based on the Group regions
            for (int ivar = 0; ivar < _dim; ivar++)
              if (ibc == ibc_mid) _normalb[in * _dim + ivar] = normal[ivar];
//             printf(" Mesh normal : %d %d  %f %f %d %d \n", _iproc, in, normal[0], normal[1], ibc, ibc_mid);

          }  // int  lbnode
        }    // el_neigh[iside] == -1
      }      // int  iside -----------------------------  End Boundary -------------------------------------
    }        //  =============== End of element loop =============================================
TRACKING_FUN( printf(" <-- MGMesh.C: MGMesh::set_node_normal  \n  ");)
    return;
  }
// ----------------------------------------------------------------------------------------------------
P_INFO( std::cout << " Mesh normal  computed  " << std::endl;)
  return;
}



// ================================================================================================
void MGMesh::get_normal_b(
    const double* xx,    // all surface coordinates <-
    const double x_c[],  // point inside the element <-
    double normal_tg[]   //  normal [0;_dim-1]+tg1[_dim;...] ->
    ) const {            // ======================================================================
// ================================================================================================
  const int el_nnodes_bd = _GeomEl.n_q[1]; 
if(_dim == 2){  //  2D -----------------------------------------------
                    //  The surface elements are such that when you go from the 1st to
                    //   the 2nd point, the outward normal is to the RIGHT
                    //  normal= normal_tg[0;_dim-1] tangent = normal_tg[_dim;2*_dim-1]
  int idx_tg = _dim;
  for (int i = 0; i < _dim; i++) normal_tg[i + idx_tg] = xx[1 + i * el_nnodes_bd] - xx[0 + i * el_nnodes_bd];
  //  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_tg[0] = 0. * normal_tg[0 + idx_tg] + 1. * normal_tg[1 + idx_tg];
  normal_tg[1] = -1. * normal_tg[0 + idx_tg] + 0. * normal_tg[1 + idx_tg];
  // if the sign is not correct then reverse it
  double sp = 0.;
  for (int kdim = 0; kdim < _dim; kdim++) sp += normal_tg[kdim] * (x_c[kdim] - xx[0 + kdim * el_nnodes_bd]);
  if (sp > 0.) {  // std::cout << " Normal inverted ! ------------------------------------  \n";
    for (int kdim = 0; kdim < _dim * _dim; kdim++) normal_tg[kdim] *= -1.;
  }
}
else if(_dim==3){  // 3D ---------------------------------------------------
  //     double tg03[DIMENSION];  // tangent plane
  //   the cross product of the two tangent vectors
  //  normal= normal_tg[0;_dim-1] tangent = normal_tg[_dim;2*_dim-1]
  //  tangent2 = normal_tg[_dim;2*_dim-1]
  int idx_tg1 = _dim;
  int idx_tg2 = 2 * _dim;
  for (int i = 0; i < _dim; i++) {
    normal_tg[i + idx_tg1] = xx[1 + i * el_nnodes_bd] - xx[0 + i * el_nnodes_bd];
    normal_tg[i + idx_tg2] = xx[2 + i * el_nnodes_bd] - xx[1 + i * el_nnodes_bd];
  }
  //   _mgutils.cross(tg01,tg03,normal_g);
  for (int i = 0; i < _dim; i++)
    normal_tg[i] = normal_tg[(i + 1) % 3 + idx_tg1] * normal_tg[(i + 2) % 3 + idx_tg2] -
                   normal_tg[(i + 2) % 3 + idx_tg1] * normal_tg[(i + 1) % 3 + idx_tg2];

  //     normal_tg[0]=normal_tg[1+idx_tg1]*normal_tg[2+idx_tg2]-normal_tg[2+idx_tg1]*normal_tg[1+idx_tg2];
  //     normal_tg[1]=normal_tg[2+idx_tg1]*normal_tg[0+idx_tg2]-normal_tg[0+idx_tg1]*normal_tg[2+idx_tg2];
  //     normal_tg[2]=normal_tg[0+idx_tg1]*normal_tg[1+idx_tg2]-normal_tg[1+idx_tg1]*normal_tg[0+idx_tg2];

  // if the sign is not correct then reverse it (normal toward outboundary)
  double sp = 0.;
  for (int kdim = 0; kdim < _dim; kdim++) sp += normal_tg[kdim] * (x_c[kdim] - xx[0 + kdim * el_nnodes_bd]);
  if (sp > 0) {
    for (int kdim = 0; kdim < _dim * _dim; kdim++) normal_tg[kdim] *= -1.;
  }
  for (int i = 0; i < _dim; i++)
    normal_tg[i + idx_tg2] = normal_tg[(i + 1) % 3] * normal_tg[(i + 2) % 3 + idx_tg1] -
                             normal_tg[(i + 2) % 3] * normal_tg[(i + 1) % 3 + idx_tg1];

}
  // normalization -----------------------
  for (int in = 0; in < _dim; in++) {
    double mm = 0.;
    for (int kdim = 0; kdim < _dim; kdim++) mm += normal_tg[kdim + in * _dim] * normal_tg[kdim + in * _dim];
    mm = sqrt(mm + 1.e-20);
    for (int idim = 0; idim < _dim; idim++) normal_tg[idim + in * _dim] /= mm;
  }
  return;
}
