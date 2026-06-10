#ifndef __mggencase__
#define __mggencase__

// std lib ----------------------
#include <string>

// configure includes ---------

// local includes -------------------
#include <libmesh/boundary_mesh.h>
#include <libmesh/mesh.h>

using namespace libMesh;

// Forwarding classes ----------
class MGUtils;
class MGGeomEl;
// class Mesh;
// class BoundaryMesh;
// ===========================================
//           MGGenCase class functions
// ===========================================

class MGGenCase {
 protected:
  // data ----------------------------------------------
  //   MGUtils& _mgutils;  ///<  MGUtils class pointer
  const Parallel::Communicator& _comm;
  int _bcmat;         ///<  boundary condition from input file
  MGGeomEl& _geomel;  ///<  MGGeomEl class pointer

  // Mesh
  int _dim;           ///< space dimension
  int _n_subdomains;  ///< # of subdomain (one for each processor)
  int _n_levels;      ///< # of levels (multigrid)
  
  int _n_local_faces;
  int _n_glob_faces;

  // node-element table parameters
  int _N_CHILD;   ///< # of fem children for each parent element
  int _N_NDV;     ///< # of param in the node table (nod_sto)
  int _dcl_nel;   ///< # of param in the element table  (elem_sto)
  int _dclb_nel;  ///< # of parameters in the boundary element table

  int** _ElPerProcPerLevel;

 public:
  MGUtils& _mgutils;  ///<  MGUtils class pointer

  // Destructor-Constructor ---------------------------------------
  MGGenCase(
      const Parallel::Communicator& comm_in, MGUtils& mgutils_in, MGGeomEl& mggeomel_in);  ///<  Constructor
  ~MGGenCase();                                                                            ///<  Destructor

  // functions -----------------------------------
  /// Case generator
  void GenCase();

  /// Compute sparse structures for  Operators
  void compute_and_print_MGOps(
      int* arg2, int arg3, int* arg5, std::vector<std::pair<int, int>> arg6, int** arg8, int* arg7,
      int* v_inv_el, int** arg9, int* min_nd[],int* max_nd[],int** elem_conn);

  /// ===========================================================================================
  /// MATRIX
   void compute_matrix(
      int* n_nodes_lev,                       // nodes for level
      int max_elnd,                           // max # of elements at each point
      int* off_el,                            // element offset
      std::vector<std::pair<int, int>> v_el,  // element ordering
      int** elem_sto,                         // element structure
      int* v_inv_nd,                          // node inverse ordering
      int** g_indexL,                         // node map
      int* min_nd[],                           // node offset
      int* max_nd[],
      int** elem_conn
  );                         
  
  
  /// Compute sparse structures for matrices
  int compute_mat_qq(int name,
      int Level, 
//       int Level_row, int Level_clmn, int ndof_ql, int ndof_ql_clmn, 
      int n_nodes, int max_elnd,
      int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto,
      int* v_inv_nd, int** g_indexL, int* min_nd,int* max_nd, 
      int offM[], 
      int aux[], int mem[], int Mat_q[], int len_q[], int len_qoff[]);

  int compute_mat_qq_ghost(int name,
      int Level, 
//       int Level_row, int Level_clmn, int ndof_ql, int ndof_ql_clmn, 
      int n_nodes, int max_elnd,
      int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto,
      int* v_inv_nd, int** g_indexL, int* min_nd,int* max_nd, 
      int offM[], 
      int aux[], int mem[], int Mat_q[], int len_q[], int len_qoff[],
      int** elem_conn);

//   // =======================================================
//   /// This function computes the quadratic matrix
//   int compute_mat_qk(
//       int Level1,        // level
//       int Level_row,     // row level
//       int Level_clmn,    // column level
//       int ndof_row_ql,   // row nodes
//       int ndof_clmn_ql,  // column nodes
//       int n_nodes,       // top nodes
//       int max_elnd, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
//       int** g_indexL, int* off_nd, int offM[], int aux[], int mem[], int Mat_q[], int len_q[],
//       int len_qoff[]);
// 
//   // =======================================================
//   /// This function computes the quadratic matrix
//   int compute_mat_lk(
//       int Level1,        // level
//       int Level_row,     // row level
//       int Level_clmn,    // column level
//       int ndof_row_ql,   // row nodes
//       int ndof_clmn_ql,  // column nodes
//       int n_nodes,       // top nodes
//       int max_elnd, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
//       int** g_indexL, int* off_nd, int offM[], int aux[], int mem[], int Mat_q[], int len_q[],
//       int len_qoff[]);
// 
//   // =======================================================
//   /// This function computes the quadratic matrix
//   int compute_mat_kq(
//       int Level1,        // level
//       int Level_row,     // row level
//       int Level_clmn,    // column level
//       int ndof_row_ql,   // row nodes
//       int ndof_clmn_ql,  // column nodes
//       int n_nodes,       // top nodes
//       int max_elnd, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
//       int** g_indexL, int* off_nd, int offM[], int aux[], int mem[], int Mat_q[], int len_q[],
//       int len_qoff[]);
//   // =======================================================
//   /// This function computes the quadratic matrix
//   int compute_mat_kl(
//       int Level1,        // level
//       int Level_row,     // row level
//       int Level_clmn,    // column level
//       int ndof_row_ql,   // row nodes
//       int ndof_clmn_ql,  // column nodes
//       int n_nodes,       // top nodes
//       int max_elnd, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
//       int** g_indexL, int* off_nd, int offM[], int aux[], int mem[], int Mat_q[], int len_q[],
//       int len_qoff[]);
//   //===================================
//   int compute_mat_kk(
//       int Level1,        // level
//       int Level_row,     // row level
//       int Level_clmn,    // column level
//       int ndof_row_ql,   // row nodes
//       int ndof_clmn_ql,  // column nodes
//       int n_nodes,       // top nodes
//       int max_elnd, int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd,
//       int** g_indexL, int* off_nd, int offM[], int aux[], int mem[], int Mat_q[], int len_q[],
//       int len_qoff[]);
  
    /// ===========================================================================================
    /// PROLUNG
    // ====================================================

  void compute_prol(
      int* n_nodes_lev,                       // nodes for level
      int* off_el,                            // element offset
      std::vector<std::pair<int, int>> v_el,  // element ordering
      int** elem_sto,                         // element structure
      int* v_inv_nd,                          // node inverse ordering
      int* v_inv_el,
      int** g_indexL,  // node map
      int* min_nd[],    // node offset
      int * max_nd[]
  ); 
  
  
  /// Computing Prolongation
  int compute_prol_qq(int name_label,
      int Level1, int Level_row, int Level_clmn, int n_nodes_c, int* off_el,
      std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd, int** g_indexL, int* min_nd, int* max_nd,
      double values_q[], int Prol_q[], int len_q[], int len_qoff[]);

//   int compute_prol_kk(
//       int Level1,       // level
//       int Level_row,    // row level
//       int Level_clmn,   // column level
//       int ndof_row_ql,  // element dofs
//       int n_nodes_f,
//       int* off_el,                            // offset elements vector
//       std::vector<std::pair<int, int>> v_el,  // element ordering
//       int** elem_sto,                         // element storage
//       int* v_inv_nd,                          // node ordering
//       int* v_inv_el,
//       int** g_indexL,     // map nodes
//       int* min_nd,        // offset nodes vector
//       int* max_nd,        // offset nodes vector
//       double values_q[],  // prolongation nonzero values
//       int Prol_q[],       // prolongation compressed pos
//       int len_q[],        // # of diagonal nonzero row entries
//       int len_qoff[]      // # of offdiagonal nonzero row entries
//   );
 /// ===========================================================================================
  
  // ====================================================================================
  void compute_rest(
      int* n_nodes_lev,                      // node level vector
      int max_elnd,                          // max # of elements at each point
      int* off_el,                           // offset element vectors
      std::vector<std::pair<int, int>> v_el  // element  ordering
      ,
      int** elem_sto,  // storage element structure
      int* v_inv_nd,   // node inverse ordering
      int* v_inv_el,
      int** g_indexL,  //  map nodes (level)
      int* min_nd[],    //  node offset vector
      int* max_nd[]    //  node offset vector
  );
  
  /// Computing Restriction
  int compute_res_qq(int iql,
      int Level1, int Level_row, int Level_clmn,  int n_nodes_f, int n_nodes_c, int max_elnd,
      int* off_el, std::vector<std::pair<int, int>> v_el, int** elem_sto, int* v_inv_nd, int** g_indexL,
      int* min_nd,int *max_nd, int mem[], double values_q[], int Prol_q[], int len_q[], int len_qoff[]);
  
  int compute_res_kk(
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
      int* v_inv_el,
      int** g_indexL,     // map nodes
      int* min_nd,        // offset nodes vector
      int* max_nd,        // offset nodes vector
      int mem[],          //  memory index
      double values_q[],  // restriction nonzero values
      int Rest_q[],       // restriction compressed pos
      int len_q[],        // # of diagonal nonzero row entries
      int len_qoff[]      // # of offdiagonal nonzero row entries
  );
  
  
  
  // PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
  // print ------------------------------------------------
  // PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
  void print_lib_mesh_h5(int** elem_sto, int n_elements, int dcl_nel);
  void print_lib_mesh(MGGeomEl& type_FEM, int** elem_sto, int n_elements, int dcl_nel);
 

  /// Print mesh -----------------------------------------------------------------
  void printMesh(
      BoundaryMesh& bd_mesht, Mesh& msht, BoundaryMesh& bd_msh0, Mesh& msh, const int type_FEM_in[]);

  /// Print mesh xdmf format
  void print_multimesh(const int* n_elements_lev, const int* n_nodes_lev);
  /// Print mesh (hdf5 format)
  void print_mesh_h5(
      int* n_nodes_lev, const int* map_mesh_in, int n_nodes, double* nod_val, int* min_nd[] , int* max_nd[],int n_elements,
      int* n_elements_lev, int bd_n_elements, int n_nodes_b, int* off_el, int* v_inv_nd, int* v_inv_el,
      int** elem_sto, int** elem_conn, std::vector<std::pair<int, int>> v_el, int** bd_elem_sto,
      std::vector<std::pair<int, int>> v_elb, int* bd_off_el, int** g_indexL,
      std::vector<std::pair<int, int>> v, int* nod_flag,
      int* mat_flag  // needed for bc and material with gambit
  );
  // print mesh
  void print_MedToMg(
      int Level, BoundaryMesh& bd_msh0,
      Mesh& msh0,  /// coarse mesh
      int n_groups, int* group_id_names, std::vector<std::pair<int, int>> v, int* v_inv_nd,
      std::vector<std::pair<int, int>> v_elb, int* nod_flag, int* mat_flag, int* off_el, int* v_inv_el,
      std::vector<std::pair<int, int>> v_el, std::vector<int> ElementsPerLevel, int** elem_sto,
      std::string filename = "default_coarse_mesh"  /// namefile
  );
  /// Print restriction operator (hdf5 format)---------------------------------------------------
  void print_op_h5(
      std::string name, int n_nodes_row, int n_nodes_cln, int count_q, int* Res_q, double* values_q,
      int* len_q, int* len_qoff, int qq);
  const char* cout;
  
  /// Print matrix operator (hdf5 format) --------------------------------------------------------
 void print_Mat(int* Mat, int nrow, int ncln);
  


 private:
  /// Print internal node structure (hdf5 format)
  void print_lib_node_h5(int** nod_sto, int n_nodes, int n_ndv);
  /// Print internal node structure (text format)
  void print_lib_node(int** nod_sto, int n_nodes, int n_ndv);
  ///  print med
  void print_med(
      int Level, BoundaryMesh& bd_msh0,
      Mesh& msh0,  /// coarse mesh
      int n_groups, int* group_id_names, std::vector<std::pair<int, int>> v,
      std::vector<std::pair<int, int>> v_elb, int* nod_flag, int* mat_flag,
      std::string filename = "default_coarse_mesh"  /// namefile
  );
};

#endif
