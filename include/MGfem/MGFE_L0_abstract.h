#pragma once

// #include "Domain_conf.h"

// #include "Printinfo_conf.h"

#include <fstream>
#include <vector>

#include "MGUtils_L0.h"
#include "MGquadrature_L0.h"
/// Class containing mathematical information about the finite element.
/**
 *
 * */


// ======================================================================================
//  FEM 
// ======================================================================================
// Note on GEOM dimension
// DIMENSION = max DIM geometry (from mesh from SALOME)
// MDIM_FEM  = max dim FEM (volume(3),surface(2),line(1)
// _dim = dimension FEM  (for MDIM_FEM=3; line=MDIM_FEM-2,surface=MDIM_FEM-1,volume=MDIM_FEM-1)
// dim paraview =3 always


// 
//     /// \f$ \xi,\eta,\chi \f$ coordinates of the HEX27 nodes
//     const int _CooH27[27 * 3] = {
//         // 0  1   2    3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
//         // 26
//         -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
//             1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//             -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//             1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//             -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//             0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//         };
//     /// \f$ \xi,\eta \f$ coordinates of the QUAD9 nodes
//     // note _CooQ9[i]=[_CooH27[i+12],_CooH27[27+i+12]
     // const int _CooQ9[9 * 2] = {
     //     /**< 1   2   3   4   5   6   7   8   9  node */
     //     -1, 1,  1, -1, 0,  1, 0, -1, 0, /**< coordinate xi */
     //         -1, -1, 1, 1,  -1, 0, 1, 0,  0  /**< coordinate eta */
     //     };
//     /// \f$ \xi \f$ coordinates of the EDGE3 nodes
//     // Note _CooE3[i]=[_CooH27[i+6]]
//     const int _CooE3[3] = {
//         // 1   2   3
//         -1, 1, 0  /// xi
//         };
// 
// 
//     const int _CooTetra10[40] = {
//         // L_i = l(j)*(a* l(k) + b)
//         //   j   k   a    b
//         0, 0, 2, -1,  // 0
//         1, 1, 2, -1,  // 1
//         2, 2, 2, -1,  // 2
//         3, 3, 2, -1,  // 3
//         1, 0, 4, 0,   // 4
//         1, 2, 4, 0,   // 5
//         2, 0, 4, 0,   // 6
//         0, 3, 4, 0,   // 7
//         1, 3, 4, 0,   // 8
//         2, 3, 4, 0    // 9
//     };
// 
//     const int _CooTetra4[16] = {
//         // L_i = a_i*xi + b_i*eta + c_i*chi + d_i
//         //   a   b   c   d
//         -1, -1, -1, 1,  // 0
//             1,  0,  0,  0,  // 1
//             0,  1,  0,  0,  // 2
//             0,  0,  1,  0   // 3
//         };
// 
//     const int _CooTriEl[9] = {
//         //  note _CooTriEl[i]=_CooTetra4
//         // L_i = a_i*xi + b_i*eta + c_i
//         //   a   b   c
//         -1, -1, 1,  // 0
//             1,  0,  0,  // 1
//             0,  1,  0   // 2
//         };
// 
// 
// 
// 
//     /// Offset for reading the generic node coordinates inside the _CooH27 array
     const int _H27Off = 27;
//     /// Offset for reading the generic node coordinates inside the _CooQ9 array
     const int _Q9Off = 9;
// 
// 
// 
// 
// };



// C++ includes ----------------------->
#include <cstddef>
#include <vector>
#include <memory>

// Femus includes ----------------------->
#include "MGFE_L0_II_enum_order.h"
// #include "MGquadrature_L0.h"


/**
 * \enum libMesh::FEContinuity defines an \p enum for finite element
 * types to libmesh_assert a certain level (or type? Hcurl?) of continuity.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum FEContinuity : int;
 * reducing header file dependencies.
 */
// enum FEContinuity : int {  DISCONTINUOUS,   C_ZERO,   C_ONE,   H_CURL,   H_DIV,SIDE_DISCONTINUOUS};

/**
 * \enum libMesh::FEFieldType defines an \p enum for finite element
 * field types - i.e. is it a scalar element, vector, tensor, etc.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum FEFieldType : int;
 * reducing header file dependencies.
 */
enum FEFieldType : int {  TYPE_SCALAR = 0,  TYPE_VECTOR};




#define TOLERANCE (1.e-6)
               
    
               
/**
 * class FEType hides (possibly multiple) FEFamily and approximation
 * orders, thereby enabling specialized finite element families.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Manages the family, order, etc. parameters for a given FE.
 */



/**
 * This class forms the foundation from which generic finite elements
 * may be derived.  In the current implementation, the templated
 * derived class \p FE offers a wide variety of commonly used finite
 * element concepts.  Check there for details.  Use the \p
 * FE0::build() method to create an object of any of the
 * derived classes.
 *
 * \note In the present design, the number of virtual members is kept
 * to a minimum for performance reasons, although this is not based on
 * rigorous profiling.
 *
 * All calls to static members of the \p FE classes should be
 * requested through the \p FEInterface.  This interface class
 * approximates runtime polymorphism for the templated finite element
 * classes.  Even internal library classes, like \p DofMap, request
 * the number of DOFs through this interface class.  This approach
 * also enables the co-existence of various element-based schemes.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * 
 */
// class FEMap; // libmesh class name

class MGFE0{

protected:

  //data ++++++++++++++++++++++++++++++++
  // std::unique_ptr<FEMap> _fe_map;
  const unsigned int _dim;/* The dimensionality of the object  */
  unsigned int _n_shape;

  //  type of computations for FE
  mutable  bool _calculations_started;// Have calculations with this object already been started?
  mutable  bool _calculate_nothing;/* Are we potentially deliberately calculating nothing? */
  mutable  bool _calculate_map;/* Are we calculating mapping functions? */
  mutable  bool _calculate_phi;/* Should we calculate shape functions?  */
  mutable  bool _calculate_dphi;/* Should we calculate shape function gradients?  */
  mutable  bool _calculate_d2phi;  /* Should we calculate shape function hessians? */
  mutable  bool _calculate_curl_phi;  /* Should we calculate shape function curls?  */
  mutable  bool _calculate_div_phi;/* Should we calculate shape function divergences? */
  mutable  bool _calculate_dphiref; /* Should we calculate reference shape function gradients?  */

  FEType _fe_type; // class FEType(const int   o = 1,   const FEFamily f = LAGRANGE) :    order(o),    family(f) 
  ElemType _elem_type; /* enum The element type the current data structures are set up for. */
  unsigned int _elem_p_level;//  element p-refinement level
  unsigned int _p_level; /* The p refinement level the current data structures are set up for.  */

  //Gauss point ---------------------------------------
  QBase * _qrule; /* A pointer to the quadrature rule employed */
  bool _shapes_on_quadrature;
  unsigned int _n_total_qp;/* The total number of quadrature points for the current configuration  */
  
  // reinit ------------------------------------->
  virtual bool shapes_need_reinit(){return false;} // needs to be re-initialized for each new element
  bool _add_p_level_in_reinit; /* Whether to add p-refinement levels in init/reinit methods  */
  
  // Constructor-Destructor -----------------------------------------------
protected:  
  /* Constructor. Protected so that this base classcannot be explicitly instantiated. */
  MGFE0 (const unsigned int dim, const FEType & fet);
public:
  virtual ~MGFE0(); /* Destructor. */
  /* Builds a specific finite element type. \returns A std::unique_ptr<FE0> to the FE object to prevent memory leaks.  */
  static std::unique_ptr<MGFE0> build (const unsigned int dim, const FEType & type);

  virtual void reinit (
      const int  elem, 
      const std::vector<double[3]> * const pts = nullptr,  
      const std::vector<double> * const weights = nullptr) {};

  virtual void reinit (
      const int elem,
      const unsigned int side,
      const double tolerance = TOLERANCE,
      const std::vector<double[3]> * const pts = nullptr,
      const std::vector<double> * const weights = nullptr) {};

  virtual void edge_reinit (const int elem,
                            const unsigned int edge,
                            const double tolerance = TOLERANCE,
                            const std::vector<double[3]> * pts = nullptr,
                            const std::vector<double> * weights = nullptr) {};

  /* Computes the reference space quadrature points on the side of an element based on the side quadrature points.  */
  virtual void side_map (const int  elem,
                         const int side,
                         const unsigned int s,
                         const std::vector<double[3]> & reference_side_points,
                         std::vector<double[3]> &       reference_points) {};

  /* the point p is located on the reference element*/
 static bool on_reference_element(const double * p, const ElemType t,
                                   const double eps = TOLERANCE);
 
 
 // ======================================================================
 // Reference canonical geometry #of nodes, nodes,
  /* \returns The reference space coordinates of \p nodes based on the element type.  */
 ElemType get_ref_elemtype(const ElemType elemType);  // reteru the reference mesh geometry
 ElemType get_ref_face_elemtype(const ElemType elemType);  // reteru the reference mesh geometry
 void get_refspace_nodes(const ElemType t,double nodes[]);// nodes[(x,y,z)_0,(x,y,z)_1.... ] 3D
 int get_refspace_n_nodes(const ElemType t);// number of nodes on element  ElemType (1= from Geom 0= from elem) 



// get -----------------------------------------------------
  /* \returns the dimension of this FE  */
  unsigned int get_dim() const{ return _dim; }
  void get_nothing() const {_calculate_nothing=true; }
  // /* \returns The \p xyz spatial locations of the quadrature points on the element.  
  virtual  const std::vector<double[3]> & get_xyz() const  { _calculate_map = false; std::abort(); }
  // /* \returns The element Jacobian times the quadrature weight for each quadrature point.
  virtual  const std::vector<double> & get_JxW() const  { _calculate_map = false;  std::abort(); }
  
  

  /* \returns The element tangents in xi-direction at the quadrature points.  */
  virtual const std::vector<double [3]> & get_dxyzdxi(int  xi_direction) const
  { _calculate_map = false;  std::abort();}// !!!!!!!!!!!!! this->_fe_map->get_dxyzdxi(); }
  /* \returns The second partial derivatives in eta(i)-zeta(j).  */
  virtual  const std::vector<double [3]> & get_d2xyz(int eta_i, int zeta_j) const
  { _calculate_map = true;std::abort();}// !!!!!!!!!!!!!!this->_fe_map->get_d2xyzdetadzeta(); }
  /* \returns The dhat(x)_i/dx_j entry in that_xe transformation matrix from physical to local coordinates. */
  virtual  const std::vector<double> & get_dzetadz(int hatx_i,int x_j) const
  {_calculate_map = true; std::abort();}// !!!!!!!!!!!!!!this->_fe_map->get_dzetadz(); }

  /* \returns The tangent vectors for face integration.   */
  virtual  const std::vector<std::vector<double[3]>> & get_tangents() const
  {_calculate_map = false; std::abort();}// !!!!!!!!!!!!!!this->_fe_map->get_tangents(); }

  /* \returns The outward pointing normal vectors for face integration.  */
  virtual  const std::vector<double> & get_normals() const
  {_calculate_map = false; std::abort();}// !!!!!!!!!!!!!!this->_fe_map->get_normals(); }

  /* \returns The curvatures for use in face integration.   */
  virtual  const std::vector<double> & get_curvatures() const
  {_calculate_map = false; std::abort();}//this->_fe_map->get_curvatures();}

 // ---------------------------------------------------------->
 /* Provides the class with the quadrature rule.  Implement this in derived classes.*/
  virtual void attach_quadrature_rule (QBase* q)=0;
  /* \returns The total number of approximation shape functions  */
  virtual unsigned int n_shape_functions () const=0;
  /* \returns The total number of quadrature points */
  virtual unsigned int n_quadrature_points () const;

  /* \returns The element type that the current shape functions.  */
  ElemType get_type()  const { return _elem_type; }
  /* \returns The p refinement level that the current shape functions have been calculated for.  */
  unsigned int get_p_level() const { return _p_level; }
  /* \returns The FE Type (approximation order and family) of the finite element.  */
  FEType get_fe_type()  const { return _fe_type; }
  /*   * \returns The approximation order of the finite element.  */
  int get_order() const  { return _fe_type.order + _p_level; }
   /* Sets the *base* FE order of the finite element.  */
  void set_fe_order(Order new_order) { _fe_type.order = new_order; }
  
  /*\returns The continuity level of the finite element.  */
  virtual FEContinuity get_continuity() const {return FEContinuity::C_ZERO;};
  /* \returns \p true if the finite element's higher order shape functions are hierarchic   */
  virtual bool is_hierarchic() const {return false;};//
  /* \returns The finite element family of this element.  */
  FEFamily get_family()  const { return _fe_type.family; }
  

  /* \returns The mapping object   * \note for InfFE, this gives a useless object.  */
  // const int * & get_fe_map() const { return nullptr;}// *_fe_map.get(); }
  // int & get_fe_map() { return nullptr;}//  *_fe_map.get(); }

  
  //Print ---------------------------------------------
 
  void print_JxW(std::ostream & os) const{}; /* _fe_map Prints the Jacobian times the weight for each quadrature point.  */
  virtual void print_phi(std::ostream & os) const {};/* _fe_map Prints the value of each shape function at each quadrature point.*/
  virtual void print_dphi(std::ostream & os) const{};// _fe_map Prints the value of each shape function's derivative
  virtual void print_d2phi(std::ostream & os) const {};// _fe_map Prints the value of each shape function's second derivatives
  void print_xyz(std::ostream & os) const{}; // _fe_map Prints the spatial location of each quadrature point
  void print_info(std::ostream & os) const{}; /* _fe_mapPrints all the relevant information about the current element.  */
  /* Same as above, but allows you to print to a stream. */
  
  // request  calculations --------------------->
  // friend std::ostream & operator << (std::ostream & os, const MGFE0 & fe){};_fe_map
  virtual void request_phi() const = 0;///* request phi calculations  */_fe_map
  virtual void request_dphi() const= 0;///* request dphi calculations  */_fe_map
  

// reinit ------------------------------------->
  /* Indicate whether to add p-refinement levels in init/reinit methods  */
  void add_p_level_in_reinit(bool value) { _add_p_level_in_reinit = value; }
  /*Whether to add p-refinement levels in init/reinit methods  */
  bool add_p_level_in_reinit() const { return _add_p_level_in_reinit; }

protected:

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEMap::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected. This needs to be implemented in the
   * derived class since this function depends on whether
   * the shape functions are vector-valued or not.
   */
  virtual void compute_shape_functions(const int /*Elem */* elem, const double & qp) {};

};





