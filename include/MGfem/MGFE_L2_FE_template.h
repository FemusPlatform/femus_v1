
 // #ifndef LIBMESH_FE_H
 // #define LIBMESH_FE_H

#pragma once 
 
// // Femus Local includes ----------------->
#include "MGFE_L0_II_enum_order.h"
#include "MGFE_L1_base.h"




// // C++ includes
// #include <cstddef>
// 
// namespace libMesh
// {
// 
// // forward declarations
// class DofConstraints;
// class DofMap;
// class QGauss;
// 
// #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
// 
// template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
// class InfFE;
// 
// #endif
// 
// 
// /**
//  * Most finite element types in libMesh are scalar-valued
//  */
// template <FEFamily T>
// struct FEOutputType
// {
//   typedef Real type;
// };
// 
// 
// /**
//  * Specialize for non-scalar-valued elements
//  */
// template<>
// struct FEOutputType<LAGRANGE_VEC>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<L2_LAGRANGE_VEC>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<HIERARCHIC_VEC>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<L2_HIERARCHIC_VEC>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<NEDELEC_ONE>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<MONOMIAL_VEC>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<RAVIART_THOMAS>
// {
//   typedef RealVectorValue type;
// };
// 
// template<>
// struct FEOutputType<L2_RAVIART_THOMAS>
// {
//   typedef RealVectorValue type;
// };
// 
// 
// /**
//  * A specific instantiation of the \p FEBase class. This
//  * class is templated, and specific template instantiations
//  * will result in different Finite Element families. Full specialization
//  * of the template for specific dimensions(\p Dim) and families
//  * (\p T) provide support for specific finite element types.
//  * The use of templates allows for compile-time optimization,
//  * however it requires that the specific finite element family
//  * and dimension is also known at compile time.  If this is
//  * too restricting for your application you can use the
//  * \p FEBase::build() member to create abstract (but still optimized)
//  * finite elements.
//  *
//  * \author Benjamin S. Kirk
//  * \date 2002-2007
//  * \brief Template class which generates the different FE families and orders.
//  */
 template <unsigned int Dim, FEFamily T>
 class MGFE : public MGFEbase{
 public:
   // // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  // static constexpr int MDIM_FEM=3; // max geom dim FEM (volume(3))
   //  //   static constexpr int NDOF_FEM=27; // max dof FEM (volume(3) LAGRANGE quad=27)
   //  //     static constexpr int NDOF_FEMB=9; // max dof FEM (volume(2) LAGRANGE quad=9)
   //  //Input parameter ----------------------------------------------------------------------------
   //  std::string femfileh5_;
     // int _dim;  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
     // int _order;///< Order of the shape functions 0=cpnst 1=linear 2=quad
     // int _NoShape;  ///< Number of shape functions of the fem
     // int  _NoElNodes;  ///< Number of nodes in the reference element (always ref quad)
     // int _GeomType;   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
     // int _FamType;    ///< Family  (geom 1: Lagrangian =0  RT =1)
   // 
   //  // Gaussian points -----------------------------------------------------------------------------
   //  int _NoGauss;  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
   //  double*  _xg;    ///< gaussian points
   //  double* _weight;              ///< Weight
   //  double* _phi_map;             ///< Shape functions
   //  double* _dphidxez_map;        ///< Shape derivative functions in gaussian points
   //  double* _dphidxx_map;         ///< Second order shape derivatives in Gaussian points
   // 
   //  // Dofs points ----------------------------------------------------------------------------------
   //  double *_xnodes;              ///< xnode (hex or tetra) in (xi,eta,seta)
   //  int *_fem_conn;              ///< dof node coords (respect to _xnodes)
   //  int     _deg;                  ///< Polinomial degree
   // 
   //  // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
   //  double* _phi_map_nodes;             ///< Shape functions
   //  double* _dphidxez_map_nodes;  ///< Shape derivative functions in nodal points
   //  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
   //Input parameter ----------------------------------------------------------------------------
  const int get_dim()      ;// {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
  const int get_order()    ;// {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
  const int get_NoShape()  ;// {return _NoShape;}  ///< Number of shape functions of the fem
  const int get_NoElNodes() ;//{return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
  const int get_GeomType();//  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
  const int get_FamType() ;//  {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
  const int get_NoGauss();//{ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  const double  get_xg(int i);//{          return  *_xg;}          ///< gaussian points
  const double  get_weight(int i);//{      return *_weight;}            ///< Weight
  const double &  get_phi();//{     return  *_phi_map;}           ///< Shape functions
  const double & get_dphidxez_map();//{return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  const double & get_dphidxx_map();//{ return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
  
   // Dofs points ----------------------------------------------------------------------------------
  const double &  get_xnodes();//{return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
  const int &  get_fem_conn();//{ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  const int     get_deg();//{return _deg};                  ///< Polinomial degree
  
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
  const double & get_phi_map_nodes();//{return   *_phi_map_nodes;   }          ///< Shape functions
  const double & get_dphidxez_map_nodes();//{return   *_dphidxez_map_nodes;}  ///< Shape derivative functions in nodal points   
    
  double compute_fe_gauss_points(const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
                              double InvJac[]    
 );
   double compute_fe_face_gauss_points(const int  gp, // gaussian point gp
   const double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
   
 );  
    
   virtual void set_par(std::string fileh5,int dim_in,int order_in,int shapes_quad,int geom_type_in,int fam_type_in)
      {_dim=dim_in, _order= order_in, _NoShape=shapes_quad, _GeomType=geom_type_in,   _FamType=fam_type_in;}
   // virtual int get_refspace_n_nodes(const ElemType itemType, const int o);
// 
   /* Constructor.  */
   explicit   MGFE(const FEType & fet);
// 
//   typedef typename
//   FEGenericBase<typename FEOutputType<T>::type>::OutputShape
//   OutputShape;
// 
  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the dimension,
   * element type, and order directly.  This allows the method to
   * be static.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static double  shape(
    const ElemType t,
    const Order o,
    const unsigned int i,
    const double  p[]);
// 
//   /**
//    * \returns The value of the \f$ i^{th} \f$ shape function at
//    * point \p p.  This method allows you to specify the dimension,
//    * element type, and order directly.  This allows the method to
//    * be static.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static OutputShape shape(const Elem * elem,
//                            const Order o,
//                            const unsigned int i,
//                            const Point & p,
//                            const bool add_p_level = true);
// 
//   /**
//    * \returns The value of the \f$ i^{th} \f$ shape function at
//    * point \p p.  This method allows you to specify the dimension and
//    * element type directly. The order is given by the FEType.
//    * This allows the method to be static.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//    static OutputShape shape(const FEType fet,
//                             const Elem * elem,
//                             const unsigned int i,
//                             const Point & p,
//                             const bool add_p_level = true);
// 
//   /**
//    * Fills \p v with the values of the \f$ i^{th} \f$
//    * shape function, evaluated at all points p.  You must specify
//    * element order directly.  \p v should already be the appropriate
//    * size.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static void shapes(const Elem * elem,
//                      const Order o,
//                      const unsigned int i,
//                      const std::vector<Point> & p,
//                      std::vector<OutputShape> & v,
//                      const bool add_p_level = true);
// 
// 
//   /**
//    * Fills \p v[i][qp] with the values of the \f$ i^{th} \f$
//    * shape functions, evaluated at all points in p.  You must specify
//    * element order directly.  \p v should already be the appropriate
//    * size.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static void all_shapes(const Elem * elem,
//                          const Order o,
//                          const std::vector<Point> & p,
//                          std::vector<std::vector<OutputShape> > & v,
//                          const bool add_p_level = true);
// 
//   /**
//    * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function at point \p p.  This method allows you to
//    * specify the dimension, element type, and order directly.
//    *
//    * On a p-refined element, \p o should be the total order of the element.
//    */
   static double shape_deriv(const ElemType t,  const Order o,  const unsigned int i, const unsigned int j, const double  p[]);
   
   
   
// 
//   /**
//    * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function.  You must specify element type, and order directly.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static OutputShape shape_deriv(const Elem * elem,
//                                  const Order o,
//                                  const unsigned int i,
//                                  const unsigned int j,
//                                  const Point & p,
//                                  const bool add_p_level = true);
// 
//   /**
//    * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function.  You must specify element type, and order (via
//    * FEType) directly.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static OutputShape shape_deriv(const FEType fet,
//                                  const Elem * elem,
//                                  const unsigned int i,
//                                  const unsigned int j,
//                                  const Point & p,
//                                  const bool add_p_level = true);
// 
//   /**
//    * Fills \p v with the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function, evaluated at all points p.  You must specify
//    * element order directly.  \p v should already be the appropriate
//    * size.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static void shape_derivs(const Elem * elem,
//                            const Order o,
//                            const unsigned int i,
//                            const unsigned int j,
//                            const std::vector<Point> & p,
//                            std::vector<OutputShape> & v,
//                            const bool add_p_level = true);
// 
//   /**
//    * Fills \p comps with dphidxi (and in higher dimensions, eta/zeta)
//    * derivative component values for all shape functions, evaluated at
//    * all points in p.  You must specify element order directly.
//    * Output component arrays in \p comps should already be the
//    * appropriate size.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static void all_shape_derivs(const Elem * elem,
//                                const Order o,
//                                const std::vector<Point> & p,
//                                std::vector<std::vector<OutputShape>> * comps[3],
//                                const bool add_p_level = true);
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
//   /**
//    * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function at the point \p p.
//    *
//    * \note Cross-derivatives are indexed according to:
//    * j = 0 ==> d^2 phi / dxi^2
//    * j = 1 ==> d^2 phi / dxi deta
//    * j = 2 ==> d^2 phi / deta^2
//    * j = 3 ==> d^2 phi / dxi dzeta
//    * j = 4 ==> d^2 phi / deta dzeta
//    * j = 5 ==> d^2 phi / dzeta^2
//    *
//    * \note Computing second derivatives is not currently supported for
//    * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
//    * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
//    * All other element types return an error when asked for second
//    * derivatives.
//    *
//    * On a p-refined element, \p o should be the total order of the element.
//    */
   static  double shape_second_deriv(
     const ElemType t,
     const Order o,
     const unsigned int i,
     const unsigned int j,
     const double  p[]);
// 
//   /**
//    * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function at the point \p p.
//    *
//    * \note Cross-derivatives are indexed according to:
//    * j = 0 ==> d^2 phi / dxi^2
//    * j = 1 ==> d^2 phi / dxi deta
//    * j = 2 ==> d^2 phi / deta^2
//    * j = 3 ==> d^2 phi / dxi dzeta
//    * j = 4 ==> d^2 phi / deta dzeta
//    * j = 5 ==> d^2 phi / dzeta^2
//    *
//    * \note Computing second derivatives is not currently supported for
//    * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
//    * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
//    * All other element types return an error when asked for second
//    * derivatives.
//    *
//    * On a p-refined element, \p o should be the base order of the
//    * element if \p add_p_level is left \p true, or can be the base
//    * order of the element if \p add_p_level is set to \p false.
//    */
//   static OutputShape shape_second_deriv(const Elem * elem,
//                                         const Order o,
//                                         const unsigned int i,
//                                         const unsigned int j,
//                                         const Point & p,
//                                         const bool add_p_level = true);
// 
//   /**
//    * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
//    * shape function at the point \p p.
//    *
//    * \note Cross-derivatives are indexed according to:
//    * j = 0 ==> d^2 phi / dxi^2
//    * j = 1 ==> d^2 phi / dxi deta
//    * j = 2 ==> d^2 phi / deta^2
//    * j = 3 ==> d^2 phi / dxi dzeta
//    * j = 4 ==> d^2 phi / deta dzeta
//    * j = 5 ==> d^2 phi / dzeta^2
//    *
//    * \note Computing second derivatives is not currently supported for
//    * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
//    * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
//    * All other element types return an error when asked for second
//    * derivatives.
//    *
//    * On a p-refined element, \p o should be the total order of the element.
//    */
//   static OutputShape shape_second_deriv(const FEType fet,
//                                         const Elem * elem,
//                                         const unsigned int i,
//                                         const unsigned int j,
//                                         const Point & p,
//                                         const bool add_p_level = true);
// 
// #endif //LIBMESH_ENABLE_SECOND_DERIVATIVES
//   /**
//    * Build the nodal soln from the element soln.
//    * This is the solution that will be plotted.
//    *
//    * On a p-refined element, \p o should be the base order of the element.
//    */
//   static void nodal_soln(const Elem * elem, const Order o,
//                          const std::vector<Number> & elem_soln,
//                          std::vector<Number> & nodal_soln,
//                          bool add_p_level = true);
// 
//   /**
//    * Build the nodal soln on one side from the (full) element soln.
//    * This is the solution that will be plotted on side-elements.
//    *
//    * On a p-refined element, \p o should be the base order of the element.
//    */
//   static void side_nodal_soln(const Elem * elem, const Order o,
//                               const unsigned int side,
//                               const std::vector<Number> & elem_soln,
//                               std::vector<Number> & nodal_soln_on_side,
//                               bool add_p_level = true);
// 
   /* \returns The number of shape functions associated with this finite element.  */
 virtual unsigned int n_shape_functions () const override;
// 
//   /* \returns The number of shape functions associated with a finite element of type \p t and approximation order \p o.
//    * On a p-refined element, \p o should be the total order of the element.
  static unsigned int n_shape_functions (const ElemType t, const Order o){ return MGFE<Dim,T>::n_dofs (t,o); }

  /**
   * \returns The number of shape functions associated with this
   * finite element.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs(const ElemType t,const Order o);

//   /**
//    * \returns The number of dofs at node \p n for a finite element
//    * of type \p t and order \p o.
//    *
//    * On a p-refined element, \p o should be the total order of the element.
//    */
 static unsigned int n_dofs_at_node(const ElemType t,const Order o, const unsigned int n);
// 
//   /**
//    * \returns The number of dofs interior to the element,
//    * not associated with any interior nodes.
//    *
//    * On a p-refined element, \p o should be the total order of the element.
//    */
   static unsigned int n_dofs_per_elem(const ElemType t, const Order o);

   /**
    * \returns The continuity level of the finite element.
    */
   FEContinuity get_continuity() const override;
// 
//   /**
//    * \returns \p true if the finite element's higher order shape functions are
//    * hierarchic
//    */
   virtual bool is_hierarchic() const override;
// 
//   /**
//    * Fills the vector di with the local degree of freedom indices
//    * associated with side \p s of element \p elem
//    *
//    * On a p-refined element, \p o should be the base order of the element.
//    */
//   static void dofs_on_side(const Elem * const elem,
//                            const Order o,
//                            unsigned int s,
//                            std::vector<unsigned int> & di,
//                            bool add_p_level=true);
//   /**
//    * Fills the vector di with the local degree of freedom indices
//    * associated with edge \p e of element \p elem
//    *
//    * On a p-refined element, \p o should be the base order of the element.
//    */
//   static void dofs_on_edge(const Elem * const elem,
//                            const Order o,
//                            unsigned int e,
//                            std::vector<unsigned int> & di,
//                            bool add_p_level=true);
// 
//   static Point inverse_map (const Elem * elem,
//                             const Point & p,
//                             const Real tolerance = TOLERANCE,
//                             const bool secure = true)
//   {
//     // libmesh_deprecated(); // soon
//     return FEMap::inverse_map(Dim, elem, p, tolerance, secure, secure);
//   }
// 
//   static void inverse_map (const Elem * elem,
//                            const std::vector<Point> & physical_points,
//                            std::vector<Point> &       reference_points,
//                            const Real tolerance = TOLERANCE,
//                            const bool secure = true)
//   {
//     // libmesh_deprecated(); // soon
//     FEMap::inverse_map(Dim, elem, physical_points, reference_points,
//                        tolerance, secure, secure);
//   }
// 
//   /**
//    * This is at the core of this class. Use this for each
//    * new element in the mesh.  Reinitializes all the physical
//    * element-dependent data based on the current element
//    * \p elem.  By default the shape functions and associated
//    * data are computed at the quadrature points specified
//    * by the quadrature rule \p qrule, but may be any points
//    * specified on the reference element specified in the optional
//    * argument \p pts.
//    */
//   virtual void reinit (const Elem * elem,
//                        const std::vector<Point> * const pts = nullptr,
//                        const std::vector<Real> * const weights = nullptr) override;
// 
//   /**
//   * This re-computes the dual shape function coefficients.
//   * The dual shape coefficients are utilized when calculating dual shape functions.
//   */
//   virtual void reinit_dual_shape_coeffs (const Elem * elem,
//                                          const std::vector<Point> & pts,
//                                          const std::vector<Real> & JxW) override;
// 
//   /**
//    * This computes the default dual shape function coefficients.
//    * The dual shape coefficients are utilized when calculating dual shape functions.
//    */
//    virtual void reinit_default_dual_shape_coeffs (const Elem * elem) override;
// 
//   /**
//    * Reinitializes all the physical element-dependent data based on
//    * the \p side of \p face.  The \p tolerance parameter is passed to
//    * the involved call to \p inverse_map().  By default the shape
//    * functions and associated data are computed at the quadrature
//    * points specified by the quadrature rule \p qrule, but may be any
//    * points specified on the reference \em side element specified in
//    * the optional argument \p pts.
//    */
//   virtual void reinit (const Elem * elem,
//                        const unsigned int side,
//                        const Real tolerance = TOLERANCE,
//                        const std::vector<Point> * const pts = nullptr,
//                        const std::vector<Real> * const weights = nullptr) override;
// 
//   /**
//    * Reinitializes all the physical element-dependent data based on
//    * the \p edge.  The \p tolerance parameter is passed to the
//    * involved call to \p inverse_map().  By default the shape
//    * functions and associated data are computed at the quadrature
//    * points specified by the quadrature rule \p qrule, but may be any
//    * points specified on the reference \em side element specified in
//    * the optional argument \p pts.
//    */
//   virtual void edge_reinit (const Elem * elem,
//                             const unsigned int edge,
//                             const Real tolerance = TOLERANCE,
//                             const std::vector<Point> * const pts = nullptr,
//                             const std::vector<Real> * const weights = nullptr) override;
// 
//   /**
//    * Computes the reference space quadrature points on the side of
//    * an element based on the side quadrature points.
//    */
//   virtual void side_map (const Elem * elem,
//                          const Elem * side,
//                          const unsigned int s,
//                          const std::vector<Point> & reference_side_points,
//                          std::vector<Point> &       reference_points) override;
// 
//   /**
//    * Computes the reference space quadrature points on the side of
//    * an element based on the edge quadrature points.
//    */
//   virtual void edge_map (const Elem * elem,
//                          const Elem * edge,
//                          const unsigned int e,
//                          const std::vector<Point> & reference_edge_points,
//                          std::vector<Point> &       reference_points);
// 
  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  void attach_quadrature_rule (QBase *q) override;
// 
// #ifdef LIBMESH_ENABLE_AMR
//   /**
//    * Computes the constraint matrix contributions (for
//    * non-conforming adapted meshes) corresponding to
//    * variable number \p var_number, using element-specific
//    * optimizations if possible.
//    */
//   static void compute_constraints (DofConstraints & constraints,
//                                    DofMap & dof_map,
//                                    const unsigned int variable_number,
//                                    const Elem * elem);
// #endif // #ifdef LIBMESH_ENABLE_AMR
// 
//   /**
//    * \returns \p true when the shape functions (for
//    * this \p FEFamily) depend on the particular
//    * element, and therefore needs to be re-initialized
//    * for each new element.  \p false otherwise.
//    */
//   virtual bool shapes_need_reinit() const override;
// 
//   static Point map (const Elem * elem,
//                     const Point & reference_point)
//   {
//     // libmesh_deprecated(); // soon
//     return FEMap::map(Dim, elem, reference_point);
//   }
// 
//   static Point map_xi (const Elem * elem,
//                        const Point & reference_point)
//   {
//     // libmesh_deprecated(); // soon
//     return FEMap::map_deriv(Dim, elem, 0, reference_point);
//   }
// 
//   static Point map_eta (const Elem * elem,
//                         const Point & reference_point)
//   {
//     // libmesh_deprecated(); // soon
//     return FEMap::map_deriv(Dim, elem, 1, reference_point);
//   }
// 
//   static Point map_zeta (const Elem * elem,
//                          const Point & reference_point)
//   {
//     // libmesh_deprecated(); // soon
//     return FEMap::map_deriv(Dim, elem, 2, reference_point);
//   }
// 
// #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
//   /**
//    * make InfFE classes friends, so that these may access
//    * the private \p map, map_xyz methods
//    */
//   template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
//   friend class InfFE;
// #endif
// 
// protected:
// 
//   /**
//    * Update the various member data fields \p phi,
//    * \p dphidxi, \p dphideta, \p dphidzeta, etc.
//    * for the current element.  These data will be computed
//    * at the points \p qp, which are generally (but need not be)
//    * the quadrature points.
//    */
//   virtual void init_shape_functions(const std::vector<Point> & qp,
//                                     const Elem * e);
// 
//   /**
//    * A default implementation for all_shape_derivs
//    */
//   static void default_all_shape_derivs (const Elem * elem,
//                                         const Order o,
//                                         const std::vector<Point> & p,
//                                         std::vector<std::vector<OutputShape>> * comps[3],
//                                         const bool add_p_level = true);
// 
// 
//   /**
//    * Init \p dual_phi and potentially \p dual_dphi, \p dual_d2phi
//    */
//   void init_dual_shape_functions(unsigned int n_shapes, unsigned int n_qp);
// 

// 
//   /**
//    * A default implementation for shapes
//    */
//   static void default_shapes (const Elem * elem,
//                               const Order o,
//                               const unsigned int i,
//                               const std::vector<Point> & p,
//                               std::vector<OutputShape> & v,
//                               const bool add_p_level = true)
//     {
//       libmesh_assert_equal_to(p.size(), v.size());
//       for (auto vi : index_range(v))
//         v[vi] = FE<Dim,T>::shape (elem, o, i, p[vi], add_p_level);
//     }
// 
//   /**
//    * A default implementation for all_shapes
//    */
//   static void default_all_shapes (const Elem * elem,
//                                   const Order o,
//                                   const std::vector<Point> & p,
//                                   std::vector<std::vector<OutputShape>> & v,
//                                   const bool add_p_level = true)
//     {
//       for (auto i : index_range(v))
//         {
//           libmesh_assert_equal_to ( p.size(), v[i].size() );
//           FE<Dim,T>::shapes (elem, o, i, p, v[i], add_p_level);
//         }
//     }
// 
//   /**
//    * A default implementation for shape_derivs
//    */
//   static void default_shape_derivs (const Elem * elem,
//                                     const Order o,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const std::vector<Point> & p,
//                                     std::vector<OutputShape> & v,
//                                     const bool add_p_level = true)
//     {
//       libmesh_assert_equal_to(p.size(), v.size());
//       for (auto vi : index_range(v))
//         v[vi] = FE<Dim,T>::shape_deriv (elem, o, i, j, p[vi], add_p_level);
//     }
// 
//   /**
//    * A default implementation for side_nodal_soln
//    */
//   static void default_side_nodal_soln(const Elem * elem, const Order o,
//                                       const unsigned int side,
//                                       const std::vector<Number> & elem_soln,
//                                       std::vector<Number> & nodal_soln_on_side,
//                                       bool add_p_level = true);
// 
//   /**
//    * An array of the node locations on the last
//    * element we computed on
//    */
//   std::vector<Point> cached_nodes;
// 
//   /**
//    * The last side and last edge we did a reinit on
//    */
   ElemType last_side;
// 
   ElemType last_edge;
// };
// 
// 
// 
// /**
//  * Clough-Tocher finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Roy Stogner
//  * \date 2004
//  */
// template <unsigned int Dim>
// class FEClough : public FE<Dim,CLOUGH>
// {
// public:
// 
//   /**
//    * Constructor. Creates a hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEClough(const FEType & fet) :
//     FE<Dim,CLOUGH> (fet)
//   {}
// };
// 
// 
// 
// /**
//  * Hermite finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Roy Stogner
//  * \date 2005
//  */
// template <unsigned int Dim>
// class FEHermite : public FE<Dim,HERMITE>
// {
// public:
// 
//   /**
//    * Constructor. Creates a hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEHermite(const FEType & fet) :
//     FE<Dim,HERMITE> (fet)
//   {}
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
//   /**
//    * 1D hermite functions on unit interval
//    */
//   static Real hermite_raw_shape_second_deriv(const unsigned int basis_num,
//                                              const Real xi);
// #endif
//   static Real hermite_raw_shape_deriv(const unsigned int basis_num,
//                                       const Real xi);
//   static Real hermite_raw_shape(const unsigned int basis_num,
//                                 const Real xi);
// };
// 
// 
// 
// /**
//  * Subdivision finite elements.
//  *
//  * Template specialization prototypes are needed for calling from
//  * inside FESubdivision::init_shape_functions
//  */
// template <>
// Real FE<2,SUBDIVISION>::shape(const Elem * elem,
//                               const Order order,
//                               const unsigned int i,
//                               const Point & p,
//                               const bool add_p_level);
// 
// template <>
// Real FE<2,SUBDIVISION>::shape_deriv(const Elem * elem,
//                                     const Order order,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const Point & p,
//                                     const bool add_p_level);
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// template <>
// Real FE<2,SUBDIVISION>::shape_second_deriv(const Elem * elem,
//                                            const Order order,
//                                            const unsigned int i,
//                                            const unsigned int j,
//                                            const Point & p,
//                                            const bool add_p_level);
// 
// #endif
// 
// class FESubdivision : public FE<2,SUBDIVISION>
// {
// public:
// 
//   /**
//    * Constructor. Creates a subdivision surface finite element.
//    * Currently only supported for two-dimensional meshes in
//    * three-dimensional space.
//    */
//   FESubdivision(const FEType & fet);
// 
//   /**
//    * This is at the core of this class. Use this for each new
//    * non-ghosted element in the mesh.  Reinitializes all the physical
//    * element-dependent data based on the current element
//    * \p elem.  By default the shape functions and associated
//    * data are computed at the quadrature points specified
//    * by the quadrature rule \p qrule, but may be any points
//    * specified on the reference element specified in the optional
//    * argument \p pts.
//    */
//   virtual void reinit (const Elem * elem,
//                        const std::vector<Point> * const pts = nullptr,
//                        const std::vector<Real> * const weights = nullptr) override;
// 
//   /**
//    * This prevents some compilers being confused by partially
//    * overriding this virtual function.
//    */
//   virtual void reinit (const Elem *,
//                        const unsigned int,
//                        const Real = TOLERANCE,
//                        const std::vector<Point> * const = nullptr,
//                        const std::vector<Real> * const = nullptr) override
//   { libmesh_not_implemented(); }
// 
//   /**
//    * Provides the class with the quadrature rule, which provides the
//    * locations (on a reference element) where the shape functions are
//    * to be calculated.
//    */
//   virtual void attach_quadrature_rule (QBase * q) override;
// 
//   /**
//    * Update the various member data fields \p phi,
//    * \p dphidxi, \p dphideta, \p dphidzeta, etc.
//    * for the current element.  These data will be computed
//    * at the points \p qp, which are generally (but need not be)
//    * the quadrature points.
//    */
//   virtual void init_shape_functions(const std::vector<Point> & qp,
//                                     const Elem * elem) override;
// 
//   /**
//    * \returns The value of the \f$ i^{th} \f$ of the 12 quartic
//    * box splines interpolating a regular Loop subdivision
//    * element, evaluated at the barycentric coordinates \p v,
//    * \p w.
//    */
//   static Real regular_shape(const unsigned int i,
//                             const Real v,
//                             const Real w);
// 
//   /**
//    * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th}
//    * \f$ of the 12 quartic box splines interpolating a regular
//    * Loop subdivision element, evaluated at the barycentric
//    * coordinates \p v, \p w.
//    */
//   static Real regular_shape_deriv(const unsigned int i,
//                                   const unsigned int j,
//                                   const Real v,
//                                   const Real w);
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
//   /**
//    * \returns The second \f$ j^{th} \f$ derivative of the
//    * \f$ i^{th} \f$ of the 12 quartic box splines interpolating
//    * a regular Loop subdivision element, evaluated at the
//    * barycentric coordinates \p v, \p w.
//    */
//   static Real regular_shape_second_deriv(const unsigned int i,
//                                          const unsigned int j,
//                                          const Real v,
//                                          const Real w);
// 
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVE
//   /**
//    * Fills the vector \p weights with the weight coefficients
//    * of the Loop subdivision mask for evaluating the limit surface
//    * at a node explicitly. The size of \p weights will be
//    * 1 + \p valence, where \p valence is the number of neighbor
//    * nodes of the node where the limit surface is to be
//    * evaluated. The weight for the node itself is the first
//    * element of \p weights.
//    */
//   static void loop_subdivision_mask(std::vector<Real> & weights,
//                                     const unsigned int valence);
// 
// 
//   /**
//    * Builds the subdivision matrix \p A for the Loop scheme. The
//    * size depends on the element's \p valence.
//    */
//   static void init_subdivision_matrix(DenseMatrix<Real> & A,
//                                       unsigned int valence);
// };
// 
// 
// 
// /**
//  * Hierarchic finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Benjamin S. Kirk
//  * \date 2002-2007
//  */
// template <unsigned int Dim>
// class FEHierarchic : public FE<Dim,HIERARCHIC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEHierarchic(const FEType & fet) :
//     FE<Dim,HIERARCHIC> (fet)
//   {}
// };
// 
// 
// 
// /**
//  * Discontinuous Hierarchic finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Truman E. Ellis
//  * \date 2011
//  */
// template <unsigned int Dim>
// class FEL2Hierarchic : public FE<Dim,L2_HIERARCHIC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEL2Hierarchic(const FEType & fet) :
//     FE<Dim,L2_HIERARCHIC> (fet)
//   {}
// };
// 
// 
// 
// /**
//  * Lagrange finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Benjamin S. Kirk
//  * \date 2002-2007
//  */
// template <unsigned int Dim>
// class FELagrange : public FE<Dim,LAGRANGE>
// {
// public:
// 
//   /**
//    * Constructor. Creates a Lagrange finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FELagrange(const FEType & fet) :
//     FE<Dim,LAGRANGE> (fet)
//   {}
// };
// 
// 
// /**
//  * Discontinuous Lagrange finite elements.
//  */
// template <unsigned int Dim>
// class FEL2Lagrange : public FE<Dim,L2_LAGRANGE>
// {
// public:
// 
//   /**
//    * Constructor. Creates a discontinuous Lagrange finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEL2Lagrange(const FEType & fet) :
//     FE<Dim,L2_LAGRANGE> (fet)
//   {}
// };
// 
// 
// /**
//  * Monomial finite elements.  Still templated on the dimension,
//  * \p Dim.
//  *
//  * \author Benjamin S. Kirk
//  * \date 2002-2007
//  */
// template <unsigned int Dim>
// class FEMonomial : public FE<Dim,MONOMIAL>
// {
// public:
// 
//   /**
//    * Constructor. Creates a monomial finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEMonomial(const FEType & fet) :
//     FE<Dim,MONOMIAL> (fet)
//   {}
// };
// 
// 
// /**
//  * The FEScalar class is used for working with SCALAR variables.
//  */
// template <unsigned int Dim>
// class FEScalar : public FE<Dim,SCALAR>
// {
// public:
// 
//   /**
//    * Constructor. Creates a SCALAR finite element
//    * which simply represents one or more
//    * extra DOFs coupled to all other DOFs in
//    * the system.
//    */
//   explicit
//   FEScalar(const FEType & fet) :
//     FE<Dim,SCALAR> (fet)
//   {}
// };
// 
// 
// /**
//  * XYZ finite elements.  These require specialization
//  * because the shape functions are defined in terms of
//  * physical XYZ coordinates rather than local coordinates.
//  *
//  * \author Benjamin S. Kirk
//  * \date 2002-2007
//  */
// template <unsigned int Dim>
// class FEXYZ : public FE<Dim,XYZ>
// {
// public:
// 
//   /**
//    * Constructor. Creates a monomial finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEXYZ(const FEType & fet) :
//     FE<Dim,XYZ> (fet)
//   {}
// 
//   /**
//    * Explicitly call base class method.  This prevents some
//    * compilers being confused by partially overriding this virtual function.
//    * \note: pts need to be in reference space coordinates, not physical ones.
//    */
//   virtual void reinit (const Elem * elem,
//                        const std::vector<Point> * const pts = nullptr,
//                        const std::vector<Real> * const weights = nullptr) override
//   { FE<Dim,XYZ>::reinit (elem, pts, weights); }
// 
//   /**
//    * Reinitializes all the physical element-dependent data based on
//    * the \p side of \p face.
//    */
//   virtual void reinit (const Elem * elem,
//                        const unsigned int side,
//                        const Real tolerance = TOLERANCE,
//                        const std::vector<Point> * const pts = nullptr,
//                        const std::vector<Real> * const weights = nullptr) override;
// 
// 
// protected:
// 
//   /**
//    * Update the various member data fields \p phi,
//    * \p dphidxi, \p dphideta, \p dphidzeta, etc.
//    * for the current element.  These data will be computed
//    * at the points \p qp, which are generally (but need not be)
//    * the quadrature points.
//    */
//   virtual void init_shape_functions(const std::vector<Point> & qp,
//                                     const Elem * e) override;
// 
//   /**
//    * After having updated the jacobian and the transformation
//    * from local to global coordinates in \p FEAbstract::compute_map(),
//    * the first derivatives of the shape functions are
//    * transformed to global coordinates, giving \p dphi,
//    * \p dphidx, \p dphidy, and \p dphidz. This method
//    * should rarely be re-defined in derived classes, but
//    * still should be usable for children. Therefore, keep
//    * it protected.
//    */
//   virtual void compute_shape_functions(const Elem * elem, const std::vector<Point> & qp) override;
// 
//   /**
//    * Compute the map & shape functions for this face.
//    */
//   void compute_face_values (const Elem * elem,
//                             const Elem * side,
//                             const std::vector<Real> & weights);
   
      // ====================================================================================
    //  CONSTRUCTOR-DESTRUCTOR  INITIALIZING AND CLEANING
    // Constructor ------------------------------------------------------------------------
    // MGFE(MGFE(const FEType & fet)
    // const FEType  fet,
    // const int  dim_in,   ///< dimension
    // const int order_in,      ///< pol approx
    // const int shapes_in,     ///< Fem dof input
    // const int geom_type_in,  ///< Geom type: 1(Quad,Hex) or 0(Tri,Tet)
    // const int fam_type_in    ///< Fem type: 0(Lagrangian) or 1 (Raviart Thomas linear)
    // );
   
    // init functions ---------------------------------------------------------------------
//     void init_qua();  ///< Generates the Lagrangian quad shape functions
//     void init_lin();  ///< Generates the Lagrangian linear shape functions
//     void init_pie();  ///< Generates the Lagrangian piecewise shape functions
//     void init_pie_rec();
//     void init_pie_tri();
// 
//     void init_qua_tri();
// //   void init_qua_Gtri();
//     void init_qua_rec();
// //   void init_qua_Grec();
// 
//     void init_lin_tri();
// //   void init_lin_Gtri();
//     void init_lin_rec();
//     void init_lin_rec_rt();
//   void init_lin_Grec();
    /// Clear data substructures -------------------------------------------------------------
    void clear();
    ~MGFE();

    // get function
    int GetFamilyType(int elem_dof, int dim);  // the function returns 0 or 1
    // -------------------------------------------------------
    // Triangle - based elements:

    // double Tri_2d_LinearPhi(int nPhi, double point[]);
    // double Tri_2d_LinearDerPhi(int nPhi, double point[], int dir);
    // double Tri_2d_QuadraticPhi(int nPhi, double point[]);
    // double Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir);
    // double Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2);
    // 
    // double Tri_3d_LinearPhi(int nPhi, double point[]);
    // double Tri_3d_LinearDerPhi(int nPhi, double point[], int dir);
    // double Tri_3d_QuadraticPhi(int nPhi, double point[]);
    // double Tri_3d_QuadraticDerPhi(int nPhi, double point[], int dir);
    // double Tri_3d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2);
    // 
    // // -------------------------------------------------------
    // // Quadrangle - based elements:
    // double Edge_Lin_Phi(int PhiCoeff, double Coordinate);
    // double Edge_Lin_DPhi(int PhiCoeff, double Coordinate);
    // double Edge_Quad_Phi(int PhiCoeff, double Coordinate);
    // double Edge_Quad_DPhi(int PhiCoeff, double Coordinate);
    // double Edge_Quad_D2Phi(int PhiCoeff, double Coordinate);
    // 
    // double Rec_Lin_Phi(int nPhi, double point[], int dimension);
    // double Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer);
    // double Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2);
    // double Rec_Quad_Phi(int nPhi, double point[], int dimension);
    // double Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer);
    // double Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2);
    // 
    // double FirstDerivateOfLocalPhi(int nPhi, double point[], int dimension, int DirDer, int FamilyType);

// // Raviart Thomas function
//     std::array<double,3> Rec_RT_Phi(int iface, double x[],int dim);
//     std::array<double,3> Rec_RT_DPhi(int iface, double x[],int dim,int dir);

    /// \NAME RETURN SHAPE DERIVATIVE FUNCTIONS AT GAUSS POINTS ---------------------------
    ///< \param[in]  <dim>    Dimension
    ///< \param[in]  <qp>     Gaussian point
    ///< \param[in]  <InvJac> Jacobean
    ///< \param[out] <dphi>   Derivative

    void get_dphi_gl_g(const int qp, const double InvJac[], double dphi[]);
    void get_ddphi_gl_g(const int qp, const double InvJac[], double ddphi[]);
    void get_dphi_gl_g(const int dim, const int qp, const double InvJac[], double dphi[], int sdim);
    void get_dphi_gl_g(const int dim, const int qp, const double InvJac[], std::vector<double>& dphi);
    void get_dphi_node(const int dim, const int qp, const double InvJac[], double dphi[]);

    void get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]);
    void get_dphi_on_given_nodeG(
    const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType);
    void get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[]);


    //  RETURN SHAPE FUNCTIONS AT GAUSS POINTS --------------------------------------------
    ///< \param[in]  <dim>    Dimension
    ///< \param[in]  <qp>     Gaussian point
    ///< \param[out] <phi>    Shape function

    void get_phi_gl_g(const int qp, double phi[]);
    void get_phi_gl_g(const int qp, std::vector<double>& phi);
    void get_phi_g_arb_el(const int dim, const int qp, double phi[], int FamilyType);
    void get_dphi_g_arb_el(const int dim, const int qp, const double InvJac[], double dphi[], int FamilyType);

    
    
      double ComputeInverseMatrix(double Matrix[], double InvMatrix[], int Dimension);
      //===================================================================
    //  Jacobian functions -------------------------------------------------------------------
  //===================================================================
    // Jacobian and inverse Jacobian at gaussian points
    double JacobianOnGauss(const int dim,const int ng,const double ElemCoords[], double Jac[]);    ///< Jacobian (3D-2D)
    double JacobianOnPoint(const int dim,  const double Point[],  double ElemCoords[], double Jac[], int FamilyType, int nShape);
    //   double NablaCoordT2D(const int ng, double x[], double  Jac[]);
//   double NablaCoordT3D(const int ng,double xyz[],double Jac[]);                   // =====================================
    // double Jac3D(const int ng, double x[], double InvJac[]);                ///< Jacobian (3D)

    double Jac(const int ng, const double x[], double InvJac[]);                  ///< Jacobian (3D-2D)
    // double Jac1D(const int ng, double x[], double InvJac[]);                ///< Jacobian (1D)
    // double Jac2D(const int ng, double x[], double InvJac[]);                ///< Jacobian (2D)
    // double Jac2D(const int ng, double x[], double InvJac[], double Jac[]);  ///< Jacobian (2D)

    // Jacobian at gaussian points for arbitrary fem element
    double JacG(
    const int ng,const double x[], double InvJac[], int FamilyType, int dimension);  ///< Jacobian (3D-2D)
    // Jacobian at nodal points
    double Jac_nodes(const int ng, const double x[], double InvJac[]);    ///< Jacobian (3D-2D)
    // double Jac1D_nodes(const int ng, double x[], double InvJac[]);  ///< Jacobian (1D)
    // double Jac2D_nodes(const int ng, double x[], double InvJac[]);  ///< Jacobian (2D)
    // double Jac3D_nodes(const int ng, double x[], double InvJac[]);  ///< Jacobian (3D)
    double JacOnGivenCanCoords( const int dim, const double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape);

    // Jacobian dim -1
    double JacSur(const int ng, const double x[], double InvJac[]) ;  ///< Boundary Jacobian (2D-1D)
    // double JacSur2D(const int ng, double x[], double InvJac[]) const;
    // double JacSur3D(const int n_gauss, double x[], double InvJac[]) const;  ///< Surface Jacobian (2D)
    //   double JacSur2D(const int ng,double x[]) const;     ///<  Line Jacobian    (1D)
    // double JacSur1D(const int ng, double x[], double InvJac[]) const;  ///<  P Jacobian    (1D)

    // functions -----------------------------------------------------------------------
    /// Compute normal normal_g[] at xx[] point
    void normal_g(const double* xx, double* normal_g) const;  ///< unit normal to the surface
    void normal_g(const double* xx, const double x_c[], double* normal_g)
    const;  ///<  and check x_c[] being interior point
    void normal_g(const double* xx, const double x_c[], double* normal_g, int& sign_normal)
    const;  ///<  and check x_c[] being interior point
    void Oxy_face(const double* xx, const double x_c[],  double *normal_g,double* tg_g,
                  double* tg03)
    const;  ///<  and check x_c[] being interior point

    // Reading - writing ---------------------------------------------------------------
    // const std::string& file1 << file name to write over
    // const std::string& namefem, << name fem geometry element
    /// Write
    void write(const std::string& filename);
 void print_rec_xdmf(const std::string& filename,const std::string& namefem,const std::string& namelem);
//     void print_nodes_xdmf(const std::string& filename,const std::string& namefem);
//     void print_gauss_xdmf(const std::string& filename,const std::string& namefem);
//     void print_function_rec_xdmf(const std::string& filename,const std::string& namefem);
//     void print_function_tetra_xdmf(const std::string& filename,const std::string& namefem);
//     void print_Dfunction_rec_xdmf(const std::string& filename,const std::string& namefem);
//     void print_Dfunction_tetra_xdmf(const std::string& filename,const std::string& namefem);

//======

// ================================================================
// Special functions
// =========================================================
    double Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim);
    double dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim);
    double ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[]);
    double ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[]);




private:
    /// Reading
    void read_c(std::istream& infile);
    /// Writing
    void write_c(std::ostream& infile);

    /// \f$ \xi,\eta,\chi \f$ coordinates of the HEX27 nodes
    const int _CooH27[27 * 3] = {
        // 0  1   2    3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
        // 26
        -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
            1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
            -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
            1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
            -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
            0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
        };
    /// \f$ \xi,\eta \f$ coordinates of the QUAD9 nodes
    // note _CooQ9[i]=[_CooH27[i+12],_CooH27[27+i+12]
    const int _CooQ9[9 * 2] = {
        /**< 1   2   3   4   5   6   7   8   9  node */
        -1, 1,  1, -1, 0,  1, 0, -1, 0, /**< coordinate xi */
            -1, -1, 1, 1,  -1, 0, 1, 0,  0  /**< coordinate eta */
        };
    /// \f$ \xi \f$ coordinates of the EDGE3 nodes
    // Note _CooE3[i]=[_CooH27[i+6]]
    const int _CooE3[3] = {
        // 1   2   3
        -1, 1, 0  /// xi
        };


    const int _CooTetra10[40] = {
        // L_i = l(j)*(a* l(k) + b)
        //   j   k   a    b
        0, 0, 2, -1,  // 0
        1, 1, 2, -1,  // 1
        2, 2, 2, -1,  // 2
        3, 3, 2, -1,  // 3
        1, 0, 4, 0,   // 4
        1, 2, 4, 0,   // 5
        2, 0, 4, 0,   // 6
        0, 3, 4, 0,   // 7
        1, 3, 4, 0,   // 8
        2, 3, 4, 0    // 9
    };

    const int _CooTetra4[16] = {
        // L_i = a_i*xi + b_i*eta + c_i*chi + d_i
        //   a   b   c   d
        -1, -1, -1, 1,  // 0
            1,  0,  0,  0,  // 1
            0,  1,  0,  0,  // 2
            0,  0,  1,  0   // 3
        };

    const int _CooTriEl[9] = {
        //  note _CooTriEl[i]=_CooTetra4
        // L_i = a_i*xi + b_i*eta + c_i
        //   a   b   c
        -1, -1, 1,  // 0
            1,  0,  0,  // 1
            0,  1,  0   // 2
        };




    /// Offset for reading the generic node coordinates inside the _CooH27 array
    const int _H27Off = 27;
    /// Offset for reading the generic node coordinates inside the _CooQ9 array
    const int _Q9Off = 9;
    
    bool phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it = 30, double tol=1.e-12) ;



   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
 };
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// 
// 
// /**
//  * FELagrangeVec objects are used for working with vector-valued
//  * finite elements
//  *
//  * \author Paul T. Bauman
//  * \date 2013
//  */
// template <unsigned int Dim>
// class FELagrangeVec : public FE<Dim,LAGRANGE_VEC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a vector Lagrange finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FELagrangeVec(const FEType & fet) :
//     FE<Dim,LAGRANGE_VEC> (fet)
//   {}
// };
// 
// 
// /**
//  * FEL2LagrangeVec objects are used for working with vector-valued
//  * finite elements
//  *
//  * \author Alexander Lindsay
//  * \date 2023
//  */
// template <unsigned int Dim>
// class FEL2LagrangeVec : public FE<Dim,L2_LAGRANGE_VEC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a vector Lagrange finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEL2LagrangeVec(const FEType & fet) :
//     FE<Dim,L2_LAGRANGE_VEC> (fet)
//   {}
// };
// 
// 
// 
// /**
//  * FEHierarchicVec objects are used for working with vector-valued
//  * high-order finite elements
//  *
//  * \author Roy H. Stogner
//  * \date 2023
//  */
// template <unsigned int Dim>
// class FEHierarchicVec : public FE<Dim,HIERARCHIC_VEC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a vector Hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEHierarchicVec(const FEType & fet) :
//     FE<Dim,HIERARCHIC_VEC> (fet)
//   {}
// };
// 
// 
// /**
//  * FEHierarchicVec objects are used for working with vector-valued
//  * high-order piecewise-continuous finite elements
//  *
//  * \author Roy H. Stogner
//  * \date 2023
//  */
// template <unsigned int Dim>
// class FEL2HierarchicVec : public FE<Dim,L2_HIERARCHIC_VEC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a vector Hierarchic finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEL2HierarchicVec(const FEType & fet) :
//     FE<Dim,L2_HIERARCHIC_VEC> (fet)
//   {}
// };
// 
// 
// /**
//  * FENedelecOne objects are used for working with vector-valued
//  * Nedelec finite elements of the first kind.
//  *
//  * \author Paul T. Bauman
//  * \date 2013
//  */
// template <unsigned int Dim>
// class FENedelecOne : public FE<Dim,NEDELEC_ONE>
// {
// public:
//   /**
//    * Constructor. Creates a Nedelec finite element of the first kind
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FENedelecOne(const FEType & fet) :
//     FE<Dim,NEDELEC_ONE> (fet)
//   {}
// };
// 
// /**
//  * FEMonomialVec objects are used for working with vector-valued
//  * discontinuous finite elements
//  *
//  * \author Alex D. Lindsay
//  * \date 2019
//  */
// template <unsigned int Dim>
// class FEMonomialVec : public FE<Dim,MONOMIAL_VEC>
// {
// public:
// 
//   /**
//    * Constructor. Creates a vector Monomial finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEMonomialVec(const FEType & fet) :
//     FE<Dim,MONOMIAL_VEC> (fet)
//   {}
// };
// 
// /**
//  * FERaviartThomas objects are used for working with vector-valued
//  * Raviart-Thomas finite elements.
//  *
//  * \author Nuno Nobre & Karthikeyan Chockalingam
//  * \date 2023
//  */
// template <unsigned int Dim>
// class FERaviartThomas : public FE<Dim,RAVIART_THOMAS>
// {
// public:
// 
//   /**
//    * Constructor. Creates a Raviart-Thomas finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FERaviartThomas(const FEType & fet) :
//     FE<Dim,RAVIART_THOMAS> (fet)
//   {}
// };
// 
// /**
//  * FEL2RaviartThomas objects are used for working with vector-valued
//  * discontinuous Raviart-Thomas finite elements, e.g. when constructing
//  * hybridized methods
//  *
//  * \author Alex D. Lindsay
//  * \date 2023
//  */
// template <unsigned int Dim>
// class FEL2RaviartThomas : public FE<Dim,L2_RAVIART_THOMAS>
// {
// public:
// 
//   /**
//    * Constructor. Creates a Raviart-Thomas finite element
//    * to be used in dimension \p Dim.
//    */
//   explicit
//   FEL2RaviartThomas(const FEType & fet) :
//     FE<Dim,L2_RAVIART_THOMAS> (fet)
//   {}
// };
// 
// =================================================================================
//  * Provide Typedefs for various element types.
// ==================================================================================
// namespace FiniteElements
// {
// Convenient definition for a 2D Clough-Tocher finite element.
// typedef FEClough<2> FEClough2D;
// Convenient definition for a 1D Hierarchic finite element.
// typedef FE<1,HIERARCHIC> FEHierarchic1D;
// Convenient definition for a 2D Hierarchic finite element.
// typedef FE<2,HIERARCHIC> FEHierarchic2D;
// Convenient definition for a 3D Hierarchic finite element.
// typedef FE<3,HIERARCHIC> FEHierarchic3D;
// Convenient definition for a 1D Discontinuous Hierarchic finite element.
// typedef FE<1,L2_HIERARCHIC> FEL2Hierarchic1D;
// Convenient definition for a 2D Discontinuous Hierarchic finite element.
// typedef FE<2,L2_HIERARCHIC> FEL2Hierarchic2D;
// Convenient definition for a 3D Discontinuous Hierarchic finite element.
// typedef FE<3,L2_HIERARCHIC> FEL2Hierarchic3D;
// Convenient definition for a 1D Lagrange finite element.
// typedef FE<1,LAGRANGE> FELagrange1D;
// Convenient definition for a 2D Lagrange finite element.
// typedef FE<2,LAGRANGE> FELagrange2D;
// Convenient definition for a 3D Lagrange finite element.
// typedef FE<3,LAGRANGE> FELagrange3D;
// Convenient definition for a 1D Discontinuous Lagrange finite element.
// typedef FE<1,L2_LAGRANGE> FEL2Lagrange1D;
// Convenient definition for a 2D Discontinuous Lagrange finite element.
// typedef FE<2,L2_LAGRANGE> FEL2Lagrange2D;
// Convenient definition for a 3D Discontinuous Lagrange finite element.
// typedef FE<3,L2_LAGRANGE> FEL2Lagrange3D;
// Convenient definition for a 1D Monomial finite element.
// typedef FE<1,MONOMIAL> FEMonomial1D;
// Convenient definition for a 2D Monomial finite element.
// typedef FE<2,MONOMIAL> FEMonomial2D;
// Convenient definition for a 3D Monomial finite element.
// typedef FE<3,MONOMIAL> FEMonomial3D;
// }
// =================================================================================
// /**
//  * Helper functions for finite differenced derivatives in cases where
//  * analytical calculations haven't been done yet.
//  */
// template <typename OutputShape>
// OutputShape fe_fdm_deriv(const Elem * elem,
//                          const Order order,
//                          const unsigned int i,
//                          const unsigned int j,
//                          const Point & p,
//                          const bool add_p_level,
//                          OutputShape(*shape_func)
//                            (const Elem *, const Order,
//                             const unsigned int, const Point &,
//                             const bool));
// 
// template <typename OutputShape>
// OutputShape fe_fdm_deriv(const ElemType type,
//                          const Order order,
//                          const unsigned int i,
//                          const unsigned int j,
//                          const Point & p,
//                          OutputShape(*shape_func)
//                            (const ElemType, const Order,
//                             const unsigned int, const Point &));
// 
// template <typename OutputShape>
// OutputShape
// fe_fdm_second_deriv(const Elem * elem,
//                     const Order order,
//                     const unsigned int i,
//                     const unsigned int j,
//                     const Point & p,
//                     const bool add_p_level,
//                     OutputShape(*deriv_func)
//                       (const Elem *, const Order,
//                        const unsigned int, const unsigned int,
//                        const Point &, const bool));
// 
// template <typename OutputShape>
// OutputShape fe_fdm_second_deriv(const ElemType type,
//                                 const Order order,
//                                 const unsigned int i,
//                                 const unsigned int j,
//                                 const Point & p,
//                                 OutputShape(*deriv_func)
//                                   (const ElemType, const Order,
//                                    const unsigned int,
//                                    const unsigned int,
//                                    const Point &));
// // =================================================================================
//
// /**
//  * Helper functions for Lagrange-based basis functions.
//  */
// void lagrange_nodal_soln(const Elem * elem,
//                          const Order order,
//                          const std::vector<Number> & elem_soln,
//                          std::vector<Number> &       nodal_soln,
//                          bool add_p_level = true);
// 
// /**
//  * Helper functions for Discontinuous-Pn type basis functions.
//  */
// unsigned int monomial_n_dofs(const ElemType t, const Order o);
// 
// /**
//  * Helper functions for rational basis functions.
//  */
// // shapes[i][j] is shape function phi_i at point p[j]
// void rational_fe_weighted_shapes(const Elem * elem,
//                                  const FEType underlying_fe_type,
//                                  std::vector<std::vector<Real>> & shapes,
//                                  const std::vector<Point> & p,
//                                  const bool add_p_level);
// 
// // shapes[i][q] is shape function phi_i at point p[q]
// // derivs[j][i][q] is dphi_i/dxi_j at p[q]
// void rational_fe_weighted_shapes_derivs(const Elem * elem,
//                                         const FEType fe_type,
//                                         std::vector<std::vector<Real>> & shapes,
//                                         std::vector<std::vector<std::vector<Real>>> & derivs,
//                                         const std::vector<Point> & p,
//                                         const bool add_p_level);
// 
// Real rational_fe_shape(const Elem & elem,
//                        const FEType underlying_fe_type,
//                        const unsigned int i,
//                        const Point & p,
//                        const bool add_p_level);
// 
// Real rational_fe_shape_deriv(const Elem & elem,
//                              const FEType underlying_fe_type,
//                              const unsigned int i,
//                              const unsigned int j,
//                              const Point & p,
//                              const bool add_p_level);
// 
// Real rational_fe_shape_second_deriv(const Elem & elem,
//                                     const FEType underlying_fe_type,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const Point & p,
//                                     const bool add_p_level);
// 
// void rational_all_shapes (const Elem & elem,
//                           const FEType underlying_fe_type,
//                           const std::vector<Point> & p,
//                           std::vector<std::vector<Real>> & v,
//                           const bool add_p_level);
// 
// template <typename OutputShape>
// void rational_all_shape_derivs (const Elem & elem,
//                                 const FEType underlying_fe_type,
//                                 const std::vector<Point> & p,
//                                 std::vector<std::vector<OutputShape>> * comps[3],
//                                 const bool add_p_level);
// 
// } // namespace libMesh
//  // =================================================================================
// #define LIBMESH_DEFAULT_VECTORIZED_FE(MyDim, MyType) \
// template<>                                           \
// void FE<MyDim,MyType>::all_shapes                    \
//   (const Elem * elem,                                \
//    const Order o,                                    \
//    const std::vector<Point> & p,                     \
//    std::vector<std::vector<OutputShape>> & v,        \
//    const bool add_p_level)                           \
// {                                                    \
//   FE<MyDim,MyType>::default_all_shapes               \
//     (elem,o,p,v,add_p_level);                        \
// }                                                    \
//                                                      \
// template<>                                           \
// void FE<MyDim,MyType>::shapes                        \
//   (const Elem * elem,                                \
//    const Order o,                                    \
//    const unsigned int i,                             \
//    const std::vector<Point> & p,                     \
//    std::vector<OutputShape> & v,                     \
//    const bool add_p_level)                           \
// {                                                    \
//   FE<MyDim,MyType>::default_shapes                   \
//     (elem,o,i,p,v,add_p_level);                      \
// }                                                    \
//                                                      \
// template<>                                           \
// void FE<MyDim,MyType>::shape_derivs                  \
//   (const Elem * elem,                                \
//    const Order o,                                    \
//    const unsigned int i,                             \
//    const unsigned int j,                             \
//    const std::vector<Point> & p,                     \
//    std::vector<OutputShape> & v,                     \
//    const bool add_p_level)                           \
// {                                                    \
//   FE<MyDim,MyType>::default_shape_derivs             \
//     (elem,o,i,j,p,v,add_p_level);                    \
// }                                                    \
//                                                      \
// template<>                                           \
// void FE<MyDim,MyType>::all_shape_derivs              \
//   (const Elem * elem,                                \
//    const Order o,                                    \
//    const std::vector<Point> & p,                     \
//    std::vector<std::vector<OutputShape>> * comps[3], \
//    const bool add_p_level)                           \
// {                                                    \
//   FE<MyDim,MyType>::default_all_shape_derivs         \
//     (elem,o,p,comps,add_p_level);                    \
// }
// 
// 
// #endif // LIBMESH_FE_H
