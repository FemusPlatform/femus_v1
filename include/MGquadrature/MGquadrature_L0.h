// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#pragma once

// Local includes
// #include "libmesh/libmesh_common.h"
// #include "libmesh/reference_counted_object.h"
// #include "libmesh/point.h"
#include "MGFE_L0_II_enum_order.h" // INVALID_ELEM
// #include "MGenum_order.h" // INVALID_ORDER

// C++ includes
#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <iostream>
#include <cassert>

// forward declarations
// class Elem;
enum QuadratureType : int;

/**
 * The \p QBase class provides the basic functionality from which
 * various quadrature rules can be derived.  It computes and stores
 * the quadrature points (in reference element space) and associated
 * weights.
 */
class QBase {
protected:

  /** Constructor.  Use the build() method instead.  */
  QBase (unsigned int dim, Order order=INVALID_ORDER);
  unsigned int _n_points;
public:

  /* Constructor destructor */
  QBase (const QBase &) = default;            // Copy constructor
  QBase (QBase &&) = default;                 // Copy/move ctor
  QBase & operator= (const QBase &) = default;// copy assignment operator
  QBase & operator= (QBase &&) = default;     // move assignment operator
  virtual ~QBase() = default;                 // destructor

  /// \returns
  virtual std::unique_ptr<QBase> clone() const;// A copy of this quadrature rule wrapped in a smart pointer.
  virtual QuadratureType type() const = 0; // The quadrature type in derived classes.

  /**
   * Builds a specific quadrature rule based on the \p name
   * string. This enables selection of the quadrature rule at
   * run-time.  The input parameter \p name must be mappable through
   * the \p Utility::string_to_enum<>() function.
   *
   * This function allocates memory, therefore a \p std::unique_ptr<QBase>
   * is returned so that the user does not accidentally leak it.
   */
  static std::unique_ptr<QBase> build (std::string_view name,
                                       const unsigned int dim,
                                       const Order order=INVALID_ORDER);

  /**
   * Builds a specific quadrature rule based on the QuadratureType.
   * This enables selection of the quadrature rule at run-time.
   *
   * This function allocates memory, therefore a \p std::unique_ptr<QBase>
   * is returned so that the user does not accidentally leak it.
   */
  static std::unique_ptr<QBase> build (const QuadratureType qt,
                                       const unsigned int dim,
                                       const Order order=INVALID_ORDER);

  /* \returns The element type we're currently using.  */
  ElemType get_elem_type() const { return _type; }
  /* \returns The p-refinement level we're currently using. */
  unsigned int get_p_level() const { return _p_level; }
  /* \returns The number of points associated with the quadrature rule.  */
  unsigned int n_points() const  {   assert (_n_points!=0);    return _n_points;  }
  /* Alias for n_points() to enable use in index_range
   * \returns The number of points associated with the quadrature rule.  */
  unsigned int size() const {    return n_points();  }
  /* \returns The spatial dimension of the quadrature rule.  */
  unsigned int get_dim() const { return _dim; }
  // /* A \p std::vector containing the quadrature point locations in reference element space. */
   /*const std::vector<Point> &*/void get_points(double points[]) const { points=_points;return ; }
  // /* \returns A \p std::vector containing the quadrature point locations
  //  * in reference element space as a writable reference. */
  // std::vector<Point> & get_points() { return _points; }
  /* A constant reference to a \p std::vector containing the quadrature weights.  */
  /*const std::vector<Real> &*/ void get_weights(double weights[]) const {  weights=_weights; return;}
  /* \returns A writable references to a \p std::vector containing the quadrature weights.  */
 /* std::vector<Real> &*/void get_weights(double weights[]) {  weights=_weights; return;}
  /* \returns The \f$ i^{th} \f$ quadrature point in reference element space.  */
  /*Point*/ void qp(double point[], const unsigned int i) const {assert (i< _n_points); for(int j=0;j<_dim;j++) point[j]=_points[i*_dim+j];  return; }
  /* \returns The \f$ i^{th} \f$ quadrature weight.  */
  double w(const unsigned int i) const {assert (i<  _n_points);    return _weights[i];  }

  /* Initializes the data structures for a quadrature rule for an element of type \p type.  */
  virtual void init (const ElemType type=INVALID_ELEM,  unsigned int p_level=0);

  /**
   * Initializes the data structures for an element potentially "cut"
   * by a signed distance function.  The array \p vertex_distance_func
   * contains vertex values of the signed distance function.  If the
   * signed distance function changes sign on the vertices, then the
   * element is considered to be cut.) This interface can be extended
   * by derived classes in order to subdivide the element and construct
   * a composite quadrature rule.
   */
  // virtual void init (const Elem & elem,
  //                    const std::vector<Real> & vertex_distance_func,
  //                    unsigned int p_level=0);

  /**
   * \returns The current "total" order of the quadrature rule which
   * can vary element by element, depending on the Elem::p_level(),
   * which gets passed to us during init().
   *
   * Each additional power of p increases the quadrature order
   * required to integrate the mass matrix by 2, hence the formula
   * below.
   *
   * \todo This function should also be used in all of the Order
   * switch statements in the rules themselves.
   */
  Order get_order() const { return static_cast<Order>(_order + 2 * _p_level); }

  /**
   * Prints information relevant to the quadrature rule, by default to
   * libMesh::out.
   */
  void print_info(std::ostream & os=std::cout) const;

  /**
   * Maps the points of a 1D quadrature rule defined by "old_range" to
   * another 1D interval defined by "new_range" and scales the weights
   * accordingly.
   */
  void scale(std::pair<double, double> old_range, std::pair<double, double> new_range);

  /* Same as above, but allows you to use the stream syntax.  */
  friend std::ostream & operator << (std::ostream & os, const QBase & q);

  /* \returns \p true if the shape functions need to be recalculated, \p false otherwise.
   * This may be required if the number of quadrature points or their position changes.  */
  virtual bool shapes_need_reinit() { return false; }

  /**
   * Flag (default true) controlling the use of quadrature rules with
   * negative weights.  Set this to false to require rules with all
   * positive weights.
   *
   * Rules with negative weights can be unsuitable for some problems.
   * For example, it is possible for a rule with negative weights to
   * obtain a negative result when integrating a positive function.
   *
   * A particular example: if rules with negative weights are not allowed,
   * a request for TET,THIRD (5 points) will return the TET,FIFTH (14 points)
   * rule instead, nearly tripling the computational effort required!
   */
  bool allow_rules_with_negative_weights;

  /**
   * The flag's value defaults to false so that one does not accidentally
   * use a nodal quadrature rule on Pyramid elements, since evaluating the
   * inverse element Jacobian (e.g. dphi) is not well-defined at the
   * Pyramid apex because the element Jacobian is zero there.
   *
   * We do not want to completely prevent someone from using a nodal
   * quadrature rule on Pyramids, however, since there are legitimate
   * use cases (lumped mass matrix) so the flag can be set to true to
   * override this behavior.
   */
  bool allow_nodal_pyramid_quadrature;

protected:


  /**
   * Initializes the 0D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  Generally this
   * is just one point with weight 1.
   *
   * \note The arguments should no longer be used for anything in
   * derived classes, they are only maintained for backwards
   * compatibility and will eventually be removed.
   */
  virtual void init_0D (const ElemType type=INVALID_ELEM, unsigned int p_level=0);
  virtual void init_1D (const ElemType type=INVALID_ELEM, unsigned int p_level=0) = 0;
  virtual void init_2D (const ElemType type=INVALID_ELEM, unsigned int p_level=0);
  virtual void init_3D (const ElemType type=INVALID_ELEM, unsigned int p_level=0);

  /* Constructs a 2D rule from the tensor product of \p q1D with
   * itself.  Used in the \p init_2D() routines for quadrilateral
   * element types.
   */
  void tensor_product_quad (QBase & q1D); //quadrilateral element types.
  void tensor_product_hex (const QBase & q1D); //hexahedral element types
  void tensor_product_prism (const QBase & q1D, const QBase & q2D); // prismatic

  /* The spatial dimension of the quadrature rule. */
  unsigned int _dim;
  /* The polynomial order which the quadrature rule is capable of integrating exactly.  */
  Order _order;
  /* The type of element for which the current values have been computed. */
  ElemType _type;
  /* The p-level of the element for which the current values have been computed.  */
  unsigned int _p_level;
  /* The locations of the quadrature points in reference element space. */
  /*std::vector<Point>*/ double *_points;
  /* The quadrature weights.  The order of the weights matches the ordering of the _points vector.  */
  /*std::vector<Real>*/ double *_weights;
};



// #endif // LIBMESH_QUADRATURE_H
