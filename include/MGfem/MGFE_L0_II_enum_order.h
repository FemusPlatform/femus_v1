// // The libMesh Finite Element Library.
// // Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
// 
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// 
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// // Lesser General Public License for more details.
// 
// // You should have received a copy of the GNU Lesser General Public
// // License along with this library; if not, write to the Free Software
// // Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// 
// 
// #ifndef LIBMESH_ENUM_ORDER_H
// #define LIBMESH_ENUM_ORDER_H
// 
#pragma once
// namespace libMesh {
// 
// /**
//  * \enum libMesh::Order defines an \p enum for polynomial orders.
//  * Fixing each label to a specific int, since \p InfFE and p refinement
//  * may cast between them.
//  *
//  * \note It is possible to safely cast integral values larger than 43
//  * to Order, provided that the cast value fits in the range of int. A
//  * potential use case for such a cast is when employing the QGrid
//  * class at extremely high order.
//  *
//  * The fixed type, i.e. ": int", enumeration syntax used here allows
//  * this enum to be forward declared as
//  * enum Order : int;
//  * reducing header file dependencies.
//  */
enum Order : int {
            CONSTANT     =  0,
            FIRST        =  1,
            SECOND       =  2,
            THIRD        =  3,
            FOURTH       =  4,
            FIFTH        =  5,
            SIXTH        =  6,
            SEVENTH      =  7,
            EIGHTH       =  8,
            NINTH        =  9,
            TENTH        = 10,
            ELEVENTH     = 11,
            TWELFTH      = 12,
            THIRTEENTH   = 13,
            FOURTEENTH   = 14,
            FIFTEENTH    = 15,
            SIXTEENTH    = 16,
            SEVENTEENTH  = 17,
            EIGHTTEENTH  = 18, // misspelling, deprecated
            EIGHTEENTH  = 18,
            NINETEENTH   = 19,
            TWENTIETH    = 20,
            TWENTYFIRST   = 21,
            TWENTYSECOND  = 22,
            TWENTYTHIRD   = 23,
            TWENTYFOURTH  = 24,
            TWENTYFIFTH   = 25,
            TWENTYSIXTH   = 26,
            TWENTYSEVENTH = 27,
            TWENTYEIGHTH  = 28,
            TWENTYNINTH   = 29,
            THIRTIETH     = 30,
            THIRTYFIRST   = 31,
            THIRTYSECOND  = 32,
            THIRTYTHIRD   = 33,
            THIRTYFOURTH  = 34,
            THIRTYFIFTH   = 35,
            THIRTYSIXTH   = 36,
            THIRTYSEVENTH = 37,
            THIRTYEIGHTH  = 38,
            THIRTYNINTH   = 39,
            FORTIETH     = 40,
            FORTYFIRST   = 41,
            FORTYSECOND  = 42,
            FORTYTHIRD   = 43,
            MAXIMUM = 255,
            // Invalid
            INVALID_ORDER};
// 
//   // Standardize this so nvc++ and clang -fsanitize=integer don't
//   // complain about all the ways in which we might do it wrong
template <typename T>  inline Order operator+(Order o, T p)  {    return static_cast<Order>(static_cast<int>(o) + int(p));  }
template <typename T>  inline Order operator-(Order o, T p)  {    return static_cast<Order>(static_cast<int>(o) - int(p));  }
template <typename T>  inline Order operator+(T p, Order o)  {    return o + p;  }
template <typename T>  inline Order & operator+=(Order &o, T p)  {    o = o + p;    return o;  }
template <typename T>  inline Order & operator-=(Order &o, T p)  {    o = o - p;    return o;  }
  
  
  
  
  
enum FEContinuity : int {  DISCONTINUOUS,   C_ZERO,   C_ONE,   H_CURL,   H_DIV,SIDE_DISCONTINUOUS};

enum FEFamily : int {
               // C0
               LAGRANGE     = 0,   HIERARCHIC   = 1, LAGRANGE_XFEM   = 40,
               // discontinuous, in local coordinates
               MONOMIAL      = 2,  L2_HIERARCHIC = 6, L2_LAGRANGE   = 7,
               // higher-order
               BERNSTEIN    = 3,  SZABAB       = 4,
               // discontinuous, in global coordinates
               XYZ          = 5,
               // C1 elements
               CLOUGH       = 21,   HERMITE      = 22,    SUBDIVISION  = 23,
               // A scalar variable that couples to
               // all other DOFs in the system
               SCALAR       = 31,
               // Vector-valued elements
               LAGRANGE_VEC   = 41,   NEDELEC_ONE    = 42,   MONOMIAL_VEC   = 43,     RAVIART_THOMAS = 44,
               L2_RAVIART_THOMAS = 45,    L2_LAGRANGE_VEC = 46,     HIERARCHIC_VEC = 47,   L2_HIERARCHIC_VEC = 48,
               // Rational basis functions
               RATIONAL_BERNSTEIN = 61,
               // Side-supported basis functions
               SIDE_HIERARCHIC = 71,
               // Invalid
               INVALID_FE   = 99};
               
               
   enum ElemType : int {
               // 1D
               EDGE2 = 0,  EDGE3 = 1,  EDGE4 = 2,
               // 2D
               TRI3 = 3,  TRI6 = 4, QUAD4 = 5,  QUAD8 = 6,  QUAD9 = 7,
               // 3D
               TET4 = 8,TET10 = 9, HEX8 = 10,  HEX20 = 11,  HEX27 = 12,  PRISM6 = 13,   PRISM15 = 14,  PRISM18 = 15,   PYRAMID5 = 16,
               PYRAMID13 = 17,   PYRAMID14 = 18,
               D1CONST =20,  D2CONST=23, D3CONST=24,
               // 0D
               NODEELEM = 27,
               // Miscellaneous Elems
               REMOTEELEM = 28,          TRI3SUBDIVISION = 29,
               // Shell Elems
               TRISHELL3 = 30,   QUADSHELL4 = 31, QUADSHELL8 = 32,
               // Elems with Tri7 (Tri with mid-face node) faces
               TRI7 = 33,  TET14 = 34,  PRISM20 = 35,  PRISM21 = 36,  PYRAMID18 = 37,  // Another shell elem
               QUADSHELL9 = 38,
               // Invalid
               INVALID_ELEM};   // should always be last          
                
                              
               
               /**
              
 * class FEType hides (possibly multiple) FEFamily and approximation
 * orders, thereby enabling specialized finite element families.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Manages the family, order, etc. parameters for a given FE.
 */
class FEType{ 
public:
   //  data -------------------------------------->  
   Order order;/* The approximation order of the element.  */
   FEFamily family; /* The type of finite element.  Valid types are \p LAGRANGE, \p HIERARCHIC, etc...  */
   ElemType fetype;
   // Constructor.  --------------------------------------> 
   FEType(const Order      o = Order::FIRST,   const FEFamily f = LAGRANGE) :    order(o),    family(f)  {}

  /* Tests equality */
  void set(ElemType fetype_in){fetype= fetype_in;}
  bool operator== (const FEType & f2) const  {  return (order == f2.order  && family == f2.family  );  }
  bool operator!= (const FEType & f2) const  {    return !(*this == f2);  }
  /*An ordering to make FEType useful as a std::map key   */
  bool operator< (const FEType & f2) const  {    
    if (order != f2.order)      return (order < f2.order);    
    if (family != f2.family)    return (family < f2.family);
    return false;
  }

  // get ---------------------------------------------->
  /**
   * \returns The default quadrature order for this \p FEType.  The
   * default quadrature order is calculated assuming a polynomial of
   * degree \p order and is based on integrating the mass matrix for
   * such an element exactly on affine elements.
   */
  Order default_quadrature_order () const{return order;};
  FEFamily get_family()const{return  family;}
  ElemType get_fetype() const{return  fetype;}
};

               
// #endif
