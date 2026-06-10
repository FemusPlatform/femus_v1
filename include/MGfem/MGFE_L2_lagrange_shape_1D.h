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
#pragma once
// 
// // Local includes
#include "MGFE_L0_II_enum_order.h" // FIRST, SECOND, etc.



// ================================================== 
inline  double fe_lagrange_1D_linear_shape(
  const unsigned int i, 
  const double xi
){// ==================================================
  assert(i< 2);
  switch (i)  {
    case 0:   return .5*(1. - xi);
    case 1:   return .5*(1. + xi);
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}


// =================================================
inline
double fe_lagrange_1D_quadratic_shape(
  const unsigned int i,
  const double xi
){// =================================================
  assert (i< 3);
  switch (i) {
    case 0: return .5*xi*(xi - 1.);
    case 1: return .5*xi*(xi + 1);
    case 2: return (1. - xi*xi);
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}


// =================================================
inline double fe_lagrange_1D_cubic_shape(
  const unsigned int i,
  const double xi
){// =================================================
  assert(i< 4);
  switch (i) {
    case 0:  return 9./16.*(1./9.-xi*xi)*(xi-1.);
    case 1:  return -9./16.*(1./9.-xi*xi)*(xi+1.);
    case 2:  return 27./16.*(1.-xi*xi)*(1./3.-xi);
    case 3:  return 27./16.*(1.-xi*xi)*(1./3.+xi);
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}


// ================================================
inline double fe_lagrange_1D_shape(
  const Order order,
  const unsigned int i,
  const double xi
){// ================================================
  switch (order)    {
    case CONSTANT:  return 1.; 
      // Lagrange linears
    case FIRST:    return fe_lagrange_1D_linear_shape(i, xi);
      // Lagrange quadratics
    case SECOND:  return fe_lagrange_1D_quadratic_shape(i, xi);
      // Lagrange cubics
    case THIRD: return fe_lagrange_1D_cubic_shape(i, xi);
    default:
      std::cerr<<"ERROR: Unsupported polynomial order = " << order; abort();return 0.;
    }
}

// FIRST DERIVATIVE
// ================================================
inline double fe_lagrange_1D_linear_shape_deriv(
  const unsigned int i,
  const unsigned int j,
  const double
){// ================================================
  assert(j== 0); assert (i< 2); // only d()/dxi in 1D!
  switch (i)  {
    case 0:  return -.5;
    case 1:  return .5;
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}

// ================================================
inline double fe_lagrange_1D_quadratic_shape_deriv(
  const unsigned int i,
  const unsigned int j,
  const double xi
){// ================================================
  assert (j== 0); assert(i< 3);// only d()/dxi in 1D!
  switch (i) {
    case 0: return xi-.5;
    case 1: return xi+.5;
    case 2: return -2.*xi;
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}

// ================================================
inline
double fe_lagrange_1D_cubic_shape_deriv(
  const unsigned int i,
  const unsigned int j,
  const double xi
){ // only d()/dxi in 1D!============================
  assert (j ==0); assert (i< 4);
  switch (i) {
    case 0:  return -9./16.*(3.*xi*xi-2.*xi-1./9.);
    case 1:  return -9./16.*(-3.*xi*xi-2.*xi+1./9.);
    case 2:  return 27./16.*(3.*xi*xi-2./3.*xi-1.);
    case 3:  return 27./16.*(-3.*xi*xi-2./3.*xi+1.);
    default:  std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}


// ================================================
inline double fe_lagrange_1D_shape_deriv(
  const Order order,
 const unsigned int i,
  const unsigned int j,
 const double xi
){// ================================================
  switch (order)    {
    case CONSTANT: return 0.;
    case FIRST: return fe_lagrange_1D_linear_shape_deriv(i, j, xi);
    case SECOND: return fe_lagrange_1D_quadratic_shape_deriv(i, j, xi);
    case THIRD: return fe_lagrange_1D_cubic_shape_deriv(i, j, xi);
    default:  
     std::cerr<<"ERROR: Unsupported polynomial order = " << order;abort(); return 0.;
    }
}


// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// fe_lagrange_1D_linear_shape_second_deriv is 0
// =====================================================
inline double fe_lagrange_1D_quadratic_shape_second_deriv(
  const unsigned int i,
  const unsigned int j,
  const double
){// ================================================
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!
  assert(j== 0);
  switch (i)   {
    case 0:  return 1.;
    case 1:  return 1.;
    case 2: return -2.;
    default:   std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}

// =====================================================
inline double fe_lagrange_1D_cubic_shape_second_deriv(
  const unsigned int i,
  const unsigned int j,
  const double xi
){// ================================================
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!
  assert(j== 0);
  switch (i) {
    case 0: return -9./16.*(6.*xi-2);
    case 1: return -9./16.*(-6*xi-2.);
    case 2: return 27./16.*(6*xi-2./3.);
    case 3:  return 27./16.*(-6*xi-2./3.);
    default:   std::cerr<<"Invalid shape function index i = " << i; return 0.;
    }
}


// ================================================
inline double fe_lagrange_1D_shape_second_deriv(
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double xi
){// ================================================
  switch (order){ // All second derivatives of linears are zero....
   case CONSTANT: return 0.;
    case FIRST: return 0.;
    case SECOND:  return fe_lagrange_1D_quadratic_shape_second_deriv(i, j, xi);
    case THIRD:   return fe_lagrange_1D_cubic_shape_second_deriv(i, j, xi);
    default:   std::cerr<<"ERROR: Unsupported polynomial order = " << order;abort();return 0.;
    } // end switch (order)
}
