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

// std include ------------>
#include <cmath>
#include <array>
// // femus Local includes ----------->
#include "MGFE_L2_FE_template.h"
#include "MGFE_L0_C.h"
#include "MGFE_L2_lagrange_shape_1D.h"


// 
// LIBMESH_DEFAULT_VECTORIZED_FE(1,LAGRANGE)
// LIBMESH_DEFAULT_VECTORIZED_FE(1,L2_LAGRANGE)
// 
// // FE class members
template <>
 MGFE<1,LAGRANGE>::MGFE (const FEType & fet) :
   MGFEbase (1,fet),
   // _dim(1),             ///< space dimension
   //  _order(0),         ///< pol approx (2=quad  1=lin 0=k)
   //  _NoShape(0),      ///< number of shapes
   //  _GeomType(0),  ///< geom type 0=hex/quad  1=tetra/tri
   //  _FamType(0), 
   last_side(INVALID_ELEM),
   last_edge(INVALID_ELEM)
 {
   // Sanity check.  Make sure the
   // Family specified in the template instantiation
   // matches the one in the FEType object
   assert (LAGRANGE == this->get_family());
 }
 
 template <>
 MGFE<1,LAGRANGE>::~MGFE(){}
 
 
  template <>  const int MGFE<1,LAGRANGE>::get_dim()       {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 template <>  const int MGFE<1,LAGRANGE>::get_order()     {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 template <>  const int MGFE<1,LAGRANGE>::get_NoShape()   {return _NoShape;}  ///< Number of shape functions of the fem
 template <>  const int MGFE<1,LAGRANGE>::get_NoElNodes() {return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 template <>  const int MGFE<1,LAGRANGE>::get_GeomType()  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 template <>  const int MGFE<1,LAGRANGE>::get_FamType()   {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
  template <>  const int     MGFE<1,LAGRANGE>::get_NoGauss()     {return _NoGauss;} 
  template <> const double   MGFE<1,LAGRANGE>::get_xg(int i)     {return  _xg[i];}          ///< gaussian points
  template <> const double   MGFE<1,LAGRANGE>::get_weight(int i) {return _weight[i];}            ///< Weight
  template <> const double & MGFE<1,LAGRANGE>::get_phi()         {return *_phi_map;}           ///< Shape functions
  template <> const double & MGFE<1,LAGRANGE>::get_dphidxez_map(){return *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  template <> const double & MGFE<1,LAGRANGE>::get_dphidxx_map() {return *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
 
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
  template <>  const double &  MGFE<1,LAGRANGE>::get_xnodes()           {return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
  template <>  const int    &  MGFE<1,LAGRANGE>::get_fem_conn()         { return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  // template <> const int     MGFE<1,LAGRANGE>::get_deg(){return _deg;}                  ///< Polinomial degree  
  template <>  const double & MGFE<1,LAGRANGE>::get_phi_map_nodes()     {return *_phi_map_nodes;   }          ///< Shape functions
  template <>  const double & MGFE<1,LAGRANGE>::get_dphidxez_map_nodes(){return *_dphidxez_map_nodes;} 
 
 
 template <>    void MGFE<1,LAGRANGE>::get_dphi_gl_g(const int qp, const double InvJac[], double dphi[]);
 template <>   void MGFE<1,LAGRANGE>::get_ddphi_gl_g(const int qp, const double InvJac[], double ddphi[]);
template <>    void MGFE<1,LAGRANGE>::get_dphi_gl_g(const int dim, const int qp, const double InvJac[], double dphi[], int sdim);
template <>    void MGFE<1,LAGRANGE>::get_dphi_gl_g(const int dim, const int qp, const double InvJac[], std::vector<double>& dphi);
template <>    void MGFE<1,LAGRANGE>::get_dphi_node(const int dim, const int qp, const double InvJac[], double dphi[]);

 template <>   void MGFE<1,LAGRANGE>::get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]);
 template <>    void MGFE<1,LAGRANGE>::get_dphi_on_given_nodeG(
    const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType){}
template <>    void MGFE<1,LAGRANGE>::get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[]){}


  

 template <>   void MGFE<1,LAGRANGE>::get_phi_gl_g(const int qp, double phi[]);
 template <>   void MGFE<1,LAGRANGE>::get_phi_gl_g(const int qp, std::vector<double>& phi);
 template <>   void MGFE<1,LAGRANGE>::get_phi_g_arb_el(const int dim, const int qp, double phi[], int FamilyType){}
  template <>  void MGFE<1,LAGRANGE>::get_dphi_g_arb_el(const int dim, const int qp, const double InvJac[], double dphi[], int FamilyType){}
 
  template <>    void  MGFE<1,LAGRANGE>::normal_g(double const*, double const*, double*) const{return;}
  template <>    void  MGFE<1,LAGRANGE>::normal_g(double const*, double const*, double*, int &) const{return;}
  template <>    void MGFE<1,LAGRANGE>::write(const std::string& filename){return ;}
 
 template <>    int MGFE<1,LAGRANGE>::GetFamilyType(int elem_dof, int dim){return 1;}
 template <>    void MGFE<1,LAGRANGE>::print_rec_xdmf(const std::string& filename,const std::string& namefem,const std::string& namelem){return ;}
 
 
  // ========================================================
 // Jacobian  ----------------------------------------------------
 // ========================================================
 template <>    double MGFE<1,LAGRANGE>::JacG(const int,const double*, double*, int, int){return 0.;}
 template <>    double MGFE<1,LAGRANGE>::JacobianOnGauss(const int,const int,const double*, double*){return 0.;}
 template <>    double MGFE<1,LAGRANGE>::JacobianOnPoint(const int,const double*, double*, double*, int, int){return 0. ;}
  // ========================================================
  template <>    double MGFE<1,LAGRANGE>::Jac(
     const int ng, 
     const double* x , 
     double* InvJac
  ){// ========================================================
      double x_xi = 0.;
      int nshape = _NoShape;
      for(int s = 0; s < nshape; s++)  x_xi += x[s] * _dphidxez_map[ s * _NoGauss + ng];
      double det=  x_xi;
      InvJac[0] = 1./(x_xi+1.e-20);  // deta dy

      return (det);
  }
  
  // quadrature rule -----------------------------------------------
 template <>
void MGFE<1,LAGRANGE>::attach_quadrature_rule (QBase *q)
 {
assert(q);
   this->_qrule = q;
   // make sure we don't cache results from a previous quadrature rule
   this->_elem_type = INVALID_ELEM;
   return;
 }
 template <>
 unsigned int MGFE<1,LAGRANGE>::n_shape_functions () const {
    return MGFE<1,LAGRANGE>::n_dofs (this->_elem_type, this->_fe_type.order + this->_p_level);
  }
 
 
// ===========================================================================
//       Shape 
// ===========================================================================
template <> double MGFE<1,LAGRANGE>::shape(
  const ElemType,
  const Order order,
  const unsigned int i,
  const double p[]
) {
   return fe_lagrange_1D_shape(order, i, p[0]);
 } // ===========================================================================
// 
// template <> Real  MGFE<1,L2_LAGRANGE>::shape(const ElemType,   const Order order,const unsigned int i,const Point & p)
// {   return fe_lagrange_1D_shape(order, i, p(0)); }
//------------
// template <> Real  MGFE<1,LAGRANGE>::shape(const Elem * elem, const Order order,  const unsigned int i,
//                            const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);   return fe_lagrange_1D_shape(order + add_p_level*elem->p_level(), i, p(0)); } 
// -------------------
// template <> Real  MGFE<1,LAGRANGE>::shape(const FEType fet,  const Elem * elem,   const unsigned int i,
//                            const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_1D_shape(fet.order + add_p_level*elem->p_level(), i, p(0)); }
// ------------------
// template <> Real  MGFE<1,L2_LAGRANGE>::shape(const Elem * elem,   const Order order,
//                               const unsigned int i,    const Point & p,   const bool add_p_level)
// { libmesh_assert(elem);  return fe_lagrange_1D_shape(order + add_p_level*elem->p_level(), i, p(0)); }
// -----------
// template <> Real  MGFE<1,L2_LAGRANGE>::shape(const FEType fet,  const Elem * elem,   const unsigned int i,
//                               const Point & p,  const bool add_p_level)
// {  libmesh_assert(elem); return fe_lagrange_1D_shape(fet.order + add_p_level*elem->p_level(), i, p(0)); }
// 
// 
// ===========================================================================
//       Shape derivative (1st)
// ===========================================================================
// ==================================================================
 template <>
 double MGFE<1,LAGRANGE>::shape_deriv(
   const ElemType, 
   const Order order, 
   const unsigned int i,
   const unsigned int j,
   const double  p[]
 ){ return fe_lagrange_1D_shape_deriv(order, i, j, p[0]);}

// template <> Real  MGFE<1,L2_LAGRANGE>::shape_deriv(const ElemType,  const Order order,
//                                     const unsigned int i,  const unsigned int j,  const Point & p)
// {   return fe_lagrange_1D_shape_deriv(order, i, j, p(0)); }
// 
// 
// 
// template <>  Real  MGFE<1,LAGRANGE>::shape_deriv(const Elem * elem,  const Order order,    const unsigned int i,
//                                  const unsigned int j, const Point & p,   const bool add_p_level)
// { libmesh_assert(elem);   return fe_lagrange_1D_shape_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,L2_LAGRANGE>::shape_deriv(const Elem * elem, const Order order,  const unsigned int i,
//                                     const unsigned int j,  const Point & p, const bool add_p_level)
// {  libmesh_assert(elem);  return fe_lagrange_1D_shape_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,LAGRANGE>::shape_deriv(const FEType fet,  const Elem * elem,   const unsigned int i,
//                                  const unsigned int j,  const Point & p,  const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_1D_shape_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }

// template <> Real  MGFE<1,L2_LAGRANGE>::shape_deriv(const FEType fet, const Elem * elem,
//                                     const unsigned int i, const unsigned int j, const Point & p,
//                                     const bool add_p_level) {  libmesh_assert(elem);
// return fe_lagrange_1D_shape_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// ===========================================================================
//       Shape derivative (2nd)
// ===========================================================================
// ================================================================================================ 
template <>
double  MGFE<1,LAGRANGE>::shape_second_deriv(
  const ElemType,
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double  p[]
) {// ================================================================================================
  return fe_lagrange_1D_shape_second_deriv(order, i, j, p[0]);
 }
//  
// template <> Real  MGFE<1,L2_LAGRANGE>::shape_second_deriv(const ElemType, const Order order,   const unsigned int i,
//                                            const unsigned int j, const Point & p){
//   return fe_lagrange_1D_shape_second_deriv(order, i, j, p(0)); }

// template <>
// Real  MGFE<1,LAGRANGE>::shape_second_deriv(const Elem * elem,    const Order order,
//                                         const unsigned int i,   const unsigned int j,
//                                         const Point & p,  const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,L2_LAGRANGE>::shape_second_deriv(const Elem * elem,     const Order order,
//                                            const unsigned int i,    const unsigned int j, const Point & p,
//                                            const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,LAGRANGE>::shape_second_deriv(const FEType fet,   const Elem * elem,      const unsigned int i,
//                                         const unsigned int j,           const Point & p,
//                                         const bool add_p_level) {
//   libmesh_assert(elem);   return fe_lagrange_1D_shape_second_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0));
// }

// template <>
// Real  MGFE<1,L2_LAGRANGE>::shape_second_deriv(const FEType fet,     const Elem * elem,  const unsigned int i,
//                                            const unsigned int j,   const Point & p,  const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// } // namespace libMesh
// 



// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // Basic QUAD/HEX
// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Edge_Quad_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = (1 - 0.5 * fabs(PhiCoeff)) * ((2 * fabs(PhiCoeff) - 1) * Coordinate * Coordinate +
//       PhiCoeff * Coordinate + (1 - fabs(PhiCoeff)));
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Edge_Quad_DPhi(int PhiCoeff, double Coordinate) {
//   double DPhiVal;
//   DPhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1) * Coordinate + PhiCoeff);
//   return DPhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Edge_Quad_D2Phi(int PhiCoeff, double /*Coordinate*/) {
//   double D2PhiVal;
//   D2PhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1));
//   return D2PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Edge_Lin_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = 0.5 * (1 + PhiCoeff * Coordinate);
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Edge_Lin_DPhi(int PhiCoeff, double /*Coordinate*/) {
//   double DPhiVal;
//   DPhiVal = 0.5 * PhiCoeff;
//   return DPhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Lin_Phi(int nPhi, double point[], int dimension) {
//   double PhiVal = 1.;
// 
//   for(int dir = 0; dir < dimension; dir++) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + dir * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Lin_Phi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiVal *= Edge_Lin_Phi(lambda_i, point[dir]);
//   }
// 
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Quad_Phi(int nPhi, double point[], int dimension) {
//   double PhiVal = 1.;
// 
//   for(int dir = 0; dir < dimension; dir++) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + dir * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Quad_Phi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiVal *= Edge_Quad_Phi(lambda_i, point[dir]);
//   }
// 
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer) {
//   int lambda_i;
// 
//   if(dimension == 1) {
//     lambda_i = _CooE3[nPhi];
//   } else if(dimension == 2) {
//     lambda_i = _CooQ9[nPhi + DirDer * _Q9Off];
//   } else if(dimension == 3) {
//     lambda_i = _CooH27[nPhi + DirDer * _H27Off];
//   } else {
//     printf("\033[1;31m MGFE1::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//   }
// 
//   double PhiDer = Edge_Lin_DPhi(lambda_i, point[DirDer]);
// 
//   for(int dir = DirDer + 1; dir < dimension + DirDer; dir++) {
//     int lambda_i;
//     int direction = dir % dimension;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + direction * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiDer *= Edge_Lin_Phi(lambda_i, point[direction]);
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer) {
//   int lambda_i;
// 
//   if(dimension == 1) {
//     lambda_i = _CooE3[nPhi];
//   } else if(dimension == 2) {
//     lambda_i = _CooQ9[nPhi + DirDer * _Q9Off];
//   } else if(dimension == 3) {
//     lambda_i = _CooH27[nPhi + DirDer * _H27Off];
//   } else {
//     printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//   }
// 
//   double PhiDer = Edge_Quad_DPhi(lambda_i, point[DirDer]);
// 
//   for(int dir = DirDer + 1; dir < dimension + DirDer; dir++) {
//     int lambda_i;
//     int direction = dir % dimension;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + direction * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiDer *= Edge_Quad_Phi(lambda_i, point[direction]);
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
//   double PhiDer;
// 
//   if(DirDer1 == DirDer2) {
//     PhiDer = 0.;
//   } else {
//     PhiDer = 1;
// 
//     for(int dir = 0; dir < dimension; dir++) {
//       int lambda_i;
//       int direction = dir % dimension;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + direction * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Lin_D2Phi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       if(dir == DirDer1 || dir == DirDer2) {
//         PhiDer *= Edge_Lin_DPhi(lambda_i, point[dir]);
//       } else {
//         PhiDer *= Edge_Lin_Phi(lambda_i, point[dir]);
//       }
//     }
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
//   double PhiDer;
// 
//   if(DirDer1 == DirDer2) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + DirDer1 * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + DirDer1 * _H27Off];
//     }
// 
//     PhiDer = Edge_Quad_D2Phi(lambda_i, point[DirDer1]);
// 
//     for(int dir = DirDer1 + 1; dir < dimension + DirDer1; dir++) {
//       int lambda_i;
//       int direction = dir % dimension;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + direction * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       PhiDer *= Edge_Quad_Phi(lambda_i, point[direction]);
//     }
//   } else {
//     PhiDer = 1;
// 
//     for(int dir = 0; dir < dimension; dir++) {
//       int lambda_i;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + dir * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Quad_D2Phi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       if(dir == DirDer1 || dir == DirDer2) {
//         PhiDer *= Edge_Quad_DPhi(lambda_i, point[dir]);
//       } else {
//         PhiDer *= Edge_Quad_Phi(lambda_i, point[dir]);
//       }
//     }
//   }
// 
//   return PhiDer;
// }
// 
// 
// 
// 
// 
// 
// 
// 
// //   // gaussian coordinates
// //   const double a = -sqrt(3. / 5.);
// //   const double b = 0.;
// //   const double c = -a;
// //
// //   const double x1D[3] = {a, b, c};
// //
// //   const double x2D[9] = {a, a, a, b, b, b, c, c, c};
// //   const double y2D[9] = {a, b, c, a, b, c, a, b, c};
// //
// //   const double x3D[27] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c};
// //   const double y3D[27] = {a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c};
// //   const double z3D[27] = {a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c};
// //
// //   // gaussian weights
// //   // 1D --------------------------------
// //   double weight1[3];
// //   weight1[0] = 5. / 9.;
// //   weight1[1] = 8. / 9.;
// //   weight1[2] = 5. / 9.;
// //
// //
// //
// //   // ****LINEAR SHAPES AND DERIVATIVES****
// //   // shape 2D                                            derivatives
// //   //                                           |_______x_______|_______y_______
// //   // phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// //   // phi1= x*(1.-y)                            |     1.-y      |      -x
// //   // phi2= x*y                                 |       y       |       x
// //   // phi3= y*(1.-x)                            |      -y       |     1.-x
// //
// //   // 1D --------------------------------
// //
// //   if(_dim==1) {
// //     for(int i = 0; i < _NoGauss; i++) { _weight[i] = weight1[i]; }
// //
// //     for(int i = 0; i < _NoGauss; i++) {
// //       phi[i] = 1.;
// //       _dphidxez_map[i] = 0.;
// //     }
// //   }
// //
// //   // 2D --------------------------------
// //   if(_dim==2) {
// //
// //     const double m = weight1[0] * weight1[1];
// //     const double l = weight1[0] * weight1[0];
// //     const double h = weight1[1] * weight1[1];
// //     const double weight_2[9] = {l, m, l, m, h, m, l, m, l};
// //
// //     for(int i = 0; i < _NoGauss; i++) { _weight[i] = weight_2[i]; }
// //
// //
// //     for(int i = 0; i < _NoGauss; i++) {
// //       const double xx = x2D[i];
// //       const double yy = y2D[i];
// //       // shape functions
// //       phi[i + 0 * _NoGauss] = 1.;
// //       // derivatives
// //       _dphidxez_map[i + (0) * _NoGauss] = 0.;  // d/dx
// //       _dphidxez_map[i + (1) * _NoGauss] = 0.;  // d/dy
// //
// //
// //     }
// //   }
// //
// //   // 3D -----------------------------------------------
// //   if(_dim == 3) {
// //
// //
// //     const double w1 = weight1[0] * weight1[0] * weight1[0];
// //     const double w2 = weight1[0] * weight1[0] * weight1[1];
// //     const double w3 = weight1[0] * weight1[1] * weight1[1];
// //     const double w4 = weight1[1] * weight1[1] * weight1[1];
// //     const double weight_3[27] = {w1, w2, w1, w2, w3, w2, w1, w2, w1, w2, w3, w2, w3, w4,
// //                                  w3, w2, w3, w2, w1, w2, w1, w2, w3, w2, w1, w2, w1
// //                                 };
// //
// //     for(int i = 0; i < _NoGauss; i++) { _weight[i] = weight_3[i]; }
// //
// //
// //     for(int i = 0; i < _NoGauss; i++) {
// //       const double xx = x3D[i];
// //       const double yy = y3D[i];
// //       const double zz = z3D[i];
// //
// //       // shape functions
// //       phi[i + (0) * _NoGauss] = 1.;
// //       // derivatives
// //       _dphidxez_map[i + (0) * _NoGauss] = 0.;  // d/dx
// //       _dphidxez_map[i + (1) * _NoGauss] = 0.;  // d/dy
// //       _dphidxez_map[i + (2) * _NoGauss] = 0.;  // d/dz
// //
// //
// //     }
// //   }  // end if 3D
// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// // Basic TRI/TETRA
// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// 
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_2d_LinearPhi(
//   int nPhi,
//   double point[]
// ) {// ====================================================================================
//   double PhiVal;
//   int a = _CooTriEl[nPhi * 3];
//   int b = _CooTriEl[nPhi * 3 + 1];
//   int c = _CooTriEl[nPhi * 3 + 2];
//   PhiVal = a * point[0] + b * point[1] + c;
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_2d_QuadraticPhi(
//   int nPhi,
//   double point[]
// ) {// ====================================================================================
//   double PhiVal;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double L1, L2;
//   L1 = Tri_2d_LinearPhi(l1, point);
//   L2 = Tri_2d_LinearPhi(l2, point);
//   PhiVal = L1 * ((1 - m) * (2. * L1 - 1) + m * 4. * L2);
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_2d_LinearDerPhi(int nPhi, double */*point[]*/, int dir) {
//   double PhiDer;
//   PhiDer = _CooTriEl[nPhi * 3 + dir];
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir) {
//   double PhiDer;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double L1, L2, dL1, dL2;
//   L1 = Tri_2d_LinearPhi(l1, point);
//   L2 = Tri_2d_LinearPhi(l2, point);
//   dL1 = Tri_2d_LinearDerPhi(l1, point, dir);
//   dL2 = Tri_2d_LinearDerPhi(l2, point, dir);
// 
//   PhiDer = dL1 * ((1 - m) * (2. * L1 - 1) + m * 4. * L2) + L1 * ((1 - m) * (2. * dL1) + m * 4. * dL2);
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2) {
//   double Phi2Der;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double dL1_2, dL2_2, dL1_1, dL2_1;
// 
//   dL1_1 = Tri_2d_LinearDerPhi(l1, point, dir1);
//   dL2_1 = Tri_2d_LinearDerPhi(l2, point, dir1);
//   dL1_2 = Tri_2d_LinearDerPhi(l1, point, dir2);
//   dL2_2 = Tri_2d_LinearDerPhi(l2, point, dir2);
// 
//   Phi2Der =
//     dL1_1 * ((1 - m) * (2. * dL1_2) + m * 4. * dL2_2) + dL1_2 * ((1 - m) * (2. * dL1_1) + m * 4. * dL2_1);
//   return Phi2Der;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_3d_LinearPhi(int nPhi, double point[]) {
//   double PhiVal;
//   int a = _CooTetra4[nPhi * 4];
//   int b = _CooTetra4[nPhi * 4 + 1];
//   int c = _CooTetra4[nPhi * 4 + 2];
//   int d = _CooTetra4[nPhi * 4 + 3];
//   PhiVal = a * point[0] + b * point[1] + c * point[2] + d;
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_3d_LinearDerPhi(int nPhi, double */*point[]*/, int dir) {
//   double PhiDer;
//   PhiDer = _CooTetra4[nPhi * 4 + dir];
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_3d_QuadraticPhi(int nPhi, double point[]) {
//   double PhiVal;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   int b = _CooTetra10[nPhi * 4 + 3];
//   double L1, L2;
//   L1 = Tri_3d_LinearPhi(l1, point);
//   L2 = Tri_3d_LinearPhi(l2, point);
//   PhiVal = L1 * (a * L2 + b);
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_3d_QuadraticDerPhi(int nPhi, double point[], int dir) {
//   double PhiDer;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   int b = _CooTetra10[nPhi * 4 + 3];
//   double L1, L2, dL1, dL2;
//   L1 = Tri_3d_LinearPhi(l1, point);
//   L2 = Tri_3d_LinearPhi(l2, point);
//   dL1 = Tri_3d_LinearDerPhi(l1, point, dir);
//   dL2 = Tri_3d_LinearDerPhi(l2, point, dir);
//   PhiDer = dL1 * (a * L2 + b) + L1 * a * dL2;
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<1,LAGRANGE>::Tri_3d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2) {
//   double PhiDer;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   double dL1_1, dL2_1, dL1_2, dL2_2;
//   dL1_1 = Tri_3d_LinearDerPhi(l1, point, dir1);
//   dL2_1 = Tri_3d_LinearDerPhi(l2, point, dir1);
//   dL1_2 = Tri_3d_LinearDerPhi(l1, point, dir2);
//   dL2_2 = Tri_3d_LinearDerPhi(l2, point, dir2);
//   PhiDer = a * (dL1_1 * dL2_2 + dL1_2 * dL2_1);
//   return PhiDer;
// }
// 
// template <> double  MGFE<1,LAGRANGE>::FirstDerivateOfLocalPhi(int nPhi, double point[], int Dimension, int DirDer, int FamilyType) {
//   double PhiDer;
// 
//   if(FamilyType == 0) {   // TRIANGULAR ELEMENTS
//     if(_order == 1) {
//       if(Dimension == 1) {
//         PhiDer = Rec_Lin_DPhi(nPhi, point, Dimension, DirDer);
//       } else if(Dimension == 2) {
//         PhiDer = Tri_2d_LinearDerPhi(nPhi, point, DirDer);
//       } else if(Dimension == 3) {
//         PhiDer = Tri_3d_LinearDerPhi(nPhi, point, DirDer);
//       }
//     }
// 
//     if(_order == 2) {
//       if(Dimension == 1) {
//         PhiDer = Rec_Quad_DPhi(nPhi, point, Dimension, DirDer);
//       } else if(Dimension == 2) {
//         PhiDer = Tri_2d_QuadraticDerPhi(nPhi, point, DirDer);
//       } else if(Dimension == 3) {
//         PhiDer = Tri_3d_QuadraticDerPhi(nPhi, point, DirDer);
//       }
//     }
//   }
// 
//   if(FamilyType == 1) {   // QUADRANGULAR ELEMENTS
//     if(_order == 1) {
//       PhiDer = Rec_Lin_DPhi(nPhi, point, Dimension, DirDer);
//     }
// 
//     if(_order == 2) {
//       PhiDer = Rec_Quad_DPhi(nPhi, point, Dimension, DirDer);
//     }
//   }
// 
//   return PhiDer;
// }

// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic piecewise continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// 
// /// /// This function generates the Lagrangian piecewise shape functions
// template <> void MGFE<1,LAGRANGE>::init_pie_rec() {  // ================================
// 
// 
// //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
//   //       1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   int hex_shift=6;  //  EDGE3  _CooE3[i]=[_CooH27[i+6]
// 
// 
//   // ---------------------------------------------------------------------------------------
//   _fem_conn[0]=_NoElNodes-1;
// 
//   // dofs  ---------------------------------------------------------------------------
//   for(int j = 0; j < _NoElNodes; j++)
//     for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
// 
//   // ---------------------------------------------------------------------------------------
// 
//   // nodes ----------------------------------------------------------------------------
//   double *pointNode=new double[_dim];
// 
//   for(int i = 0; i < _NoElNodes; i++) {  // points
//     for(int dir = 0; dir < _dim; dir++)  pointNode[dir] =   _xnodes[i+_NoElNodes*dir];
// 
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]   =1.;
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =0;     //
//       }
//     }
//   }
// 
//   delete[]pointNode;
//   // ---------------------------------------------------------------------------------------
//  // gaussian points --------------------------------------------------------->
//  ElemType  elem_name=_fe_type.fetype; 
//   Order elem_order= _fe_type.order ;
//   double *point=new double[1];
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv(elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;   
// //   // Gauss integration points (rectangular coordinates) ---------------------------------------------------
// //   // basic gaussian points coordinates
// //   const double a = -sqrt(3. / 5.);
// //   const double b = 0.;
// //   const double c = -a;
// //   const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
// //       a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
// //       a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
// //     };
// //   // gaussian weights 1D
// //   double weight1[3];
// //   weight1[0] = 5. / 9.;
// //   weight1[1] = 8. / 9.;
// //   weight1[2] = 5. / 9.;
// // 
// //   double *point=new double[_dim];
// //   int iw[3];
// // 
// //   for(int i = 0; i < _NoGauss; i++) {
// //     iw[0]=i%3;
// //     iw[1]= (i/3)%3;
// //     iw[2]= i/9;
// //     _weight[i] =1.;
// // 
// //     for(int dir=0; dir<_dim; dir++) {
// //       _weight[i] *= weight1[iw[dir]];        // gausssian weights
// //       _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
// //       point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
// //     }
// //     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
// //       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
// //       phi[j*_NoGauss+i] = 1.; // shape functions
// //       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
// //         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = 0. ;
// //         for(int dir2 = 0; dir2 < _dim; dir2++)
// //           _dphidxx_map[(j + (dir * _dim + dir2) *  _NoShape) * _NoGauss + i] =0.;
// // 
// // //           Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
// //       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// //     }
// //   }
// // 
// //   delete[]point;
// //   //--------------------------------------------------------------------------------
// // //--------------------------------------------------------------------------------
// //     // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
// //    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
// //   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
// //     std::string name =   "const_1D";//+_var_names_D[ivar];
// //     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
// //      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
// //   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
// //                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
// //   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
// //   H5Sclose(dataspace);
// //   H5Dclose(dataset);
// //     
// //       H5Fclose(file);
// //   
// //   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh  
// 
//   return;
// }
// // ===================================================================================
// /// /// This function generates the Lagrangian piecewise shape functions
// template <> void MGFE<1,LAGRANGE>::init_pie_tri(
// ) {  // ================================================================================
// 
//   /* //        ********************************************
//     //                               TRI 3
//     //     ********************************************
//     //                              2
//     //        |\
//     //        | \
//     //        |  \
//     //   r=l1 |   \ 1-s-r=l0
//     //        |    \
//     //        |     \
//     //        |______\
//     //        0  s=l2 1
//     */
//   // gaussian coordinates
//   // 1D --------------------------------
//   const double x[3] = {-sqrt(3. / 5.), 0., sqrt(3. / 5.)};
//   // 2D --------------------------------
//   const double alphabeta[4] = {1. - 2. * (2. / 7. + sqrt(15.) / 21.), 1. - 2. * (2. / 7. - sqrt(15.) / 21.),
//       2. / 7. + sqrt(15.) / 21., 2. / 7. - sqrt(15.) / 21.
//     };
//   const double r[7] = {1. / 3.,      alphabeta[2], alphabeta[0], alphabeta[2],
//       alphabeta[3], alphabeta[1], alphabeta[3]
//     };
//   const double s[7] = {1. / 3.,      alphabeta[2], alphabeta[2], alphabeta[0],
//       alphabeta[3], alphabeta[3], alphabeta[1]
//     };
// 
// //   int n_nodes=NDOF_FEM;
//   /*  if(_dim==1) {n_nodes=3;} //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//     if(_dim==2) {n_nodes=NDOF_FEMB;}//  QUAD9  _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)] */
//   const double x10[11*3] = {0,1,0,0,  .5,.5,0.,0.,0.,.5,  0.33333333333333333333333,
//       0,0,1,0,  0.,.5,.5,0.,.5,0.,  0.33333333333333333333333,
//       0,0,0,1,  0.,0.,0.,.5,.5,.5,  0.33333333333333333333333
//     };
// // 10 point = quad tetra  + center point
//   // ****LINEAR SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates              derivatives
//   //                                      |_______r_______|_______s_______
//   // phi0= 1-r-s                          |      -1.      |      -1.
//   // phi1= r                              |       1.      |       0.
//   // phi2= s                              |       0.      |       1.
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
//     // dofs -----------------------------------------------------------------------
//     _fem_conn[0]=2;
//     // ---------------------------------------------------------------------------------------
//     // nodes ----------------------------------------------------------------------------
//     _xnodes[0]=0.;    _xnodes[1]=1.;   _xnodes[2]=0.5;
// 
//     double *pointNode=new double[_dim];
//     for(int i = 0; i < _NoElNodes; i++) {  // points
//       for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
//       for(int j = 0; j < _NoShape; j++) {  // shapse
//         //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//         phi_nodes[j*_NoElNodes+i]   = 1.;
//         // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//         for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =0.;     // dphi/dxi
//         }
//       }
//     }
//     delete[]pointNode;
//     
//     
//  // gaussian points --------------------------------------------------------->
//  ElemType  elem_name=_fe_type.fetype; 
//   Order elem_order= _fe_type.order ;
//   double *point=new double[1];
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv(elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;            
// // // Gaussian points ---------------------------------------------------------------
// //     // gaussian weights
// //     _weight[0] = 5. / 18.;    _weight[1] = 8. / 18.;    _weight[2] = 5. / 18.;
// //     for(int i = 0; i < _NoGauss; i++) {  _xg[i] = 0.5*(x[i]+1);      phi[i] =1. ;      _dphidxez_map[i] = 0.; }
// //   
// //     // delete[] point;
// // //--------------------------------------------------------------------------------
// //     // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
//    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
//     std::string name =   "const_1D";//+_var_names_D[ivar];
//     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
//      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
//   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
//                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
//   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//     
//       H5Fclose(file);
//   
//   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh  
//   }
// 
// // ================================================================================================
// /// This function generates the Lagrangian quad shape functions
// template <> void MGFE<1,LAGRANGE>::init_pie() {  // ================================
// 
//   if(_GeomType==1) init_pie_rec();
//   if(_GeomType==0) init_pie_tri();
//   return;
// }

// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic piecewise linear functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// template <> void MGFE<1,LAGRANGE>::init_lin_rec(
// 
// ) {// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// 
//   //    ********************************************
//   //                               Edge 2
//   //     ********************************************
//   //
//   //              ___________|     o---> xi
//   //             0             1
//   //   (xi) is the reference element coordinate associated with
//   //   the given numbering.
//   //   // ****LINEAR SHAPES AND DERIVATIVES****
//   //   // shape 2D                                            derivatives
//   //   //                                           |_______x_______|_______y_______
//   //   // phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
//   //   // phi1= x*(1.-y)                            |     1.-y      |      -x
//   //   // phi2= x*y                                 |       y       |       x
//   //   // phi3= y*(1.-x)                            |      -y       |     1.-x
//   //
//   //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1,
//   //     -1, 1,  1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   int hex_shift= 4;   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   
// 
//   // ---------------------------------------------------------------------------------------
// 
//   // shape order  ---------------------------------------------------------------------------
//   for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
//   // nodes ----------------------------------------------------------------------------
//  _xnodes[0]=0.;    _xnodes[1]=1.;   _xnodes[2]=0.5;// central point
//   // // ---------------------------------------------------------------------------------------
//   // 
//   // // nodes ----------------------------------------------------------------------------
//   // for(int j = 0; j < _NoElNodes; j++)for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
//   // for(int dir = 0; dir < _dim; dir++) _xnodes[_NoElNodes-1+_NoElNodes*dir]=0.;
// 
//   // 1D function and derivatives at points(3)
//   double *pointNode=new double[1];
//   for(int i = 0; i < _NoElNodes; i++) {  // points
//     pointNode[0] = _xnodes[i];
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]   = Rec_Lin_Phi(j, pointNode,_dim);
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = Rec_Lin_DPhi(j,pointNode,_dim,dir);     // dphi/dxi
//       }
//     }
//   }
// 
//   delete[]pointNode; // ---------------------------------------------------------------------------------------
// 
// // gaussian points --------------------------------------------------------->
//  ElemType  elem_name=_fe_type.fetype; 
//   Order elem_order= _fe_type.order ;
//   double *point=new double[1];
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv(elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;       
//    phi=phi;
//   // // Gauss integration points (rectangular coordinates) ---------------------------------------------------
//   // // basic gaussian points coordinates
//   // const double a = -std::sqrt(3. / 5.);
//   // const double b = 0.0;
//   // const double c = -a;
//   // const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
//   //     a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
//   //     a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
//   //   };
//   // // gaussian weights 1D
//   // double weight1[3];
//   // weight1[0] = 5. / 9.;
//   // weight1[1] = 8. / 9.;
//   // weight1[2] = 5. / 9.;
//   // 
//   // double *point=new double[_dim];
//   // int iw[3];
//   // 
//   // for(int i = 0; i < _NoGauss; i++) {
//   //   iw[0]=i%3;    iw[1]= (i/3)%3;    iw[2]= i/9;    _weight[i] =1.;
//   // 
//   //   for(int dir=0; dir<_dim; dir++) {
//   //     _weight[i] *= weight1[iw[dir]];        // gausssian weights
//   //     _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
//   //     point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
//   //   }
//   // 
//   //   for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//   //     // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
//   //     phi[j*_NoGauss+i] = fe_lagrange_1D_shape(FIRST,j,point[0]);
//   //     // Rec_Lin_Phi(j, point, _dim); // shape functions
//   // 
//   //     for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//   //       _dphidxez_map[(j+dir*_NoShape)*_NoGauss+i] = fe_lagrange_1D_shape_deriv(FIRST,j,dir, point[0]);
//   //       // Rec_Lin_DPhi(j, point, _dim, dir);
//   //       for(int dir2 = 0; dir2 < _dim; dir2++) // 2nd derivatives
//   //         _dphidxx_map[(j+(dir*_dim+dir2)*_NoShape)*_NoGauss+i] =
//   //          fe_lagrange_1D_shape_second_deriv(FIRST,j,dir +dir2, point[0]);
//   //         // =0.;
//   //     }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   //   }
//   // }
//   // 
//   // delete[]point;
//   // //--------------------------------------------------------------------------------
//   //   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
// //    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
// //   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
// //     std::string name =   "lin_1D";//+_var_names_D[ivar];
// //     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
// //      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
// //   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
// //                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
// //   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
// //   H5Sclose(dataspace);
// //   H5Dclose(dataset);
// //     
// //       H5Fclose(file);
//   
//   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh     
//   
//  
// 
//   return;
// }
// 
// 
// template <> void MGFE<1,LAGRANGE>::init_lin_rec_rt() {  return;}
// 
// 
// // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// template <> void MGFE<1,LAGRANGE>::init_lin_tri(
// ) {// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
//   /* //        ********************************************
//   //                               TRI 3
//   //     ********************************************
//   //                              2
//   //        |\
//   //        | \
//   //        |  \
//   //   r=l1 |   \ 1-s-r=l0
//   //        |    \
//   //        |     \
//   //        |______\
//   //        0  s=l2 1
//   */
//   // gaussian coordinates
//   // 1D --------------------------------
//   const double x[3] = {-sqrt(3. / 5.), 0., sqrt(3. / 5.)};
//   // // 2D --------------------------------
//   // const double alphabeta[4] = {1. - 2. * (2. / 7. + sqrt(15.) / 21.), 1. - 2. * (2. / 7. - sqrt(15.) / 21.),
//   //     2. / 7. + sqrt(15.) / 21., 2. / 7. - sqrt(15.) / 21.
//   //   };
//   // const double r[7] = {1. / 3.,      alphabeta[2], alphabeta[0], alphabeta[2],
//   //     alphabeta[3], alphabeta[1], alphabeta[3]
//   //   };
//   // const double s[7] = {1. / 3.,      alphabeta[2], alphabeta[2], alphabeta[0],
//   //     alphabeta[3], alphabeta[3], alphabeta[1]
//   //   };
//   // // 3D --------------------------------
//   // const double p[14] = {0.31088591926330060980,
//   //     0.31088591926330060980,
//   //     1. - 3. * 0.31088591926330060980,
//   //     0.31088591926330060980,
//   //     0.092735250310891226402,
//   //     0.092735250310891226402,
//   //     1. - 3. * 0.092735250310891226402,
//   //     0.092735250310891226402,
//   //     0.5 - 0.045503704125649649492,
//   //     0.5 - 0.045503704125649649492,
//   //     0.045503704125649649492,
//   //     0.045503704125649649492,
//   //     0.5 - 0.045503704125649649492,
//   //     0.045503704125649649492
//   //   };
//   // const double q[14] = {
//   //   0.31088591926330060980,  1. - 3. * 0.31088591926330060980, 0.31088591926330060980,
//   //   0.31088591926330060980,  0.092735250310891226402,          1. - 3. * 0.092735250310891226402,
//   //   0.092735250310891226402, 0.092735250310891226402,          0.5 - 0.045503704125649649492,
//   //   0.045503704125649649492, 0.045503704125649649492,          0.5 - 0.045503704125649649492,
//   //   0.045503704125649649492, 0.5 - 0.045503704125649649492
//   // };
//   // const double t[14] = {0.31088591926330060980,        0.31088591926330060980,
//   //     0.31088591926330060980,        1. - 3. * 0.31088591926330060980,
//   //     0.092735250310891226402,       0.092735250310891226402,
//   //     0.092735250310891226402,       1. - 3. * 0.092735250310891226402,
//   //     0.045503704125649649492,       0.045503704125649649492,
//   //     0.5 - 0.045503704125649649492, 0.045503704125649649492,
//   //     0.5 - 0.045503704125649649492, 0.5 - 0.045503704125649649492
//   //   };
// 
// 
// //   int n_nodes=NDOF_FEM;
//   /*  if(_dim==1) {n_nodes=3;} //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//     if(_dim==2) {n_nodes=NDOF_FEMB;}//  QUAD9  _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)] */
//   const double x10[11*3] = {0,1,0,0,  .5,.5,0.,0.,0.,.5,  0.33333333333333333333333,
//       0,0,1,0,  0.,.5,.5,0.,.5,0.,  0.33333333333333333333333,
//       0,0,0,1,  0.,0.,0.,.5,.5,.5,  0.33333333333333333333333
//     };
// // 10 point = quad tetra  + center point
// 
// 
//   // ****LINEAR SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates              derivatives
//   //                                      |_______r_______|_______s_______
//   // phi0= 1-r-s                          |      -1.      |      -1.
//   // phi1= r                              |       1.      |       0.
//   // phi2= s                              |       0.      |       1.
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
//     ElemType  elem_name=_fe_type.fetype; 
//     Order elem_order= _fe_type.order ;
//     // dofs -----------------------------------------------------------------------
//     for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
//     double *xnodes=new double [_NoElNodes*_dim];
//     double *point=new double[_dim];
//     // ---------------------------------------------------------------------------------------
//     // nodes ----------------------------------------------------------------------------
//     
//  
//  
//   
//   get_refspace_nodes(elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
//     // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
//     for(int i = 0; i < _NoElNodes; i++) {  // points
//       for(int dir = 0; dir < _dim; dir++) { 
//         _xnodes[i+_NoElNodes*dir]=xnodes[dir+_dim*_NoElNodes]; 
//         point[dir] = _xnodes[i+_NoElNodes*dir];
//       }
// 
//       for(int j = 0; j < _NoShape; j++) {  // shapse
//         //  phi_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//         phi_nodes[j*_NoElNodes+i]   =  fe_lagrange_1D_shape(elem_order, j,*point);//Rec_Quad_Phi(j, pointNode,1);
//         // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//         for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = fe_lagrange_1D_shape_deriv(elem_order, j,dir,* point);//Rec_Quad_DPhi(j,pointNode,1,dir);     // dphi/dxi
//         }
//       }
//     }
// 
//     delete[] xnodes;
// 
// // gaussian points --------------------------------------------------------->
// 
//  
// 
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv(elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;          
// // // Gaussian points ---------------------------------------------------------------
// //     // gaussian weights
// //     _weight[0] = 5. / 18.;    _weight[1] = 8. / 18.;    _weight[2] = 5. / 18.;
// // 
// // 
// //     for(int i = 0; i < _NoGauss; i++) {
// //       _xg [i] = x[i];
// //       phi[i] = 0.5*(1.-_xg[i]);
// //       phi[i + _NoGauss] =0.5*(1. +_xg[i]);
// //       _dphidxez_map[i] = -1.;
// //       _dphidxez_map[i + _NoGauss] = 1.;
// //     }
// //     
// //   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
//    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
//     std::string name =   "lin_1D";//+_var_names_D[ivar];
//     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
//      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
//   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
//                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
//   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//     
//       H5Fclose(file);
//   
//   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh     
//     
//     
//   return;
// }



// // ================================================================================================
// ///  This function generates the Lagrangian linear shape functions
// template <> void MGFE<1,LAGRANGE>::init_lin() {  // ================================
// 
//   switch(_GeomType) {
//   case 1: // quad/hex ---------------------
//     if(_FamType==0)  init_lin_rec();
//     if(_FamType==1){      init_lin_rec_rt();    }
//     if(_FamType>1) {      std::cout<<" MGFE1::init_lin error: _FamType";      exit(1);    }   break;
// 
//   case 0:   // tri/tetra ---------------------
//     if(_FamType==0)  init_lin_tri();
//     else {      std::cout<<" MGFE1::init_lin error: _FamType";     exit(3);    }    break;
//   default:     std::cout<<" MGFE1::init_lin error: _GeomType";    exit(3);
//   }
//   return;
// }
// 
// 
// 
// // ================================================================================================
// //    QUAD
// // ============================================================================
// template <> void MGFE<1,LAGRANGE>::init_qua_rec(
// ) {//==========================================================================
// 
// ///
// //     ********************************************
// //                               EDGE 3
// //     ********************************************
// ///
// ///  *        |____________|     o---> xi
// ///  *       0       2      1
//   /// * (xi) is the reference element coordinate associated with
// /// * the given numbering.
// ///
//  
// // dofs  ---------------------------------------------------------------------------
//   for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
//   // nodes ----------------------------------------------------------------------------
//   _xnodes[0]= -1; _xnodes[1]= 1;_xnodes[2]= 0.;
//   
//   double *pointNode=new double[1];
//   for(int i = 0; i <  _NoElNodes; i++) {  // points
//     pointNode[0] = _xnodes[i];
//    
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       phi_nodes[j*_NoElNodes+i]   = Rec_Quad_Phi(j, pointNode,_dim);
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = Rec_Quad_DPhi(j,pointNode,_dim,dir);     // dphi/dxi
//       }
//     }
//   }
//   delete[]pointNode;
// 
// // gaussian points --------------------------------------------------------->
//  ElemType  elem_name=_fe_type.fetype; 
//   Order elem_order= _fe_type.order ;
//   double *point=new double[1];
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape(elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv(elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;      
//   
//    phi=phi;
//   
//   
// // // basic gaussian points coordinates  ---------------------------------------------------------------------
// //    _xg[0] =-sqrt(3. / 5.);  _xg[2]  =0.; _xg[1] =-_xg[0]; // gaussian points
// //    _weight[0]=5. / 9.; _weight[2] = 8. / 9.;  _weight[1]= _weight[0];// gaussian weights 1D
// //   
// //  
// //  double *point=new double[1];
// //   for(int i = 0; i < _NoGauss; i++) {
// //     
// //       point[i] =_xg[i];  //gaussian points
// // 
// //     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
// //       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
// //       phi[j*_NoGauss+i] = Rec_Quad_Phi(j, point, _dim); // shape functions
// // 
// //       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
// //         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = Rec_Quad_DPhi(j, point, _dim, dir);
// // 
// //         for(int dir2 = 0; dir2 < _dim; dir2++)
// //           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
// //             Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
// //       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// //     }
// //   } //--------------------------------------------------------------------------------
//   
//    // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
// //    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
// //   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
// //     std::string name =   "quad_1D";//+_var_names_D[ivar];
// //     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
// //      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
// //   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
// //                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
// //   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
// //   H5Sclose(dataspace);
// //   H5Dclose(dataset);
// //     
// //       H5Fclose(file);
//   
//   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh  
//  
//   return;
// }
// 
// 
// // ================================================================================
// template <> void MGFE<1,LAGRANGE>::init_qua_tri() {
//   /*           ********************************************
//   //                               TRI 6
//   //     ********************************************
//   //                              2
//   //        |\
//   //        | \
//   //        |  \
//   //         r=l1 |   \ 1-s-r=l0
//   //        |    \
//   //        |     \
//   //              |______\
//   //             0  s=l2 1
//   */
//   // gaussian coordinates
//   // 1D --------------------------------
//   const double x[3] = {-sqrt(3. / 5.), 0., sqrt(3. / 5.)};
// 
//  
//   // ****QUADRATIC SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates                  derivatives
//   //                                           |_______r_______|_______s_______
//   // phi0= 2.*(1-r-s)*((1.-r-s)-0.5)           | -3.+4.s+4.r   |  -3.+4.*s+4.*r
//   // phi1= 2.*r*(r-0.5)                        |    4.*r-1     |      0.
//   // phi2= 2.*s*(s-0.5)                        |    0.         |   4.*s-1
//   // phi3= 4.*r*(1.-r-s)                       |  4.-4.*s-8.*r |   -4.*r
//   // phi4= 4.*r*s                              |   4.*s        |   4.*r
//   // phi5= 4.*s*(1.-r-s)                       |  -4.*s        |   4-4.*r-8.*s
// 
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
// 
// 
//     // dofs -----------------------------------------------------------------------
//     for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
//     // ---------------------------------------------------------------------------------------
//     // nodes ----------------------------------------------------------------------------
//     _xnodes[0]=0.;    _xnodes[1]=1.;   _xnodes[2]=0.5;
//     double *pointNode=new double[1];
//     for(int i = 0; i < _NoElNodes; i++) {  // points
//       pointNode[i] = _xnodes[i];
// 
//       for(int j = 0; j < _NoShape; j++) {  // shapse
//         //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//         _phi_map_nodes[j*_NoElNodes+i]   =  Rec_Quad_Phi(j, pointNode,1);
//         // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//         for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = Rec_Quad_DPhi(j,pointNode,1,dir);     // dphi/dxi
//         }
//       }
//     }
// 
//     delete[]pointNode;
// 
// // gaussian points --------------------------------------------------------->
//  ElemType  elem_name=_fe_type.fetype; 
//   Order elem_order= _fe_type.order ;
//   double *point=new double[1];
//  
//  for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv( elem_order, j,dir, point[0]);
//         for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =
//           fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
// 
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } 
//   // <----------------------------------------------------------------------------
// 
//    delete[]point;          
//     
//   //   // Gauss ------------------------------------------------------------------------
//   //   // weights
//   //   _weight[0]=5./9.;    _weight[1]=8./9.;    _weight[2]=5./9.;
//   //   for(int i = 0; i <  _NoGauss; i++) {
//   //     _xg[i]=x[i];
//   //     // shape functions
//   //     for(int nPhi = 0; nPhi < _NoShape; nPhi++) {
//   //       int off = nPhi *  _NoGauss;
//   //       phi[i + off] = Rec_Quad_Phi(nPhi, _xg, _dim);
//   //       _dphidxez_map[i + off] = Rec_Quad_DPhi(nPhi, _xg, _dim, 0);
//   //       _dphidxx_map[i + off] = Rec_Quad_D2Phi(nPhi, _xg, _dim, 0, 0);
//   //     }
//   //   }
//   // // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
//   // ==================================================================
//    // file hdf5
//    // hid_t file = H5Fcreate("/home/ricky/software/numericplatform_v0.0.1/PLAT_USERS/ns/femus2_xfem/DATA/femq.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
// 
//    hid_t file = H5Fopen(femfileh5_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//   hsize_t dimsf[2];  dimsf[0] =  _NoShape*_NoGauss;  dimsf[1] = 1;
//     std::string name =   "quad_1D";//+_var_names_D[ivar];
//     // _mgutils_D.print_Ihdf5(file_id, var_name, dimsf, phi);
//      hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
//   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
//                             H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
//   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//     
//       H5Fclose(file);
//   
//   // hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh  
//   return;
// }
// 
// 
// // ================================================================================================
// /// This function generates the Lagrangian quad shape functions
// template <> void MGFE<1,LAGRANGE>::init_qua() {  // ================================
//  
//   if(_GeomType==1) init_qua_rec();
//   if(_GeomType==0) init_qua_tri();
//   return;
// }

// ================================================================================================
// ================================================================================================
// ================================================================================================


template <>   double MGFE<1,LAGRANGE>::JacSur(
  const int ng, 
  const double x[], 
  double InvJac[]
)  {// ========================================================
       double det=0.; double x_xi = 0.;
      int nshape = _NoShape;
       for(int  kdim = 0;  kdim < 2;  kdim++){ // FEM1dim ->JacSurf 2D -> 2
        x_xi = 0.;
      for(int s = 0; s < nshape; s++)  x_xi += x[s+ kdim *_NoGauss] * _dphidxez_map[ s * _NoGauss +   ng];
       det +=  x_xi* x_xi;
       }
      det= std::sqrt(det);
      InvJac[0] = 1./(det+1.e-20);  // deta dy

      return (det);
  }
  
  
  
// =======================================================
/// This function computes the  derivatives and Jacobian at the gaussian point ng
template <>   double MGFE<1,LAGRANGE>::Jac_nodes(
const int ng,          // Gaussian point <-
const double x[],            // coordinates  <-
double InvJac[]        // Jacobean
) {                        // ===================================================
  double x_xi = 0.;
  int nshape = _NoShape;
  //   int offset=_NoShape*_NoGauss;
  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    x_xi += x[s] * _dphidxez_map_nodes[sng];
  }
  double det = x_xi;
  double idet = 1. / det;
  InvJac[0] = x_xi * idet;  // deta dy
  return (det);
}





// ==============================================================
template<> double  MGFE<1,LAGRANGE>::JacOnGivenCanCoords(
  const int dim, 
  const double ElemCoords[], 
  double CanCoords[], 
  double InvJac[], 
  int FamilyType, 
  int nShape
) {// ===========================================================
//   double* LocDPhi = new double[dim];
//   double* CoordsDer = new double[dim * dim];
// 
//   for(int dir = 0; dir < dim * dim; dir++) { CoordsDer[dir] = 0.; }
// 
//   for(int s = 0; s < nShape; s++) {
//     for(int dir = 0; dir < dim; dir++) {
//       LocDPhi[dir] = FirstDerivateOfLocalPhi(s, CanCoords, dim, dir, FamilyType);
//     }
// 
//     for(int dir1 = 0; dir1 < dim; dir1++)
//       for(int dir2 = 0; dir2 < dim; dir2++) {
//         CoordsDer[dir1 * dim + dir2] += ElemCoords[s + dir1 * nShape] * LocDPhi[dir2];
//       }
//   }
// 
//   double det = CoordsDer[0];
// 
//  
// 
//   delete[] LocDPhi;
//   delete[] CoordsDer;
  // return det; 
  return 0.;
  
}

















// ======================================
/// This function computes the normal at the gauss point
template <>
void MGFE<1,LAGRANGE>::normal_g(
const double* xx,  // all surface coordinates <-
double* normal_g   // normal ->
) const {          // ======================================
  int NDOF_FEMB=1;
   // coordinates
  double xx3D[3 * NDOF_FEMB];
  normal_g[0]=1.;
 

  return;
}

// =========================================================================================
template <>
void  MGFE<1,LAGRANGE>::Oxy_face(
const double* xx,   ///< coordinates
const double x_c[], ///< central point
double* normal_g,   ///<  i_n,i_t1,i_t2 face Oxyz
double* tg01,
double* tg03
) const {            // ======================================
int NDOF_FEMB=1;
  // coordinates
  double xx3D[3 * NDOF_FEMB];
  normal_g[0]=1.;
  for(int i = 0; i < 1 * NDOF_FEMB; i++) { xx3D[i] = xx[i]; tg01[i] =0.;tg01[3] =0.;}
  // if the sign is not correct then reverse it
  if(normal_g[0] * (x_c[0] - xx3D[0])  >    0) {    normal_g[0] *= -1;    }

  return;
}

// =============================================================================
template <> void MGFE<1,LAGRANGE>::set_par(
  std::string fileh5,
  int dim_in,
  int order_in,
  int shapes_quad,
  int geom_type_in,
  int fam_type_in
) {// =============================================================================
  _dim=dim_in, _order= order_in,  
  _NoShape=shapes_quad,   _GeomType=geom_type_in,     _FamType=fam_type_in;
   femfileh5_=fileh5;
     ///< family type  0=Lagrangian(C^0) 1=Raviart Thomas (L^2-C^0)
//   // Fem order
//   // =============================
//   //  ------------------------------------------------------------------------
//   // gauss points (_NoGauss) ------------------------------------------------------
//     _NoGauss=0;
//   assert(_qrule);
//  _NoGauss=_qrule->n_points();
//  if(_NoGauss==0) {
//    std::cout << "Quadrature not assigned"; abort();
// }
//    // _NoGauss=3; 
//   ElemType  elem_name=_fe_type.fetype; 
//   
//   // number of nodes (_NoElNodes) in the reference EDGE3 element
//   _NoElNodes= get_refspace_n_nodes(elem_name);// number of nodes on element  ElemType (it should be from Geom) todo3;  
  // gauss points (_NoGauss) ------------------------------------------------------
   _NoGauss=0;  assert(_qrule); _NoGauss=_qrule->n_points();
 if(_NoGauss==0) {   std::cout << "Quadrature not assigned"; abort();}

  // number of nodes (_NoElNodes) in the reference quad element
   ElemType  elem_name=_fe_type.fetype; 
   ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
   Order elem_order= _fe_type.order ;
    _order= (int)elem_order;
  // number of nodes (_NoElNodes) in the reference EDGE3 element
    _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
   
  
  //  -----------------------------------------------------------------------------------
  // x node  coordinates   (_xnodes[xp1,xp2,xp3,...yp1,yp2,yp3,....zp1,zp2,zp3,...] )----
  _xnodes=new double[ _NoElNodes];
  // shape functions at nodes (_phi_map_nodes[sh1(x1),sh1(x2),....  sh2(x1),sh2(x2),.... ])
  _phi_map_nodes = new double[_NoShape* _NoElNodes];
  // 1st derivatives at nodes x _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map_nodes = new double[ _NoShape * _NoElNodes];
  // element dofs  (reference to nodes) --------------------------------------------------
  _fem_conn=new int[_NoShape];
  // gaussian (g) coordinates (_xg[xg1,xg2,xg3,...yg1,yg2,yg3,....zg1,zg2,zg3,...] )-----
  _xg=new double[_NoGauss];
  // weights [wg1,wg2,wg3,.......]
  _weight = new double[_NoGauss];
  // shape functions phi[sh1(xg1),sh1(xg2),....  sh2(xg1),sh2(xg2),.... ]
  _phi_map = new double[ _NoGauss*_NoShape ];
  // 1st derivatives at g _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map = new double[ _NoShape* _NoGauss];
  // 2nd derivatives at g _ddphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxx_map = new double[_NoShape * _NoGauss];
  //   -------------------------------------------------------------------------------------
 // dofs -----------------------------------------------------------------------
 for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
 double *xnodes=new double [_NoElNodes*_dim];
 double *point=new double[_dim];
    //
    // ---------------------------------------------------------------------------------------
    // nodes ----------------------------------------------------------------------------
  get_refspace_nodes(geom_elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
    // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
    for(int i = 0; i < _NoElNodes; i++) {  // points
      for(int dir = 0; dir < _dim; dir++) { 
        _xnodes[i+_NoElNodes*dir]=xnodes[dir+_dim*i]; 
        point[dir] = _xnodes[i+_NoElNodes*dir];
      }

    for(int j = 0; j < _NoShape; j++) {  // shapse
      //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
      _phi_map_nodes[j*_NoElNodes+i]=fe_lagrange_1D_shape( elem_order, j, point[0]);//Rec_Quad_Phi(j, point,_dim);

      // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = fe_lagrange_1D_shape_deriv( elem_order, j,dir, point[0]);// Rec_Quad_DPhi(j,point,_dim,dir);     // dphi/dxi
      }
    }
  }
   delete[]xnodes;
// gaussian points --------------------------------------------------------->
   
 for(int i = 0; i < _NoGauss; i++) {
    _weight[i] =_qrule->w(i); // weight
     _qrule->qp(point,i);     // ith-gaussian point ->point(i)
    for(int j = 0; j < _NoShape; j++) {     // loop over test function id
      _phi_map[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
      
       
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv( elem_order, j,dir, point[0]);
        for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
          _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
      }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
  } 
  delete []point;
 
//   if(_FamType==0) {  // 0= Lagrangian  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     // lagrangian -> n_point= n_shape =_order DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
//     //            -> pol degree=_order/2 (int div)
//     _deg=0;// DDDDDDDDDDDDDDDDDDDDDDDDD
// 
//     switch(_order) {  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
//     case 0:
//       init_pie();  break;
//     case 1:
//       init_lin();  break;
//     case 2:
//       init_qua();  break;
//     default:
//       std::cout<<"MGFE:"<<" Dim "<<_dim<<" case "<<_order<<" Not implemented";  exit(3);
//     }// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// 
//   }//  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
// 
//   if(_FamType==1) {  //1= Raviart FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     if(_GeomType==1) {  //HEX
//       _deg=1;   if(_dim==1)  _deg=0;
//       switch(_order) {
//       case 0:
//         init_pie();   break;
//       case 1:
//         init_lin_rec_rt();   break;
//       default:
//         std::cout<<"MGFE::MGFE:"<<" Dim "<<_dim<<" case "<<_order<<"RT Not implemented"; exit(3);
//       }
//     }
//     if(_GeomType==0) {  // TETRA GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//       std::cout<< " Not implemented";
//       exit(0);
//     }//  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//   }//  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//   
//   if(_FamType==2) {  //2= discontinuous FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     std::cout<< " Not implemented";
//     exit(0);
//   } //  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#ifdef PRINT_INFO
  std::cout << "\n MGFE: elem_name(Dim): "           <<elem_name<<"("<< _dim <<")" ;
  std::cout << "\n       gaussian points: "          << _NoGauss ;
  std::cout << "\n       number of shape functions " << _NoShape ;
  std::cout << "\n       polynomial order "          << elem_order;
  std::cout << " \n";
#endif
  return;
}



// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic Raviart Thomas non continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// //===============================================================================================
// template <>
// std::array<double,3 >  MGFE<1,LAGRANGE>::Rec_RT_Phi(int iface, double x[],int dim) {
// //===============================================================================================
//   switch(iface) {
//   case 0:
//     return {1.,0.,0.};    break;
//   case 1:
//     return {1.,0.,0.};    break;
//   default:
//     std::cout<< "MGFE1::RT_quad Error iface "; return  {0.,0.,0.};  exit(3);
//   }
//   return  {0.,0.,0.};
// }
// //===============================================================================================
// 
// 
// template <>
// std::array<double,3>  MGFE<1,LAGRANGE>::Rec_RT_DPhi(int /*iface*/, double x[],int dim, int /*dir*/) {
//   return  {0.,0.,0.};  
// }



// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// Basic special operators
// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO



template <> double MGFE<1,LAGRANGE>::ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[]) {
  double det = Matrix[0] ;
  double idet = 1. / det;
  InvMatrix[0] =  idet;   // dxi dx
  return det;
}
template <> double MGFE<1,LAGRANGE>::ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[]) {
  double det = 1.;
  InvMatrix[0] = 1;   // dxi dx
  return det;
}


//===============================================================================================
template <> double MGFE<1,LAGRANGE>::ComputeInverseMatrix(
  double Matrix[], 
  double InvMatrix[], 
  int dim) {
  double det = 0.;

  if(dim == 1)   det = ComputeInverseMatrix_vol(Matrix,InvMatrix);
  if(dim == 0) det = ComputeInverseMatrix_sur(Matrix,InvMatrix);
  
  return det;
}
//===============================================================================================







//===============================================================================================
template <>  double MGFE<1,LAGRANGE>::Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim) {
  double sum=0.; /*int ndof=8;*/

  for(int k=0; k<ndof; k++) {
    double val=u[k+j];

    for(int kdim=0; kdim<idim; kdim++) {
      val *= (1.+ _CooH27[k+kdim*27]*xi[kdim]);
    }

    sum += val;
  }

  return sum/ndof;
}
//===============================================================================================
template <>  double MGFE<1,LAGRANGE>::dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim) {
  double sum=0.; /*int ndof=8;*/

  for(int k=0; k<ndof; k++) {
    double val=u[k+j];

    for(int kdim=0; kdim<i; kdim++) {
      val *= (1.+_CooH27[k+kdim*27]*xi[kdim]);
    }

    val *= _CooH27[k+i*27];

    for(int kdim=i+1; kdim<idim; kdim++) {
      val *= (1.+_CooH27[k+kdim*27]*xi[kdim]);
    }

    sum += val;
  }

  return sum/ndof;
}
//========================================================================================
// =========================================================
//  GET phi dphi ddphi
//==========================================================
//========================================================================================
template <> void MGFE<1,LAGRANGE>::get_dphi_on_given_node(
  const int dim, 
  double ElemCoords[], 
  double CanPos[], 
  double dphi[]
) {
  // double* InvJac = new double[dim * dim];
  // double* gradphi_g = new double[dim];
  // int nShape = _NoShape;
  // JacOnGivenCanCoords(dim, ElemCoords, CanPos, InvJac, _GeomType, nShape);
  // 
  // for(int eln = 0; eln < nShape; eln++) {
  //   for(int idim = 0; idim < dim; idim++) {
  //     gradphi_g[idim] = FirstDerivateOfLocalPhi(eln, CanPos, dim, idim, _GeomType);
  //   }
  // 
  //   for(int idim = 0; idim < dim; idim++) {
  //     double sum = 0.;
  // 
  //     for(int jdim = 0; jdim < dim; jdim++) { sum += InvJac[jdim + idim * dim] * gradphi_g[jdim]; }
  // 
  //     dphi[eln + idim * nShape] = sum;
  //   }
  // }
  // 
  // delete[] gradphi_g;
  // delete[] InvJac;
  return;
}//======================================================================
template <> double   MGFE<1,LAGRANGE>::compute_fe_face_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz_bd[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
 ){
  
   for(int l=0;l<_NoShape;l++)       {
    phi[l]=_phi_map[l* _NoGauss  +gp] ;// shape
   dphidxi[l]=_dphidxez_map[l*_NoGauss + gp];  phi[l]=dphidxi[l];
}
 
   // int NDOF_FEMB=3;
 
  // Values to compute at gaussian points
  double dxdxi = 0.;    double dydxi = 0.;  // d(x,y)d(xi,eta)
  // const int nshape = _NoShape;

  //   int offset=_NoShape[0]*_NoGauss[0];
  for(int s = 0; s < _NoShape; s++) {
    // int sng = s * _NoGauss + gp;
    // double dphidxi = _dphidxez_map[s * _NoGauss + gp];
    dxdxi += xyz_bd[s] * dphidxi[s];
    dydxi += xyz_bd[s +  _NoShape] * dphidxi[s];
  }

  // surface weighted jacobean
  double det = sqrt(dxdxi * dxdxi + dydxi * dydxi);
  double InvJac0 = 1. / det;
  
  for(int l=0;l<_NoShape;l++)  {  dphidx[l]=InvJac0* dphidxi[l];  }
  
  
  return  det;
}

template <> double   MGFE<1,LAGRANGE>::compute_fe_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
                                                      double InvJac_in[]          
 ){
     // Compute the shape function values (and derivatives)
//   // at the Quadrature points qp. 
   
  //  ElemType  elem_name=_fe_type.fetype; 
  //  ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
  //  Order elem_order= _fe_type.order ;
  //   _order= (int)elem_order;
  // // number of nodes (_NoElNodes) in the reference EDGE3 element
  // _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  // _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
  
  for(int l=0;l<_NoShape;l++)       {
    phi[l]=_phi_map[l* _NoGauss  +gp] ;// shape
   dphidxi[l]=_dphidxez_map[(l)*_NoGauss + gp]; 
}
  double InvJac=InvJac_in[0];
  double x_xi = 0.; 
  for(int s = 0; s < _NoShape; s++) {
    x_xi += xyz[s] *  _dphidxez_map[s * _NoGauss + gp];
  }
  double det = x_xi; InvJac =1./det;
  for(int l=0;l<_NoShape;l++)  {  dphidx[l]=InvJac* dphidxi[l];  dphi[l]=dphidx[l]; }
  return det;
   
   
 }

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<1,LAGRANGE>::get_phi_gl_g(
const int qp,    // gaussian point <-
double phi[]     // shape functions ->
) {                  // =================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<1,LAGRANGE>::get_phi_gl_g(
const int qp,             // gaussian point   <-
std::vector<double>& phi  // shape functions  ->
) { // ==============================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];

  return;
}

// // ==================================================================
// void MGFE1::get_phi_gl_g(
//   double point[],           // gaussian point <
//   double phi[],
//   int fam_type,
//   int sdim
// ) { // ==============================================================
// 
//   const int el_nnodes = _NoShape;   // # of shape functions
// 
//   for(int eln = 0; eln < el_nnodes; eln++) {
// 
//     phi[eln] = LocalPhi(eln, point, sdim, fam_type) ;
// 
//   }
// 
//   return;
// }
// void MGFE1::get_phi_g_arb_el(
// const int kdim,    // dimension <-
// const int qp,      // gaussian point <-
// double phi[],      // shape functions ->
// int FamilyType) {  // =================================================
//   for(int ish = 0; ish < _GNoShape[FamilyType][kdim - 1]; ish++) {
//     phi[ish] = _Gphi_map1[(kdim - 1) + FamilyType * 3][ish * _GNoGauss1[FamilyType][kdim - 1] + qp];
//   }
//   return;
// }
// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ]
template <> void MGFE<1,LAGRANGE>::get_dphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[]           // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  std::vector<double> gradphi_g(_dim);  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < _dim; idim++) gradphi_g[idim]=_dphidxez_map[lqp+idim*goffset];

    for(int idim = 0; idim < _dim; idim++) {
      double sum=0.; for(int jdim=0;jdim<_dim;jdim++) sum +=InvJac[jdim+idim*_dim]*gradphi_g[jdim];
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// void MGFE1::get_dphi_g_arb_el(
// const int kdim,         // dimension <-
// const int qp,           // gaussian point <
// const double InvJac[],  // Jacobean
// double dphi[],          // global derivatives ->
// int FamilyType) {       // =========================================
//
//   const int el_nnodes = _GNoShape[FamilyType][kdim - 1];             // # of shape functions
//   const int el_ngauss = _GNoGauss1[FamilyType][kdim - 1];            // # of gauss points
//   const int goffset = el_nnodes * _GNoGauss1[FamilyType][kdim - 1];  // gauss offset
//
//   double gradphi_g[DIMENSION];  // temp grad phi
//
//   for(int eln = 0; eln < el_nnodes; eln++) {
//     int lqp = eln * el_ngauss + qp;
//     for(int idim = 0; idim < _dim; idim++) {
//       gradphi_g[idim] = _Gdphidxez_map1[kdim - 1 + 3 * FamilyType][lqp + idim * goffset];
//     }
//     for(int idim = 0; idim < _dim; idim++) {
//       double sum = 0.;
//       for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }
//       dphi[eln + idim * el_nnodes] = sum;
//     }
//   }
//
//   return;
// }
// =============================================================
/// Shape functions 2nd derivatives at the gauss point qp
/// ddphi[id*_dim*_dim+i*_dim+j]= Tij
///    Txx Txy Txz
///    Tyx Tyy Tyz
///    Tzx Tzy Tzz
///
///     ddphi[Txx(0),Txy(0),Txz(0),Tyx(0),Tyy(0),Tyz(0),Tzx(0),Tzy(0),Tzz(0),
///                 ....................................
///                 Txx(id),Txy(id),Txz,Tyx(id),Tyy,Tyz,Tzx(id),Tzy(id),Tzz(id)]
///   tensor  order  at the gauss point qp, node inode
template <> void MGFE<1,LAGRANGE>::get_ddphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double ddphi[]          // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  // double ddphi_loc[NDOF_FEM * DIMENSION * DIMENSION];
   double ddphi_loc[3 * 1 * 1];

  for(int eln = 0; eln < el_nnodes; eln++) {   // LOOP OVER NODES ===================================
    const int shift = eln * _dim * _dim;       // offset for node derivatives
    const int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim * _dim; idim++) {
      ddphi_loc[shift + idim] = _dphidxx_map[lqp + idim * goffset];
    }

    // double Hess_tmp[DIMENSION * DIMENSION], JacTf[DIMENSION * DIMENSION];
      double Hess_tmp[1 * 1], JacTf[1 * 1];

    for(int init = 0; init < _dim; init++)
      for(int jnit = 0; jnit < _dim; jnit++) {
        Hess_tmp[init * _dim + jnit] = 0.;
        ddphi[shift + init * _dim + jnit] = 0.;
        JacTf[init * _dim + jnit] = InvJac[jnit * _dim + init];
      }

    for(int ii = 0; ii < _dim; ii++)
      for(int jj = 0; jj < _dim; jj++)
        for(int ss = 0; ss < _dim; ss++) {
          Hess_tmp[ii * _dim + jj] += ddphi_loc[shift + ii * _dim + ss] * JacTf[ss * _dim + jj];
        }

    // calculting d2 phi / dxi dxj
    for(int ii = 0; ii < _dim; ii++)
      for(int jj = 0; jj < _dim; jj++)
        for(int ss = 0; ss < _dim; ss++) {
          ddphi[shift + ii * _dim + jj] += InvJac[ii * _dim + ss] * Hess_tmp[ss * _dim + jj];
        }
  }  // END LOOP OVER NODES =====================================================================

  // Pay attention, different order with respect to first order derivatives!
  //   delete[]ddphi_loc;
  return;
}
// =============================================================
/// Shape functions derivatives at the nodal points
template <> void MGFE<1,LAGRANGE>::get_dphi_node(
const int /*kdim*/,         // dimension <-
const int node,         // nodal point <-
const double InvJac[],  // Jacobean <-
double dphi[]           // global derivatives at the nodal point ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  // double gradphi_g[DIMENSION];  // temp grad phi
  double gradphi_g[1];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + node;

    for(int idim = 0; idim < _dim; idim++) {
      gradphi_g[idim] = _dphidxez_map_nodes[lqp + idim * goffset];
    }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }

      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// // =============================================================
// /// Shape functions derivatives at the nodal points
// void MGFE1::get_dphi_arb_node(
// std::vector<double> NodeCoord,
// const int order,
// double InvJac[],
// double dphi[]) {  // =========================================
//   const int kdim = NodeCoord.size();
//
//   const int el_nnodes = _NoShape;             // # of shape functions
//   const int el_ngauss = _NoGauss1[kdim - 1];            // # of gauss points
//   const int goffset = el_nnodes * _NoGauss1[kdim - 1];  // gauss offset
//
//   double gradphi_g[DIMENSION];  // temp grad phi
//
//   int el_nodes = (order == 1) ? 4 : 9;
//   if(kdim == 3) { el_nodes *= (order == 1) ? 2 : 3; }
//
//   for(int eln = 0; eln < el_nodes; eln++) {
//     gradphi_g[0] = 1.;
//     gradphi_g[1] = 1.;
//     gradphi_g[kdim - 1] = 1.;
//     for(int idim = 0; idim < kdim; idim++) {   // derivative of eln-th phi ind idim direction
//       const int idim2 = (idim + 1) % kdim;
//       gradphi_g[idim] =
//       0.5 * _CooQ9[eln + idim * _Q9Off] * (0.5 * (1. + NodeCoord[idim2] * _CooQ9[eln + idim2 * _Q9Off]));
//     }
//     for(int idim = 0; idim < _dim; idim++) {
//       double sum = 0.;
//       for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }
//       dphi[eln + idim * el_nodes] = sum;
//     }
//   }
//   return;
// }

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<1,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,         // dimension <-
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[],          // global derivatives ->
int sdim) {             // =========================================

  const int el_nnodes = _NoShape;   // # of shape functions
  const int el_ngauss = _NoGauss;  // # of guass points
  const int goffset = el_nnodes * el_ngauss;  // gauss offset

  double gradphi_g[3];  // temp grad phi

  for(int idim = 0; idim < 3; idim++) { gradphi_g[idim] = 0.; }
  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < sdim; idim++) { gradphi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }
    for(int idim = 0; idim < sdim; idim++) {
      double sum = 0.;
      for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// // =============================================================
// /// Shape functions derivatives at the gauss point qp
// template <> void MGFE<1,LAGRANGE>::get_dphi_gl_g(
// const int /*kdim*/,         // dimension <-
// double point[],           // gaussian point <
// const double InvJac[],  // Jacobean
// double dphi[],          // global derivatives ->
// int sdim,
// int fam_type) {             // =========================================
// 
//   const int el_nnodes = _NoShape;   // # of shape functions
//   const int el_ngauss = _NoGauss;  // # of guass points
//   const int goffset = el_nnodes * el_ngauss;  // gauss offset
// 
//   double gradphi_g[3];  // temp grad phi
// 
//   for(int idim = 0; idim < 3; idim++) { gradphi_g[idim] = 0.; }
// 
//   for(int eln = 0; eln < el_nnodes; eln++) {
// 
//     for(int idim = 0; idim < sdim; idim++) { gradphi_g[idim] = FirstDerivateOfLocalPhi(eln, point, sdim, idim, fam_type); }
// 
//     for(int idim = 0; idim < sdim; idim++) {
//       double sum = 0.;
// 
//       for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }
// 
//       dphi[eln + idim * el_nnodes] = sum;
//     }
//   }
// 
//   return;
// }

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<1,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,            // dimension <-
const int qp,              // gaussian point <-
const double InvJac[],     // Jacobean
std::vector<double>& dphi  // global derivatives ->
) {                            // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of guass points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double dphidxi_g[1];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim; idim++) { dphidxi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * dphidxi_g[jdim]; }

      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

template <>  bool MGFE<1,LAGRANGE>::phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it, double tol) {return false;}


